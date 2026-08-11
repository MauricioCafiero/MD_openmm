#!/usr/bin/env python
"""Resume the WT-METAD production after the original run died at step ~11.944M
(47.8% of 25M), with no OpenMM checkpoint saved.

Warm restart: take positions + the NPT box from the last traj.dcd frame, re-init
velocities to T (the Langevin thermostat re-thermalizes within ps), and put PLUMED
into RESTART mode so it reloads the existing ~23,892 HILLS as the bias and APPENDS
new ones (HILLS/COLVAR are not truncated). Runs the remaining steps into NEW
reporter files (traj_resume.dcd, energy_resume.csv) so the originals stay intact;
also writes checkpoint_resume.bin periodically so a future death is recoverable by
exact reload instead of another warm restart.

HILLS/COLVAR ARE appended (PLUMED RESTART), so sum_hills at the end uses the full
bias history. The OpenMM step counter restarts at 0 in this context, so the
appended energy/traject have their own step numbering -- concatenate for analysis.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import mdtraj as md
from openmm import app, openmm as mm, unit, XmlSerializer
from openmmplumed import PlumedForce

from openmm_md.config import Config
from openmm_md.dynamics import get_platform, _build_with_fallback

TOTAL_STEPS = 25_000_000
TRAJ_INTERVAL = 1000      # matches the original run's DCDReporter stride
REPORT_INTERVAL = 500     # matches the original run's StateDataReporter stride
CHECKPOINT_INTERVAL = 50_000  # ~every 100 ps; cheap, single overwrite file


def load_system(system_xml: Path) -> mm.System:
    with open(system_xml) as f:
        return XmlSerializer.deserializeSystem(f.read())


def add_barostat(system: mm.System, cfg: Config):
    if cfg.pressure > 0:
        system.addForce(
            mm.MonteCarloBarostat(
                cfg.pressure * unit.atmosphere,
                cfg.temperature * unit.kelvin,
                cfg.barostat_interval,
            )
        )


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--system", type=Path, default=Path("../outputs/system.xml"))
    ap.add_argument("--topology", type=Path, default=Path("../outputs/complex.pdb"))
    ap.add_argument("--plumed", type=Path, default=Path("plumed.dat"))
    ap.add_argument("--out-dir", type=Path, default=Path("metad_run"))
    ap.add_argument("--platform", default="auto")
    ap.add_argument("--dry-run", action="store_true",
                    help="build + step(1) to validate setup (RESTART HILLS reload, "
                         "state, one integrator step), then exit without running. "
                         "No reporter output or hill deposition (strides >1).")
    args = ap.parse_args()

    cfg = Config()
    out = args.out_dir
    pdb = app.PDBFile(str(args.topology))
    topology = pdb.topology

    # ---- warm restart: last DCD frame -> positions + NPT box ----
    traj = md.load_dcd(str(out / "traj.dcd"), top=str(args.topology))
    last_step = traj.n_frames * TRAJ_INTERVAL
    remaining = TOTAL_STEPS - last_step
    pos_nm = traj[-1].xyz[0]
    box_nm = traj[-1].unitcell_vectors[0]
    positions = pos_nm * unit.nanometer
    box_vectors = [v * unit.nanometer for v in box_nm]
    print(f"[resume] warm restart from traj.dcd frame {traj.n_frames} "
          f"(step ~{last_step}/{TOTAL_STEPS}); remaining = {remaining} steps", flush=True)

    # ---- production system + PLUMED in RESTART mode ----
    # RESTART makes PLUMED reload the existing HILLS as the bias and APPEND new
    # hills / COLVAR rows instead of truncating. Inject it at the top of the script
    # (don't modify the on-disk plumed.dat). Keep HILLS/COLVAR on absolute paths so
    # they land in the out-dir regardless of cwd.
    plumed_script = "RESTART\n" + args.plumed.read_text()
    out_abs = out.resolve()
    plumed_script = plumed_script.replace("FILE=HILLS", f"FILE={out_abs}/HILLS")
    plumed_script = plumed_script.replace("FILE=COLVAR", f"FILE={out_abs}/COLVAR")

    prod_sys = load_system(args.system)
    add_barostat(prod_sys, cfg)
    prod_sys.addForce(PlumedForce(plumed_script))
    print(f"[resume] attached PlumedForce (RESTART) from {args.plumed}", flush=True)

    plat = get_platform(args.platform)
    prod_sim, name = _build_with_fallback(topology, prod_sys, cfg, plat.getName())
    prod_sim.context.setPositions(positions)
    prod_sim.context.setPeriodicBoxVectors(*box_vectors)
    prod_sim.context.setVelocitiesToTemperature(cfg.temperature * unit.kelvin)

    # NEW reporter files (OpenMM reporters don't append; originals stay intact).
    prod_sim.reporters.append(
        app.StateDataReporter(
            str(out / "energy_resume.csv"), REPORT_INTERVAL,
            step=True, time=True, potentialEnergy=True, kineticEnergy=True,
            totalEnergy=True, temperature=True, volume=True, separator=",",
        )
    )
    prod_sim.reporters.append(
        app.DCDReporter(str(out / "traj_resume.dcd"), TRAJ_INTERVAL))
    # single overwrite file, always the latest full Context state -> a future death
    # is recoverable by Simulation.loadCheckpoint instead of a warm restart
    prod_sim.reporters.append(
        app.CheckpointReporter(str(out / "checkpoint_resume.bin"), CHECKPOINT_INTERVAL))

    print(f"[resume] WT-METAD resume {remaining} steps on {name} "
          f"(HILLS/COLVAR appended; traj->traj_resume.dcd, energy->energy_resume.csv)", flush=True)
    if args.dry_run:
        prod_sim.step(1)  # triggers PLUMED RESTART HILLS reload + one step; writes nothing
        print("[resume] dry-run OK (setup + RESTART + 1 step); not running production.", flush=True)
        return
    prod_sim.step(remaining)

    final = prod_sim.context.getState(getPositions=True)
    with open(out / "final.pdb", "w") as f:
        app.PDBFile.writeFile(topology, final.getPositions(), f, keepIds=True)
    print(f"[resume] done. Reconstruct F(d) with: "
          f"plumed sum_hills --hills {out}/HILLS --mintozero --outfile {out}/Fes.dat", flush=True)


if __name__ == "__main__":
    main()