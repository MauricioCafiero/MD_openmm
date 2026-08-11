#!/usr/bin/env python
"""Well-tempered metadynamics on the rotaxane shuttle coordinate via PLUMED.

Loads the solvated rotaxane system built by ``omd build-multimol``
(outputs/system.xml + outputs/complex.pdb), attaches the OpenMM PLUMED plugin
(openmm-plumed) with the generated plumed.dat, and runs:

  phase 1 -- unbiased minimize + short NPT equilibration (no PLUMED bias), so the
             box / water settle without the metad hills seeing the relaxation;
  phase 2 -- WT-MetaD production with the PLUMED bias on, started from the
             equilibrated state. HILLS + COLVAR land in the out-dir.

PLUMED uses 1-based atom indices matching the complex.pdb atom order (ROD then
WHL then waters/ions); make_plumed.py writes those from the topology.

Requires: ``conda install -c conda-forge openmm-plumed`` into the openmm-md env
(not done automatically -- install before running).
"""
from __future__ import annotations

import argparse
from pathlib import Path

from openmm import app, openmm as mm, unit, XmlSerializer

from openmm_md.config import Config
from openmm_md.dynamics import get_platform, _new_sim, _build_with_fallback


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
    ap.add_argument("--steps", type=int, default=5_000_000,
                    help="biased production steps (default 5M = 10 ns at 2 fs)")
    ap.add_argument("--equil-steps", type=int, default=5000,
                    help="unbiased NPT equilibration steps before the bias turns on")
    ap.add_argument("--platform", default="auto", help="OpenCL | CPU | auto")
    ap.add_argument("--report-interval", type=int, default=500)
    ap.add_argument("--traj-interval", type=int, default=1000)
    args = ap.parse_args()

    cfg = Config()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pdb = app.PDBFile(str(args.topology))
    topology, positions = pdb.topology, pdb.positions

    # ---- phase 1: unbiased minimize + NPT equilibration (no PlumedForce) ----
    equi_sys = load_system(args.system)
    add_barostat(equi_sys, cfg)
    plat = get_platform(args.platform)
    equi_sim, equi_name = _build_with_fallback(topology, equi_sys, cfg, plat.getName())
    equi_sim.context.setPositions(positions)
    print(f"[metad] minimizing + {args.equil_steps}-step unbiased equil on {equi_name} ...")
    equi_sim.minimizeEnergy(tolerance=cfg.minimize_tolerance)
    equi_sim.context.setVelocitiesToTemperature(cfg.temperature * unit.kelvin)
    equi_sim.step(args.equil_steps)
    st = equi_sim.context.getState(getPositions=True, getVelocities=True,
                                    getEnergy=True)
    print(f"[metad] equil KE={st.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole):.1f} kJ/mol, "
          f"PE={st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.1f} kJ/mol")

    # ---- phase 2: biased production with the PLUMED WT-MetaD force ----
    from openmmplumed import PlumedForce  # noqa: import here so the equil runs even if the plugin is missing
    plumed_script = Path(args.plumed).read_text()
    # PLUMED writes HILLS / COLVAR relative to the cwd, not the out-dir -> point
    # them at absolute paths inside the out-dir so a run from anywhere lands them
    # alongside the trajectory.
    out_abs = args.out_dir.resolve()
    plumed_script = plumed_script.replace("FILE=HILLS", f"FILE={out_abs}/HILLS")
    plumed_script = plumed_script.replace("FILE=COLVAR", f"FILE={out_abs}/COLVAR")
    prod_sys = load_system(args.system)
    add_barostat(prod_sys, cfg)
    prod_sys.addForce(PlumedForce(plumed_script))
    print(f"[metad] attached PlumedForce from {args.plumed}")

    prod_sim = _new_sim(topology, prod_sys, cfg, equi_name)
    prod_sim.context.setPositions(st.getPositions())
    prod_sim.context.setVelocities(st.getVelocities())
    prod_sim.context.setPeriodicBoxVectors(*st.getPeriodicBoxVectors())

    prod_sim.reporters.append(
        app.StateDataReporter(
            str(args.out_dir / "energy.csv"), args.report_interval,
            step=True, time=True, potentialEnergy=True, kineticEnergy=True,
            totalEnergy=True, temperature=True, volume=True, separator=",",
        )
    )
    prod_sim.reporters.append(
        app.DCDReporter(str(args.out_dir / "traj.dcd"), args.traj_interval))

    print(f"[metad] WT-METAD production {args.steps} steps on {equi_name} "
          f"(HILLS + COLVAR written to {args.out_dir}) ...")
    prod_sim.step(args.steps)

    final = prod_sim.context.getState(getPositions=True)
    with open(args.out_dir / "final.pdb", "w") as f:
        app.PDBFile.writeFile(topology, final.getPositions(), f, keepIds=True)
    print(f"[metad] done. traj: {args.out_dir/'traj.dcd'}  "
          f"CV: {args.out_dir/'COLVAR'}  hills: {args.out_dir/'HILLS'}")


if __name__ == "__main__":
    main()