"""Run MD: minimize -> equilibrate (with restraints) -> production.

Platform probing: OpenMM has no MPS/Metal platform. On Apple Silicon the pip/conda
builds expose Reference, CPU, and OpenCL. get_platform("auto") tries OpenCL (GPU)
first, then CPU, then Reference, verifying each can actually create a context.

Hybrid platform (this Mac): Apple's OpenCL GPU device is single-precision only, and
single precision cannot stably run the *restrained* minimize+equilibration on a
large solvated system -- L-BFGS gets stuck at a catastrophic energy and the first
restrained dynamics steps NaN. CPU (double) handles it fine, and OpenCL single is
stable for production once the restraints are released (k=0). So when restraints
are on and the production platform is single-precision OpenCL, we equilibrate on
CPU and produce on OpenCL, transferring the equilibrated state across.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from openmm import app, openmm as mm, unit, XmlSerializer

from .config import Config

# Standard amino acids + common AMBER protonation variants + terminal caps.
PROTEIN_RES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
    "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CYX", "CYM", "ASH", "GLH", "LYN",
    "ACE", "NME", "NH2", "FOR",
}


def protein_heavy_indices(topology) -> list[int]:
    idx = []
    for atom in topology.atoms():
        if atom.residue.name in PROTEIN_RES and atom.element is not None and atom.element.symbol != "H":
            idx.append(atom.index)
    return idx


def get_platform(requested: str = "auto") -> mm.Platform:
    """Pick a platform, verifying it can actually create a context.

    OpenMM has no MPS/Metal platform; on Apple Silicon the conda build exposes
    Reference, CPU, and OpenCL. We probe *bare* (no Precision default set --
    setting an unsupported Precision default, e.g. mixed/double, poisons the
    OpenCL singleton on this Mac where the GPU device lacks double precision).
    Precision is applied per-context when the Simulation is built.
    """
    if requested != "auto":
        p = _try_platform(requested)
        if p is None:
            raise RuntimeError(f"requested platform {requested!r} is not available")
        return p
    for name in ("OpenCL", "CPU", "Reference"):
        p = _try_platform(name)
        if p is not None:
            return p
    raise RuntimeError("no OpenMM platform could create a context")


def _try_platform(name):
    try:
        p = mm.Platform.getPlatformByName(name)
        s = mm.System()
        s.addParticle(1.0)
        integ = mm.VerletIntegrator(0.001)
        ctx = mm.Context(s, integ, p)
        del ctx
        return p
    except Exception:
        return None


def _new_sim(topology, system, cfg: Config, platform_name: str):
    """Build a Simulation on ``platform_name`` (Precision=single for OpenCL)."""
    integ = mm.LangevinMiddleIntegrator(
        cfg.temperature * unit.kelvin,
        cfg.friction / unit.picosecond,
        cfg.timestep * unit.femtosecond,
    )
    integ.setConstraintTolerance(1e-6)
    platform = mm.Platform.getPlatformByName(platform_name)
    props = {"Precision": cfg.opencl_precision} if platform_name == "OpenCL" else None
    return app.Simulation(topology, system, integ, platform, platformProperties=props)


def _build_with_fallback(topology, system, cfg: Config, platform_name: str):
    """Build a Simulation, falling back OpenCL->CPU if it fails on the real system.

    The 1-particle platform probe can pass while a real PME system fails (Apple's
    deprecated OpenCL lacks features/precision some systems need). Returns
    (simulation, actual_platform_name).
    """
    try:
        return _new_sim(topology, system, cfg, platform_name), platform_name
    except Exception as e:
        if platform_name == "OpenCL" and cfg.platform == "auto":
            print(f"[run] OpenCL failed on the real system ({str(e)[:70]}); falling back to CPU")
            return _new_sim(topology, system, cfg, "CPU"), "CPU"
        raise


def run(
    system_xml: Path,
    topology_pdb: Path,
    out_dir: Path,
    cfg: Config | None = None,
    steps: int | None = None,
):
    """minimize -> equilibrate (restraints on) -> release restraints -> production.

    With restraints active on single-precision OpenCL, equilibrate on CPU and
    produce on OpenCL (state transferred across); otherwise one platform throughout.
    """
    cfg = cfg or Config()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    steps = steps if steps is not None else cfg.production_steps

    with open(system_xml) as f:
        system = XmlSerializer.deserializeSystem(f.read())
    pdb = app.PDBFile(str(topology_pdb))
    topology, positions = pdb.topology, pdb.positions

    # NPT barostat (added before context creation)
    if cfg.pressure > 0:
        system.addForce(
            mm.MonteCarloBarostat(
                cfg.pressure * unit.atmosphere,
                cfg.temperature * unit.kelvin,
                cfg.barostat_interval,
            )
        )

    # Positional restraints on protein heavy atoms. k is a *global* parameter,
    # so we can zero it out for production without rebuilding the system.
    # A protein-free system (single-molecule solvation) has no heavy atoms to
    # restrain -- in that case restraints are a no-op and we must NOT take the
    # hybrid CPU-equil path (it exists only because restrained equilibration on
    # single-precision OpenCL is unstable; with nothing restrained there is no
    # such issue, so we run one platform throughout).
    restraint = None
    if cfg.restrain_protein:
        heavy = protein_heavy_indices(topology)
        if heavy:
            restraint = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
            restraint.addGlobalParameter("k", cfg.restraint_weight * unit.kilojoule_per_mole / unit.nanometer**2)
            restraint.addPerParticleParameter("x0")
            restraint.addPerParticleParameter("y0")
            restraint.addPerParticleParameter("z0")
            for i in heavy:
                p = positions[i].value_in_unit(unit.nanometer)
                restraint.addParticle(i, [float(p[0]), float(p[1]), float(p[2])])
            system.addForce(restraint)
            print(f"[run] restraining {len(heavy)} protein heavy atoms")
        else:
            print("[run] no protein heavy atoms found; restraints disabled")

    # ---- choose platforms ----
    prod_platform = get_platform(cfg.platform)
    opencl_single = prod_platform.getName() == "OpenCL" and cfg.opencl_precision == "single"
    cpu_equil = restraint is not None and opencl_single
    prod_name = prod_platform.getName()
    requested_equi = "CPU" if cpu_equil else prod_name

    # ---- Phase 1: minimize + equilibrate (restraints on) ----
    equi_sim, equi_name = _build_with_fallback(topology, system, cfg, requested_equi)
    # if the *requested* equilibration platform fell back to CPU (OpenCL failed on
    # the real PME system), produce on CPU too so we reuse the one good context.
    if equi_name != requested_equi:
        prod_name = "CPU"
    if cpu_equil:
        print(f"[run] OpenCL(single) can't run restrained equilibration; "
              f"equilibrating on CPU, producing on {prod_name}")
    else:
        print(f"[run] platform: {equi_name}"
              + (f" (precision={cfg.opencl_precision})" if equi_name == "OpenCL" else ""))

    equi_sim.context.setPositions(positions)
    print("[run] minimizing ...")
    equi_sim.minimizeEnergy(maxIterations=cfg.minimize_max_iter, tolerance=cfg.minimize_tolerance)
    equi_sim.context.setVelocitiesToTemperature(cfg.temperature * unit.kelvin)
    print(f"[run] equilibrating {cfg.equilibrate_steps} steps (restraints on) on {equi_name} ...")
    equi_sim.step(cfg.equilibrate_steps)

    # ---- Phase 2: production (restraints off) ----
    if prod_name == equi_name:
        prod_sim = equi_sim
    else:
        try:
            prod_sim = _new_sim(topology, system, cfg, prod_name)
            st = equi_sim.context.getState(getPositions=True, getVelocities=True)
            prod_sim.context.setPositions(st.getPositions())
            prod_sim.context.setVelocities(st.getVelocities())
            prod_sim.context.setPeriodicBoxVectors(*st.getPeriodicBoxVectors())
        except Exception as e:
            print(f"[run] {prod_name} failed for production ({str(e)[:70]}); using CPU")
            prod_sim = equi_sim
            prod_name = "CPU"

    if restraint is not None:
        prod_sim.context.setParameter("k", 0.0)

    prod_sim.reporters.append(
        app.StateDataReporter(
            str(out_dir / "energy.csv"), cfg.report_interval,
            step=True, time=True, potentialEnergy=True, kineticEnergy=True,
            totalEnergy=True, temperature=True, volume=True, density=True,
            separator=",",
        )
    )
    prod_sim.reporters.append(app.DCDReporter(str(out_dir / "traj.dcd"), cfg.traj_interval))
    prod_sim.reporters.append(app.CheckpointReporter(str(out_dir / "checkpoint.chk"), cfg.checkpoint_interval))

    print(f"[run] production {steps} steps (restraints off) on {prod_name} ...")
    prod_sim.step(steps)

    state = prod_sim.context.getState(getPositions=True)
    with open(out_dir / "final.pdb", "w") as f:
        app.PDBFile.writeFile(topology, state.getPositions(), f, keepIds=True)

    print(f"[run] done. trajectory: {out_dir/'traj.dcd'}  final: {out_dir/'final.pdb'}")
    return out_dir / "traj.dcd"