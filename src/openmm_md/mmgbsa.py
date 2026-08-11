"""MM/GBSA binding free energy via AmberTools MMPBSA.py.

Builds an Amber prmtop for the solute (receptor = protein + kept cofactors, ligand
= LIG) whose atom order *exactly* matches the production trajectory, then runs
MMPBSA.py single-trajectory GB. Reuses the production force field parameters
(ff14SB + GAFF2 via openmmforcefields, AM1-BCC) -- the same parameters the
trajectory was produced with -- so the analysis is fully consistent with the run.

Why this route (not tleap, not a straight ParmEd convert):
  * ParmEd direct-convert of the solvated system.xml fails (a known ParmEd/OpenMM
    incompatibility: with constraints=HBonds, the constrained bonds carry no
    HarmonicBondForce entry -> bond.type is None -> the prmtop writer crashes on
    `.used`). Building a fresh *solute-only* system with constraints=None avoids
    it.
  * tleap re-parameterization hits an OpenMM<->tleap impedance: tleap reorders
    protein H to its ff14SB template order and disagrees on histidine tautomers
    (the production HIE residues carry no HD1; tleap's HIE template rejects the
    extra naming), and mdtraj-suffixed non-protein atom names (C1x) don't match
    antechamber unit names. Re-protonating in tleap would desync from the
    trajectory's own H. Building the prmtop from the production solute topology
    instead keeps the production tautomers and order.
  * The ligand/cofactor residues must be re-added from OpenFF (correct bond
    orders): re-reading a PDB loses bond orders, so GAFF template matching by
    molecular formula would fail. The existing cofactor pipeline (_residue_to_rdmol
    + _explicit_h_rdmol_to_openmm, and _rdkit_to_openmm for the ligand) supplies
    them, preserving atom order.

Stage split: ``prep`` (build prmtop + convert the trajectory to netcdf + write
the MMPBSA input) is CPU-light and safe to run during a live MD job; ``run`` (the
per-frame MMPBSA.py energy evaluation) is CPU-heavy and gated behind ``--run``.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

# Ensure the env's bin/ is on PATH even when this Python was launched by absolute
# path without `conda activate`. GAFFTemplateGenerator assigns AM1-BCC charges via
# openff-toolkit's AmberToolsToolkitWrapper, which calls the `antechamber`/`sqm`
# binaries; if they aren't on PATH the wrapper reports "toolkit not available" and
# createSystem aborts. (The run step also shells out to `MMPBSA.py` from bin/.)
_env_bin = os.path.dirname(sys.executable)
if _env_bin and _env_bin not in os.environ.get("PATH", "").split(os.pathsep):
    os.environ["PATH"] = _env_bin + os.pathsep + os.environ.get("PATH", "")
# conda-forge ambertools sets AMBERHOME=$CONDA_PREFIX via an activation script; the
# `MMPBSA.py` run step needs it to find $AMBERHOME/dat/mmpbsa. Set it from the env
# prefix (parent of bin/) if activation didn't.
if not os.environ.get("AMBERHOME"):
    _env_prefix = os.path.dirname(_env_bin)
    if os.path.isdir(os.path.join(_env_prefix, "dat", "mmpbsa")):
        os.environ["AMBERHOME"] = _env_prefix

from openmm import app, unit
from openmmforcefields.generators import GAFFTemplateGenerator


def _ensure_ambertools_toolkit() -> None:
    """Make sure openff-toolkit's AmberToolsToolkitWrapper (AM1-BCC via sqm) is
    in the GLOBAL_TOOLKIT_REGISTRY.

    openff-toolkit builds that registry at first import. If the env bin/ wasn't on
    PATH at that moment (e.g. Python launched by absolute path, so antechamber/sqm
    didn't resolve), the wrapper reported "not available" and was skipped -- which
    makes GAFFTemplateGenerator's AM1-BCC charge assignment fail with "No
    registered toolkits can provide assign_partial_charges". The PATH guard above
    has now put the env bin on PATH, so registering the wrapper succeeds. No-op if
    it's already registered (normal `conda activate` + `omd` flow)."""
    from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY
    if not any(t.toolkit_name == "AmberTools"
               for t in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits):
        from openff.toolkit.utils.ambertools_wrapper import AmberToolsToolkitWrapper
        GLOBAL_TOOLKIT_REGISTRY.register_toolkit(
            AmberToolsToolkitWrapper(), exception_if_unavailable=False)


_ensure_ambertools_toolkit()

from .build_system import (
    FREE_ION_RESNAMES, WATER_RESNAMES,
    discover_hetero, parse_hetnam, _ligand_heavy_coords_ang,
    _residue_to_rdmol, _explicit_h_rdmol_to_openmm, _rdkit_to_openmm,
)
from .config import Config
from .prepare_ligand import load_ligand_rdkit


def _run(cmd: list[str], log: Path) -> None:
    """Run a command, streaming combined stdout/stderr to `log`; raise on failure."""
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "w") as f:
        f.write("$ " + " ".join(cmd) + "\n\n")
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        f.write(r.stdout.decode("utf-8", "replace"))
    if r.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed (rc={r.returncode}); see {log}")


def build_solute_prmtop(
    topology_pdb: Path, protein_pdb: Path, ligand_sdf: Path, out_dir: Path, cfg: Config,
) -> tuple[Path, Path]:
    """Build a solute-only prmtop in trajectory atom order.

    ``topology_pdb`` is the production solute topology (e.g. traj_wrapped.pdb) --
    it carries the production protein H / histidine tautomers in the trajectory's
    atom order. ``protein_pdb`` is the original (pre-solvation) protein PDB with
    CONECT records, used by the cofactor pipeline to extract each cofactor with
    correct bond orders and bound coordinates. The ligand comes from its SDF.

    LIG + cofactor residues are deleted from the solute topology and re-added from
    OpenFF (correct bond orders, same atom order as production); addHydrogens is
    NOT called (the solute PDB already has the correct protein H). constraints=None
    so every bond carries a type (ParmEd's prmtop writer needs that). NoCutoff for
    the solute-only system (no PBC; MMPBSA recomputes the MM energy anyway).
    """
    import parmed

    out_dir = Path(out_dir)
    sol = app.PDBFile(str(topology_pdb))            # production solute (right H, right order)
    top, pos = sol.topology, sol.positions

    # cofactor discovery + extraction use the ORIGINAL protein PDB (CONECT + coords)
    opdb = app.PDBFile(str(protein_pdb))
    opos = opdb.positions.value_in_unit(unit.angstrom)
    lig = load_ligand_rdkit(ligand_sdf)
    hetnam = parse_hetnam(protein_pdb)
    keep, _drop, _metals, _renames = discover_hetero(
        opdb.topology, opos, _ligand_heavy_coords_ang(lig), cfg, hetnam
    )

    lig_offmol, lig_top, lig_pos = _rdkit_to_openmm(lig, "LIG")
    offmols = [lig_offmol]
    cof_adds = []  # (topology, positions) per cofactor residue copy, in PDB order
    for r in opdb.topology.residues():
        if r.name in keep:
            off, ft, fp = _explicit_h_rdmol_to_openmm(
                _residue_to_rdmol(opdb.topology, opos, r, keep[r.name]), r.name
            )
            offmols.append(off)
            cof_adds.append((ft, fp))
    print(f"[mmgbsa] cofactor residues to re-add: {len(cof_adds)} ({list(keep)})")

    # drop water/ions/LIG/cofactors from the solute topology; re-add cofactors +
    # ligand from OpenFF. Protein residues are left untouched (their H/tautomers
    # and atom order are the production ones = the trajectory's).
    del_atoms = [
        a for r in top.residues() for a in r.atoms()
        if r.name in WATER_RESNAMES or r.name in FREE_ION_RESNAMES
        or r.name == "LIG" or r.name in keep
    ]
    mod = app.Modeller(top, pos)
    mod.delete(del_atoms)
    for ft, fp in cof_adds:
        mod.add(ft, fp)
    mod.add(lig_top, lig_pos)
    print(f"[mmgbsa] solute modeller: {mod.topology.getNumAtoms()} atoms")

    ff = app.ForceField(cfg.protein_ff, cfg.water_ff)
    gaff = GAFFTemplateGenerator(molecules=offmols, forcefield=cfg.ligand_ff)
    ff.registerTemplateGenerator(gaff.generator)
    system = ff.createSystem(
        mod.topology, nonbondedMethod=app.NoCutoff,
        constraints=None, rigidWater=False,
    )

    st = parmed.openmm.load_topology(
        mod.topology, system, xyz=mod.positions.value_in_unit(parmed.unit.angstroms)
    )
    n_none = sum(1 for b in st.bonds if b.type is None)
    if n_none:
        raise RuntimeError(f"{n_none} bonds have no type (prmtop write would crash)")
    prmtop = out_dir / "complex.prmtop"
    inpcrd = out_dir / "complex.inpcrd"
    st.save(str(prmtop))
    st.save(str(inpcrd))
    # ParmEd leaves the GB radii (RADII) at 0.0 when converting from OpenMM (which
    # uses sigma, not radii) -> MMPBSA's GB term divides by zero and returns NaN.
    # Reload as an AmberParm (load_topology gives a plain Structure with no
    # parm_data) and set mbondi2 (Bondi) radii, the set GBn2 (igb=8) expects.
    from parmed.tools import actions
    parm = parmed.load_file(str(prmtop))
    actions.changeRadii(parm, "mbondi2").execute()
    prmtop.unlink(missing_ok=True)
    parm.save(str(prmtop))
    print(f"[mmgbsa] prmtop + inpcrd -> {prmtop} ({len(parm.atoms)} atoms, "
          f"radii={parm.parm_data['RADIUS_SET']})")
    return prmtop, inpcrd


def strip_receptor_ligand_prmtops(complex_prmtop: Path, out_dir: Path) -> tuple[Path, Path]:
    """Split the complex prmtop into receptor-only and ligand-only prmtops.

    MMPBSA.py binding mode (stability=False) requires -rp/-lp prmtops whose atom
    counts sum to the complex's (complex.natom == receptor.natom + ligand.natom);
    the receptor_mask/ligand_mask in the input file then tell it how to map each
    complex trajectory frame onto them. Receptor = !:LIG (protein + cofactors),
    ligand = :LIG.
    """
    import parmed
    from parmed.tools import actions

    out_dir = Path(out_dir)
    receptor = out_dir / "receptor.prmtop"
    ligand = out_dir / "ligand.prmtop"
    for mask, out in [(":LIG", receptor), ("!:LIG", ligand)]:
        p = parmed.load_file(str(complex_prmtop))
        actions.strip(p, mask).execute()
        out.unlink(missing_ok=True)
        p.save(str(out))
    print(f"[mmgbsa] receptor + ligand prmtops -> {receptor}, {ligand}")
    return receptor, ligand


def convert_traj_nc(traj: Path, topology: Path, out_dir: Path) -> Path:
    """Re-save the wrapped trajectory as NetCDF (MMPBSA.py's preferred input).

    No reordering: the prmtop was built in the trajectory's atom order, so a plain
    format copy suffices (xtc -> nc)."""
    import mdtraj as md

    t = md.load(str(traj), top=str(topology))
    out_nc = Path(out_dir) / "traj.nc"
    t.save(str(out_nc))
    print(f"[mmgbsa] trajectory -> {out_nc} ({t.n_atoms} atoms, {t.n_frames} frames)")
    return out_nc


def write_mmpbsa_input(out_dir: Path, cfg: Config, n_frames: int, igb: int = 8) -> Path:
    """MMPBSA.py input for single-trajectory GB (GBn2, igb=8).

    Receptor = everything not LIG (protein + cofactors in the solute-only prmtop);
    ligand = LIG. PB is a later follow-up (a &pb block); not written here.
    """
    out_dir = Path(out_dir)
    inp = out_dir / "mmgbsa.in"
    inp.write_text(
        "Sample MM/GBSA input -- single trajectory, GBn2\n"
        "&general\n"
        "   verbose=2,\n"
        f"   startframe=1, endframe={n_frames}, interval=1,\n"
        "   receptor_mask=\"!:LIG\", ligand_mask=\":LIG\",\n"
        "/\n"
        "&gb\n"
        f"   igb={igb}, saltcon={cfg.ionic_strength},\n"
        "/\n"
    )
    print(f"[mmgbsa] MMPBSA.py input -> {inp}")
    return inp


def prep(
    protein_pdb: Path, ligand_sdf: Path, traj: Path, topology: Path,
    out_dir: Path, cfg: Config | None = None,
) -> dict[str, Path]:
    """Build prmtop + netcdf trajectory + MMPBSA input. Does NOT run MMPBSA.py."""
    import mdtraj as md

    cfg = cfg or Config()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    prmtop, inpcrd = build_solute_prmtop(topology, protein_pdb, ligand_sdf, out_dir, cfg)
    receptor, ligand = strip_receptor_ligand_prmtops(prmtop, out_dir)
    traj_nc = convert_traj_nc(traj, topology, out_dir)
    t = md.load(str(traj), top=str(topology))
    n_frames, _traj_n_atoms = t.n_frames, t.n_atoms
    inp = write_mmpbsa_input(out_dir, cfg, n_frames)

    # sanity: prmtop atom count must equal trajectory atom count (same source ->
    # exact order, no remap). Element sequence match was verified during dev.
    import parmed
    parm = parmed.load_file(str(prmtop))
    n = len(parm.atoms)
    if n != _traj_n_atoms:
        raise RuntimeError(
            f"prmtop has {n} atoms but trajectory has {_traj_n_atoms}; "
            "atom-order mismatch -- do NOT feed traj.nc to MMPBSA without a remap"
        )
    print(f"[mmgbsa] prmtop atoms: {n}; ligand residue present: "
          f"{any(r.name == 'LIG' for r in parm.residues)}")
    return {"prmtop": prmtop, "inpcrd": inpcrd, "receptor": receptor,
            "ligand": ligand, "traj_nc": traj_nc, "input": inp}


def run(out_dir: Path, frames: tuple[int, int] | None = None) -> Path:
    """Run MMPBSA.py GB over the prepared trajectory. CPU-heavy -- gated by --run.

    ``frames`` optionally overrides startframe/endframe (e.g. (1,1) for a 1-frame
    sanity dry-run)."""
    out_dir = Path(out_dir)
    inp = out_dir / "mmgbsa.in"
    if frames is not None:
        import re
        txt = inp.read_text()
        txt = re.sub(r"startframe=\d+", f"startframe={frames[0]}", txt)
        txt = re.sub(r"endframe=\d+", f"endframe={frames[1]}", txt)
        inp.write_text(txt)
    # Single-trajectory binding mode: -cp complex + -rp/-lp pre-stripped prmtops
    # (receptor.natom + ligand.natom == complex.natom) + the trajectory (-y). The
    # receptor_mask/ligand_mask in the input file map each complex frame onto the
    # receptor/ligand prmtops. Without -rp/-lp MMPBSA falls into stability
    # (complex-only) mode and ignores the masks.
    _run(
        ["MMPBSA.py", "-i", str(inp), "-o", str(out_dir / "FINAL_RESULTS_MMPBSA.dat"),
         "-cp", str(out_dir / "complex.prmtop"),
         "-rp", str(out_dir / "receptor.prmtop"),
         "-lp", str(out_dir / "ligand.prmtop"),
         "-y", str(out_dir / "traj.nc")],
        log=out_dir / "mmpbsa.log",
    )
    print(f"[mmgbsa] MMPBSA.py done -> {out_dir/'FINAL_RESULTS_MMPBSA.dat'}")
    return out_dir / "FINAL_RESULTS_MMPBSA.dat"