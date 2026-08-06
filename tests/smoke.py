"""Smoke tests for the openmm-md pipeline.

Run directly:  ``python tests/smoke.py``
or with pytest if installed.

Always runs:
  * ligand generation from SMILES (ethanol) -> SDF
  * platform probing (auto)

Runs the full prep-protein -> build -> run -> analyze pipeline only when the
user provides real inputs via env vars:
  OMD_PROTEIN_PDB  path to a protein PDB
  OMD_LIGAND_SDF   path to a ligand SDF (pre-docked pose)  OR
  OMD_LIGAND_SMILES  SMILES to build the ligand from
  OMD_SITE  optional "x y z" angstrom centroid for a SMILES ligand
  OMD_COFACTORS  optional comma-separated cofactor resnames to keep + GAFF-type
                 (resolved against build_system.KNOWN_COFACTORS, e.g. A3P)
"""
from __future__ import annotations

import os
import tempfile
from pathlib import Path


def test_ligand_from_smiles():
    from openmm_md.prepare_ligand import smiles_to_sdf, load_ligand_rdkit
    from rdkit import Chem

    with tempfile.TemporaryDirectory() as d:
        out, energy = smiles_to_sdf("CCO", Path(d) / "etoh.sdf", num_confs=5)
        assert out.exists()
        mol = load_ligand_rdkit(out)
        assert mol.GetNumAtoms() > 3, "ethanol should have >3 atoms incl. H"
        assert mol.GetNumConformers() == 1
        assert Chem.MolToSmiles(Chem.RemoveHs(mol)) == "CCO"
    print("  [ok] ligand_from_smiles")


def test_platform_probe():
    from openmm_md.dynamics import get_platform
    p = get_platform("auto")
    assert p is not None
    print(f"  [ok] platform_probe -> {p.getName()}")


def test_full_pipeline():
    pdb = os.environ.get("OMD_PROTEIN_PDB")
    sdf = os.environ.get("OMD_LIGAND_SDF")
    smi = os.environ.get("OMD_LIGAND_SMILES")
    if not pdb or not (sdf or smi):
        print("  [skip] full_pipeline (set OMD_PROTEIN_PDB + OMD_LIGAND_SDF/OMD_LIGAND_SMILES to run)")
        return
    site = None
    if os.environ.get("OMD_SITE"):
        site = [float(v) for v in os.environ["OMD_SITE"].replace(",", " ").split()]

    from openmm_md.config import Config
    from openmm_md.build_system import build_system
    from openmm_md.dynamics import run as run_dynamics
    from openmm_md.analyze import analyze

    cfg = Config()
    cfg.equilibrate_steps = 200
    cfg.production_steps = 500          # short smoke run
    cfg.report_interval = 50
    cfg.traj_interval = 50
    cfg.checkpoint_interval = 250
    cfg.padding = 1.0

    cof = os.environ.get("OMD_COFACTORS", "")
    if cof:
        from openmm_md.build_system import KNOWN_COFACTORS
        cfg.cofactors = {}
        for name in cof.split(","):
            name = name.strip()
            if name not in KNOWN_COFACTORS:
                raise KeyError(f"OMD_COFACTORS: {name!r} not in KNOWN_COFACTORS")
            cfg.cofactors[name] = KNOWN_COFACTORS[name]
        cfg.auto_cofactors = False  # explicit-only: keep the smoke run offline/deterministic

    # Persist artifacts to outputs/smoke/ (override with OMD_OUT_DIR) so they can
    # be inspected after the run, instead of a throwaway temp dir.
    d = Path(os.environ.get("OMD_OUT_DIR", "outputs/smoke"))
    if d.exists():
        import shutil
        shutil.rmtree(d)
    d.mkdir(parents=True, exist_ok=True)
    lig = Path(sdf) if sdf else (d / "lig.sdf")
    if smi:
        from openmm_md.prepare_ligand import smiles_to_sdf, translate_to, write_sdf
        _, _ = smiles_to_sdf(smi, lig, num_confs=5)
        if site is not None:
            from openmm_md.prepare_ligand import load_ligand_rdkit
            m = load_ligand_rdkit(lig); translate_to(m, site); write_sdf(m, lig)
    sys_xml, top = build_system(pdb, lig, d, cfg=cfg, site_xyz_ang=site)
    assert sys_xml.exists() and top.exists()
    traj = run_dynamics(sys_xml, top, d, cfg=cfg, steps=cfg.production_steps)
    assert traj.exists()
    png = analyze(traj, top, d)
    assert png.exists()
    print(f"  [ok] full_pipeline (artifacts in {d})")


def test_single_molecule():
    """Solvate a single small molecule (no protein) and run a short MD + analysis.

    Always runs (no inputs needed): builds paracetamol from SMILES, solvates it
    alone in TIP3P + ions, runs a 200-step production, and analyzes. Exercises the
    build_mol_system -> dynamics (restraints auto-disabled) -> analyze (single-
    molecule branch) path end to end.
    """
    from openmm_md.config import Config
    from openmm_md.build_system import build_mol_system
    from openmm_md.dynamics import run as run_dynamics
    from openmm_md.analyze import analyze

    cfg = Config()
    cfg.equilibrate_steps = 100
    cfg.production_steps = 200
    cfg.report_interval = 20
    cfg.traj_interval = 20
    cfg.checkpoint_interval = 100
    cfg.padding = 1.0
    cfg.restrain_protein = True  # must be auto-disabled (no protein) by dynamics.run
    cfg.pressure = 0.0  # NVT: a ~230-water box is small enough that NPT pressure
    # noise makes the barostat shrink the box below 2*nonbonded_cutoff. NVT
    # validates the build/run/analyze path; real single-molecule runs should use
    # a larger box (padding >= 1.5 nm) and/or staged NVT->NPT equilibration.
    cfg.platform = "CPU"  # tiny system; also avoids OpenCL contention with any
    # concurrent protein/ligand production run sharing the GPU

    d = Path(os.environ.get("OMD_MOL_OUT_DIR", "outputs/smoke_mol"))
    if d.exists():
        import shutil
        shutil.rmtree(d)
    d.mkdir(parents=True, exist_ok=True)

    lig = d / "lig.sdf"
    from openmm_md.prepare_ligand import smiles_to_sdf
    smiles_to_sdf("CC(=O)Nc1ccc(O)cc1", lig, num_confs=3)  # paracetamol

    sys_xml, top = build_mol_system(lig, d, cfg=cfg)
    assert sys_xml.exists() and top.exists()
    traj = run_dynamics(sys_xml, top, d, cfg=cfg, steps=cfg.production_steps)
    assert traj.exists()
    png = analyze(traj, top, d)
    assert png.exists()
    # single-molecule analysis writes solute RMSD + Rg, not protein CA RMSD
    import csv
    with open(d / "rmsd.csv") as f:
        cols = next(csv.reader(f))
    assert "solute_rmsd_nm" in cols and "solute_rg_nm" in cols, cols
    print(f"  [ok] single_molecule (artifacts in {d})")


def test_auto_cofactors():
    """Auto cofactor discovery: parse HETNAM, classify hetero by overlap with the
    docked ligand. On SULT1A3 (dimer + PAP + L-dopamine), PAP (A3P) is kept and
    L-dopamine (LDP) is dropped from BOTH monomers (per-resname rule: any residue
    of a resname overlapping the docked SDF drops the whole resname). Offline --
    A3P is a built-in cofactor and LDP is dropped, so no PubChem lookup is needed.
    """
    from openmm import app, unit
    from openmm_md.build_system import (
        parse_hetnam, discover_hetero, _ligand_heavy_coords_ang,
    )
    from openmm_md.config import Config
    from openmm_md.prepare_ligand import load_ligand_rdkit

    root = Path(__file__).resolve().parent.parent
    pdb_path = root / "data" / "protein" / "sult1a3_2A3R.pdb"
    sdf_path = root / "data" / "ligands" / "sult1a3_2A3R_c0.sdf"
    if not pdb_path.exists() or not sdf_path.exists():
        print("  [skip] auto_cofactors (SULT1A3 data files not found)")
        return

    het = parse_hetnam(pdb_path)
    assert "A3P" in het and "LDP" in het, het
    assert "DIPHOSPHATE" in het["A3P"].upper()
    assert "DOPAMINE" in het["LDP"].upper()

    pdb = app.PDBFile(str(pdb_path))
    pos_ang = pdb.positions.value_in_unit(unit.angstrom)
    lig = _ligand_heavy_coords_ang(load_ligand_rdkit(sdf_path))
    cfg = Config()  # auto_cofactors=True, threshold 5.0 A by default
    keep, drop, metals, renames = discover_hetero(pdb.topology, pos_ang, lig, cfg, het)
    assert "A3P" in keep and "A3P" not in drop, (keep, drop)
    assert "LDP" in drop and "LDP" not in keep, (keep, drop)
    assert metals == [] and renames == {}, (metals, renames)  # SULT1A3 has no metal
    print(f"  [ok] auto_cofactors (keep={list(keep)}, drop={list(drop)})")


def test_structural_metal():
    """A bound structural Zn2+ is kept (not dropped) and locked into its crystal
    Cys2His2 coordination site: typed +2 by the water FF, with 4 metal-ligand bonds
    + 6 ligand-metal-ligand angles + nonbonded exclusions, and the coordinating
    Cys -> CYM (thiolate) / His -> HID (coordinating NE2 free) re-protonated.
    Build-only (no MD) on the 1ZNF zinc-finger fixture. CPU to avoid OpenCL contention.
    """
    from openmm import app, unit, XmlSerializer
    from openmm_md.config import Config
    from openmm_md.build_system import build_system
    from openmm_md.prepare_ligand import smiles_to_sdf
    from collections import Counter

    root = Path(__file__).resolve().parent.parent
    pdb_path = root / "data" / "protein" / "1znf_zinc_finger.pdb"
    if not pdb_path.exists():
        print("  [skip] structural_metal (1ZNF fixture not found)")
        return

    d = Path(os.environ.get("OMD_METAL_OUT_DIR", "outputs/smoke_metal"))
    if d.exists():
        import shutil
        shutil.rmtree(d)
    d.mkdir(parents=True, exist_ok=True)

    # dummy ligand (ethanol) placed at the protein centroid -- far from the Zn site;
    # the ligand only matters for the build signature, not the metal handling.
    pdb = app.PDBFile(str(pdb_path))
    pos = pdb.positions.value_in_unit(unit.angstrom)
    import numpy as np
    heavy = [pos[a.index] for a in pdb.topology.atoms() if a.element and a.element.symbol != "H"]
    centroid = np.mean([[p[0], p[1], p[2]] for p in heavy], axis=0)
    lig = d / "lig.sdf"
    smiles_to_sdf("CCO", lig, num_confs=3)

    cfg = Config()
    cfg.padding = 1.0
    cfg.platform = "CPU"
    sys_xml, top = build_system(pdb_path, lig, d, cfg=cfg, site_xyz_ang=centroid.tolist())
    assert sys_xml.exists() and top.exists()

    # coordinating residues re-protonated in the built complex: 2 CYS -> CYM (thiolate,
    # no HG), and the 2 NE2-coordinating HIS forced to the HID tautomer (H on ND1,
    # NE2 free). addHydrogens keeps the HIS *name* as HIS (it protonates by composition,
    # matching the existing pipeline), so check protonation by H placement not name.
    built = app.PDBFile(str(top))
    rn = Counter(r.name for r in built.topology.residues())
    assert rn.get("CYM", 0) == 2, f"expected 2 CYM (thiolate), got {rn.get('CYM')}"
    assert rn.get("ZN", 0) == 1, f"expected 1 ZN kept, got {rn.get('ZN')}"
    for r in built.topology.residues():
        if r.name == "HIS" and r.id in ("19", "23"):
            names = {a.name for a in r.atoms()}
            assert "HD1" in names and "HE2" not in names, f"HIS{r.id} not HID-protonated: {names}"

    # Zn charge +2, coordination bonds/angles added
    system = XmlSerializer.deserialize(sys_xml.read_text())
    nb = next(f for f in system.getForces() if f.__class__.__name__ == "NonbondedForce")
    zn_idx = next(a.index for r in built.topology.residues() if r.name == "ZN" for a in r.atoms())
    q, _, _ = nb.getParticleParameters(zn_idx)
    assert abs(q.value_in_unit(unit.elementary_charge) - 2.0) < 1e-6, "Zn should be +2"
    n_bonds = n_angles = 0
    for f in system.getForces():
        if f.__class__.__name__ == "HarmonicBondForce":
            for i in range(f.getNumBonds()):
                p1, p2, _, _ = f.getBondParameters(i)
                if zn_idx in (p1, p2):
                    n_bonds += 1
        elif f.__class__.__name__ == "HarmonicAngleForce":
            for i in range(f.getNumAngles()):
                p1, p2, p3, _, _ = f.getAngleParameters(i)
                if zn_idx == p2:  # apex is the metal
                    n_angles += 1
    assert n_bonds == 4, f"expected 4 Zn-ligand bonds, got {n_bonds}"
    assert n_angles == 6, f"expected 6 ligand-Zn-ligand angles, got {n_angles}"
    print(f"  [ok] structural_metal (Zn+2 kept, {n_bonds} bonds + {n_angles} angles, "
          f"2 CYM + 2 HID; artifacts in {d})")


def test_zinc_finger_md():
    """Full build -> run -> analyze on the 1ZNF zinc finger + a placed ligand, with a
    bound Zn2+ locked in its Cys2His2 site by the bonded model. The key robustness
    check: the Zn-ligand distances stay ~2.3 A through dynamics (the metal does not
    drift out of its crystal coordination site). CPU, short run, no OpenCL contention.
    """
    from openmm import app, unit
    import numpy as np
    from openmm_md.config import Config
    from openmm_md.build_system import build_system
    from openmm_md.dynamics import run as run_dynamics
    from openmm_md.analyze import analyze
    from openmm_md.prepare_ligand import smiles_to_sdf

    root = Path(__file__).resolve().parent.parent
    pdb_path = root / "data" / "protein" / "1znf_zinc_finger.pdb"
    if not pdb_path.exists():
        print("  [skip] zinc_finger_md (1ZNF fixture not found)")
        return

    d = Path(os.environ.get("OMD_ZNMD_OUT_DIR", "outputs/smoke_znmd"))
    if d.exists():
        import shutil
        shutil.rmtree(d)
    d.mkdir(parents=True, exist_ok=True)

    # place a proxy ligand just outside the protein (no steric clash with the finger)
    pdb = app.PDBFile(str(pdb_path))
    pos = pdb.positions.value_in_unit(unit.angstrom)
    xs = [pos[a.index][0] for a in pdb.topology.atoms() if a.element and a.element.symbol != "H"]
    ys = [pos[a.index][1] for a in pdb.topology.atoms() if a.element and a.element.symbol != "H"]
    site = [max(xs) + 5.0, float(np.mean(ys)), float(np.mean(ys))]
    lig = d / "lig.sdf"
    smiles_to_sdf("c1ccccc1O", lig, num_confs=3)  # phenol as a proxy ligand

    cfg = Config()
    cfg.equilibrate_steps = 100
    cfg.production_steps = 300
    cfg.report_interval = 30
    cfg.traj_interval = 30
    cfg.checkpoint_interval = 150
    cfg.padding = 1.0
    cfg.platform = "CPU"

    sys_xml, top = build_system(pdb_path, lig, d, cfg=cfg, site_xyz_ang=site)
    traj = run_dynamics(sys_xml, top, d, cfg=cfg, steps=cfg.production_steps)
    assert traj.exists()
    png = analyze(traj, top, d)
    assert png.exists()

    # robustness check: Zn still coordinated (~2.3 A) in the last frame. mdtraj's
    # Topology.residues is a one-shot generator and Residue.atoms is a property, so
    # materialize the residue list and iterate .atoms directly (no call).
    import mdtraj as md
    t = md.load(str(traj), top=str(top))
    residues = list(t.topology.residues)
    zn = next(a for r in residues if r.name == "ZN" for a in r.atoms)
    # the coordinating Cys were renamed CYM (thiolate) in the built topology; HIS keeps
    # its name (protonation by composition -> HID tautomer, NE2 free).
    ligand_atoms = {("CYM", 3, "SG"), ("CYM", 6, "SG"), ("HIS", 19, "NE2"), ("HIS", 23, "NE2")}
    ligs = [a for r in residues for a in r.atoms
            if (r.name, int(r.resSeq), a.name) in ligand_atoms]
    assert len(ligs) == 4, f"expected 4 Zn ligands, found {len(ligs)}"
    last = t.xyz[-1]  # nm
    dists = [float(np.linalg.norm(last[zn.index] - last[a.index]) * 10.0) for a in ligs]  # A
    assert all(1.8 < d_ < 3.0 for d_ in dists), f"Zn drifted: {dists}"
    print(f"  [ok] zinc_finger_md (Zn-ligand {[round(d_,2) for d_ in dists]} A in last frame; "
          f"metal held; artifacts in {d})")


def main():
    for fn in (test_ligand_from_smiles, test_platform_probe, test_single_molecule,
               test_auto_cofactors, test_structural_metal, test_zinc_finger_md,
               test_full_pipeline):
        print(f"running {fn.__name__} ...")
        fn()
    print("ALL SMOKE TESTS PASSED")


if __name__ == "__main__":
    main()