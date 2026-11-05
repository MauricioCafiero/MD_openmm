"""Explicit organic-solvent boxes (acetone, DMSO, ...) via packmol.

OpenMM's ``Modeller.addSolvent`` only knows water models, so for any other solvent
the box is packed with packmol instead: the solvent molecule is built from SMILES
(RDKit ETKDGv3 + MMFF), GAFF-typed by the same GAFFTemplateGenerator as the
solutes, and N copies are packed into a cube around the fixed, centred solute.

Box edge = the same rule ``addSolvent`` uses for a cube: max(2 * solute radius +
padding, 2 * padding). N = fill * (experimental number density) * (box volume - solute
volume); fill < 1 eases packing and the NPT run restores the density. No ions are
added: the solvents here are aprotic, and a charged solute is left un-neutralized
(PME applies a uniform background charge) with a warning.

The packed box is NOT at the right density yet -- equilibrate under NPT (``omd run``
logs density in energy.csv) and check it against experiment before production.
Only cubic boxes are supported, since packmol packs a rectangular region.

To add a solvent: one entry in ``SOLVENTS`` (SMILES, residue name, experimental
density in g/cm^3 near 298 K, molar mass in g/mol).
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
from openmm import app, unit
from openff.toolkit import Molecule as OFFMol
from rdkit import Chem
from rdkit.Chem import AllChem

SOLVENTS = {
    "acetone": ("CC(=O)C", "ACN", 0.785, 58.08),
    "dmso": ("CS(=O)C", "DMS", 1.100, 78.13),
}
AVOGADRO = 6.02214076e23
HEAVY_ATOM_VOLUME_NM3 = 0.0185  # rough solute volume per heavy atom, for the count only


def solvent_offmol(name: str):
    """(OFFMol with one 3D conformer, residue name) for a solvent in SOLVENTS."""
    smi, resname, _, _ = SOLVENTS[name]
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xC0FFEE
    AllChem.EmbedMolecule(mol, params)
    AllChem.MMFFOptimizeMolecule(mol)
    return OFFMol.from_rdkit(mol, hydrogens_are_explicit=True), resname


def _write_pdb(path: Path, elements, coords_ang):
    """Minimal PDB for packmol (HETATM + element column); atom order is preserved."""
    with open(path, "w") as f:
        for i, (el, xyz) in enumerate(zip(elements, coords_ang)):
            f.write(f"HETATM{i % 100000:5d} {(el + str(i % 1000)):<4s} MOL A   1    "
                    f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {el:>2s}\n")
        f.write("END\n")


def _read_pdb_coords(path: Path) -> np.ndarray:
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in open(path) if l.startswith(("ATOM", "HETATM"))])


def pack_organic_solvent(modeller: app.Modeller, solvent: str, cfg, work_dir: Path):
    """Pack ``solvent`` around the solute already in ``modeller``.

    Returns (new Modeller with solute + solvent and the periodic box set, solvent
    OFFMol to register with the GAFF generator, number of solvent molecules).
    packmol's inputs/log are left in ``work_dir/packmol/`` for inspection.
    """
    if solvent not in SOLVENTS:
        raise ValueError(f"unknown solvent {solvent!r}; known: water, {', '.join(sorted(SOLVENTS))}")
    if cfg.box_shape != "cube":
        raise ValueError(f"--solvent {solvent} needs --box-shape cube (packmol packs a "
                         f"rectangular region); got {cfg.box_shape!r}")
    packmol = shutil.which("packmol")
    if packmol is None:
        raise RuntimeError("packmol not found on PATH (conda install -c conda-forge packmol)")

    _, _, rho, mw = SOLVENTS[solvent]
    s_off, resname = solvent_offmol(solvent)
    s_top = s_off.to_topology().to_openmm()
    for r in s_top.residues():
        r.name = resname
    s_xyz = s_off.conformers[0].m_as("angstrom")
    s_el = [a.symbol for a in s_off.atoms]

    top = modeller.topology
    xyz = np.array(modeller.positions.value_in_unit(unit.angstrom))
    el = [a.element.symbol for a in top.atoms()]
    n_heavy = sum(1 for e in el if e != "H")
    # Same box-width rule as Modeller.addSolvent (OpenMM 8.4), so --padding means the
    # same thing for every solvent: width = max(2*radius + padding, 2*padding), radius =
    # largest atom distance from the bounding-box centre.
    pad = cfg.padding * 10.0
    centre = (xyz.max(axis=0) + xyz.min(axis=0)) / 2
    radius = float(np.linalg.norm(xyz - centre, axis=1).max())
    L = max(2 * radius + pad, 2 * pad)
    xyz = xyz - centre + L / 2   # centre the solute in the box

    v_box = (L / 10) ** 3
    v_sol = HEAVY_ATOM_VOLUME_NM3 * n_heavy
    n_exp = rho / mw * AVOGADRO * (v_box - v_sol) * 1e-21
    n = int(round(cfg.solvent_fill * n_exp))
    print(f"[solvent] {solvent} ({resname}): cubic box {L:.2f} A, solute ~{v_sol:.2f} nm^3 "
          f"-> {n_exp:.0f} molecules at {rho} g/cm^3; packing {n} ({cfg.solvent_fill:.0%})")

    pk = Path(work_dir) / "packmol"
    pk.mkdir(parents=True, exist_ok=True)
    _write_pdb(pk / "solute.pdb", el, xyz)
    _write_pdb(pk / "solvent1.pdb", s_el, s_xyz)
    m = 1.0  # keep packed atoms 1 A off the faces so periodic images don't clash
    (pk / "packmol.inp").write_text(
        f"tolerance 2.0\nfiletype pdb\noutput packed.pdb\nseed 2026\n\n"
        f"structure solute.pdb\n  number 1\n  fixed 0. 0. 0. 0. 0. 0.\nend structure\n\n"
        f"structure solvent1.pdb\n  number {n}\n  inside cube {m} {m} {m} {L - 2 * m:.3f}\n"
        f"end structure\n")
    with open(pk / "packmol.inp") as fin, open(pk / "packmol.log", "w") as flog:
        subprocess.run([packmol], stdin=fin, stdout=flog, cwd=pk, check=True)
    if "Success!" not in (pk / "packmol.log").read_text():
        raise RuntimeError(f"packmol did not converge -- see {pk / 'packmol.log'}")
    packed = _read_pdb_coords(pk / "packed.pdb")
    na, ns = s_off.n_atoms, len(el)
    if len(packed) != ns + n * na:
        raise RuntimeError(f"packmol wrote {len(packed)} atoms, expected {ns + n * na}")

    out = app.Modeller(top, unit.Quantity(packed[:ns] / 10, unit.nanometer))
    for k in range(n):
        a = ns + k * na
        out.add(s_top, unit.Quantity(packed[a:a + na] / 10, unit.nanometer))
    out.topology.setPeriodicBoxVectors(np.eye(3) * L / 10 * unit.nanometer)
    return out, s_off, n
