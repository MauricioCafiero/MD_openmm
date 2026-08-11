#!/usr/bin/env python
"""Assemble an OpenMM-ready 2-fragment SDF for the rot2htpuma rotaxane.

This is a *thin assembler*, not a geometry builder. The threaded, well-placed
geometry already exists on disk (built by the Rotaxanes project's pipeline:
``build_rotaxane.py`` -> ``conformer_search.py`` -> ``optimize_uma.py`` ->
``displace_wheel.py``). That pipeline outputs a plain XYZ, which has the right
3D coordinates but **no bonds** -- and OpenMM/GAFF2 needs bonds + a per-fragment
topology. So this script does only what that XYZ lacks:

  1. Build the rod and wheel as bonded RDKit mols from the two SMILES, reusing
     the Rotaxanes repo's exact embed (AddHs + ETKDGv3 seed 0xC0FFEE + MMFF) so
     the RDKit atom order **matches the XYZ atom order** (verified once: 114
     atoms, rod 58 + wheel 56, element order identical).
  2. Overlay the XYZ's coordinates onto the bonded mols (rod first, then wheel),
     so we keep the relaxed/displaced geometry *and* the bond graph.
  3. Write a 2-fragment SDF (rod as mol 0, wheel as mol 1) -- the input
     ``omd build-multimol`` consumes to give each fragment its own GAFF2
     template + residue (ROD, WHL) and solvate them as two non-bonded solutes.

The threading / PCA alignment / steric-offset minimization / GFN2 relaxation /
well-placement are NOT reimplemented here -- they live in the Rotaxanes repo
and are already baked into the XYZ we read.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# Reused verbatim from $PYLOC/Rotaxanes/code/build_rotaxane.py -- same SMILES,
# same AddHs + ETKDGv3(0xC0FFEE) + MMFF, so the atom order matches the displaced
# XYZ written by that repo's pipeline. Changing the seed or the optimizer
# would re-order atoms and break the --from-xyz overlay.
EMBED_SEED = 0xC0FFEE


def read_smiles(path: Path) -> tuple[str, str]:
    """Parse ``rod:`` / ``wheel:`` lines into (rod_smiles, wheel_smiles)."""
    rod = wheel = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("rod:"):
                rod = line[len("rod:"):].strip()
            elif line.startswith("wheel:"):
                wheel = line[len("wheel:"):].strip()
    if not rod or not wheel:
        raise ValueError(f"Could not find rod:/wheel: SMILES in {path}")
    return rod, wheel


def embed_mol(mol: Chem.Mol) -> Chem.Mol:
    """AddHs + ETKDGv3(seed) + MMFF -> a bonded mol with a 3D conformer.

    Mirrors ``mol_to_xyz_block`` in the Rotaxanes repo (AddHs first, then embed,
    then MMFF with a UFF fallback). Returns the mol (not just the element list)
    so the bond graph + conformer survive for SDF export.
    """
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = EMBED_SEED
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError("RDKit 3D embedding failed")
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        AllChem.UFFOptimizeMolecule(mol)
    return mol


def read_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    """Plain XYZ -> (element list, Nx3 coords in Angstrom)."""
    with open(path) as f:
        n = int(f.readline())
        f.readline()
        elems, coords = [], []
        for _ in range(n):
            parts = f.readline().split()
            elems.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return elems, np.array(coords, dtype=float)


def overlay_positions(mol: Chem.Mol, coords_ang: np.ndarray, label: str) -> Chem.Mol:
    """Set the mol's conformer to the given Angstrom coords (element-order
    guard is done by the caller in main(), since the rod/wheel split is known
    only there)."""
    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    if coords_ang.shape[0] != n:
        raise ValueError(f"{label}: {n} atoms in mol vs {coords_ang.shape[0]} coords from XYZ")
    for i in range(n):
        conf.SetAtomPosition(i, (float(coords_ang[i, 0]),
                                  float(coords_ang[i, 1]),
                                  float(coords_ang[i, 2])))
    return mol


def write_sdf(mols: list[Chem.Mol], path: Path, names: list[str]):
    """Write a multi-mol SDF (one MolToMolFile block per mol)."""
    with open(path, "wb") as fh:
        for mol, name in zip(mols, names):
            mol = Chem.Mol(mol)
            mol.SetProp("_Name", name)
            fh.write(Chem.MolToMolBlock(mol).encode())
            fh.write(f">  <_Name>\n{name}\n\n".encode())
            fh.write(b"$$$$\n")
    print(f"wrote {len(mols)}-fragment SDF -> {path}  "
          f"({sum(m.GetNumAtoms() for m in mols)} atoms)")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--smiles", type=Path, default=Path("rot2htpuma.txt"),
                   help="rod:/wheel: SMILES file (default: rot2htpuma.txt).")
    p.add_argument("--from-xyz", type=Path,
                   default=Path("../rot2htpuma_displaced_pinu.xyz"),
                   help="threaded/displaced XYZ whose coordinates to overlay "
                        "(bonds still come from the SMILES). "
                        "Default: ../rot2htpuma_displaced_pinu.xyz")
    p.add_argument("--out", type=Path, default=Path("rot2htpuma.sdf"),
                   help="output 2-fragment SDF (default: rot2htpuma.sdf).")
    return p.parse_args()


def main():
    args = parse_args()
    rod_smi, wheel_smi = read_smiles(args.smiles)
    print(f"rod SMILES  : {rod_smi}")
    print(f"wheel SMILES: {wheel_smi}")

    rod = embed_mol(Chem.MolFromSmiles(rod_smi))
    wheel = embed_mol(Chem.MolFromSmiles(wheel_smi))
    print(f"rod atoms  : {rod.GetNumAtoms()}   wheel atoms: {wheel.GetNumAtoms()}")

    elems, coords = read_xyz(args.from_xyz)
    print(f"xyz        : {len(elems)} atoms  ({args.from_xyz})")

    # Element-sequence guard: the XYZ must be the same build (same embed) so the
    # per-atom order matches. A mismatch means the XYZ is from a different
    # molecule / build and the overlay would put bonds on the wrong atoms.
    built = [a.GetSymbol() for a in rod.GetAtoms()] + \
            [a.GetSymbol() for a in wheel.GetAtoms()]
    if built != elems:
        # show the first divergence so the user knows the XYZ is incompatible
        for i, (a, b) in enumerate(zip(built, elems)):
            if a != b:
                raise SystemExit(
                    f"atom-order mismatch with {args.from_xyz} at index {i} "
                    f"(SMILES mol={a!r}, XYZ={b!r}); the XYZ is not from this "
                    f"SMILES/embed -- rebuild it with the Rotaxanes pipeline or "
                    f"pass a matching XYZ.")
        raise SystemExit(f"length mismatch: SMILES {len(built)} vs XYZ {len(elems)}")
    print("atom order matches XYZ -- overlaying coordinates onto bonded mols")

    n_rod = rod.GetNumAtoms()
    rod = overlay_positions(rod, coords[:n_rod], "rod")
    wheel = overlay_positions(wheel, coords[n_rod:], "wheel")

    # quick steric sanity: min rod-wheel heavy-atom distance (should be > ~2 A;
    # a dethreaded or overlapping geometry would be < 1.5 A).
    rh = np.array([list(rod.GetConformer().GetAtomPosition(i))
                   for i, a in enumerate(rod.GetAtoms()) if a.GetAtomicNum() != 1])
    wh = np.array([list(wheel.GetConformer().GetAtomPosition(i))
                   for i, a in enumerate(wheel.GetAtoms()) if a.GetAtomicNum() != 1])
    d = np.linalg.norm(rh[:, None, :] - wh[None, :, :], axis=-1)
    print(f"min rod-wheel heavy-atom distance: {d.min():.2f} A "
          f"(threads OK if > ~2 A; < 1.5 A = overlap / dethreaded)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_sdf([rod, wheel], args.out, ["rod", "wheel"])


if __name__ == "__main__":
    sys.exit(main())