"""Ligand preparation.

Two input modes:
  * SDF passthrough  -- a pre-docked pose is supplied; we only sanitize + re-save.
  * SMILES -> SDF    -- embed 5 conformers (ETKDGv3), MMFF-optimize each, keep the
                        lowest-energy conformer, write SDF.

A centroid-placement helper translates a SMILES-built ligand so its centroid sits
at a user-supplied coordinate (angstrom), used when no pre-docked pose exists.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def smiles_to_sdf(smiles: str, output_sdf: Path, num_confs: int = 5, seed: int = 42):
    """SMILES -> lowest-energy MMFF conformer -> SDF.

    Returns (output_sdf, energy_kcal_mol) of the chosen conformer.
    """
    output_sdf = Path(output_sdf)
    output_sdf.parent.mkdir(parents=True, exist_ok=True)

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"could not parse SMILES: {smiles!r}")
    mol = Chem.AddHs(mol)  # explicit H needed for embedding + MMFF + GAFF

    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
    if not cids:
        raise RuntimeError("RDKit 3D embedding failed")

    # MMFF optimize every conformer. Returns (result, energy) per conformer,
    # where result == 0 means the optimization converged.
    results = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, maxIters=500)
    energies = [e if (r == 0) else float("inf") for (r, e) in results]
    best_pos = int(np.argmin(energies))
    best = cids[best_pos]            # conformer id (EmbedMultipleConfs ids need not be 0)
    best_energy = energies[best_pos]

    Chem.SanitizeMol(mol)
    # write the chosen conformer by its actual id -- Chem.Mol(mol, confId=best)
    # keeps the conformer at id `best` (it is NOT renumbered to 0), so writing
    # confId=0 would raise "Bad Conformer Id" whenever best != 0.
    Chem.MolToMolFile(mol, str(output_sdf), confId=best, kekulize=True)
    print(f"[prep-ligand] SMILES -> {output_sdf} (conf {best}, "
          f"MMFF energy = {best_energy:.3f} kcal/mol)")
    return output_sdf, best_energy


def load_ligand_rdkit(ligand_sdf: Path) -> Chem.Mol:
    """Read an SDF/MOL into an RDKit Mol, normalizing bad valence fields.

    Some exporters (notably OpenBabel) write MDL valence / hydrogen-count fields
    that set ``NoImplicit`` on atoms *and* omit the corresponding explicit H
    atoms. RDKit then perceives radical electrons instead of implicit hydrogens,
    which OpenFF rejects. We read without sanitizing, clear ``NoImplicit`` so
    implicit H is recomputed from bond orders, then sanitize. The 3D coordinates
    of the supplied atoms (i.e. the docked pose) are preserved.
    """
    mol = Chem.MolFromMolFile(str(ligand_sdf), removeHs=False, sanitize=False)
    if mol is None:
        # not an MDL-valence issue; try a plain sanitized read for a clearer error
        mol = Chem.MolFromMolFile(str(ligand_sdf), removeHs=False, sanitize=True)
        if mol is None:
            raise ValueError(f"could not read SDF: {ligand_sdf}")
        return mol
    for a in mol.GetAtoms():
        a.SetNoImplicit(False)
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol)
    return mol


def centroid_angstrom(mol: Chem.Mol, conf_id: int = 0) -> np.ndarray:
    conf = mol.GetConformer(conf_id)
    pos = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
        dtype=float,
    )
    return pos.mean(axis=0)


def translate_to(mol: Chem.Mol, target_xyz_ang, conf_id: int = 0) -> Chem.Mol:
    """Translate the ligand so its centroid sits at ``target_xyz_ang`` (angstrom)."""
    target = np.asarray(target_xyz_ang, dtype=float)
    delta = target - centroid_angstrom(mol, conf_id)
    conf = mol.GetConformer(conf_id)
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (p.x + delta[0], p.y + delta[1], p.z + delta[2]))
    return mol


def write_sdf(mol: Chem.Mol, output_sdf: Path, conf_id: int = 0) -> Path:
    output_sdf = Path(output_sdf)
    output_sdf.parent.mkdir(parents=True, exist_ok=True)
    Chem.MolToMolFile(mol, str(output_sdf), confId=conf_id, kekulize=True)
    return output_sdf