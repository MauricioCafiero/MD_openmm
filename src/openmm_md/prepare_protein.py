"""Protein preparation with PDBFixer."""
from __future__ import annotations

from pathlib import Path

from openmm import app
from pdbfixer import PDBFixer


def prepare_protein(
    input_pdb: Path,
    output_pdb: Path,
    ph: float = 7.0,
    add_missing_residues: bool = False,
) -> Path:
    """Repair a protein PDB and write a clean, hydrogenated structure.

    By default only missing *atoms* are patched (not whole missing loops), which
    is what you want for a structured binding-site protein. Set
    ``add_missing_residues=True`` to let PDBFixer model missing loops too.
    """
    input_pdb = Path(input_pdb)
    output_pdb = Path(output_pdb)
    output_pdb.parent.mkdir(parents=True, exist_ok=True)

    fixer = PDBFixer(filename=str(input_pdb))
    fixer.findMissingResidues()
    if not add_missing_residues:
        # patch atoms only; don't invent terminal loops
        fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=ph)

    with open(output_pdb, "w") as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

    n_res = sum(1 for r in fixer.topology.residues())
    n_atoms = fixer.topology.getNumAtoms()
    print(f"[prep-protein] {input_pdb.name}: {n_res} residues, {n_atoms} atoms -> {output_pdb}")
    return output_pdb