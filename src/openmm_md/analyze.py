"""Post-MD analysis: protein/ligand RMSD, per-residue RMSF, energy plots."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import mdtraj as md
import matplotlib

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt

# Residue names mdtraj/OpenMM give to common ions (from amber14/tip3p.xml). The
# solute-only wrapped trajectory drops these along with water.
_ION_RESNAMES = {"NA", "CL", "K", "MG", "CA", "ZN", "FE", "CU", "MN", "BR", "F", "I", "LI", "RB", "CS", "SR", "BA"}


def _solute_selection(top):
    """Protein + docked ligand (LIG) + non-standard non-water non-ion residues
    (cofactors like A3P). mdtraj misreads OpenFF "UNK" small molecules as amino
    acids, but build_system stamps the dock ligand "LIG" and cofactors with
    their PDB resname, so we key off resname rather than the `protein` flag."""
    keep = []
    for r in top.residues:
        is_water = r.name in ("HOH", "WAT", "TIP3", "SOL")
        if is_water or r.name in _ION_RESNAMES:
            continue
        # standard amino acids (protein) and any non-water/non-ion residue
        # (ligand LIG, cofactor A3P, ...) all stay.
        keep.extend(a.index for a in r.atoms)
    return keep


def _write_wrapped(traj, out_dir):
    """Reassemble molecules across PBC, keep solute only, center it."""
    wrapped = traj.image_molecules(inplace=False)
    sel = _solute_selection(wrapped.topology)
    if not sel:
        raise ValueError("no solute atoms found to wrap")
    if len(sel) < wrapped.n_atoms:
        wrapped = wrapped.atom_slice(sel)
    wrapped = wrapped.center_coordinates()
    # Solute-only xtc (compact, scalable to long runs) + a matching first-frame
    # PDB to use as the topology when loading in PyMOL/VMD:
    #   load traj_wrapped.pdb, complex ; load traj_wrapped.xtc, complex
    wrapped.save(str(out_dir / "traj_wrapped.xtc"))
    wrapped[0].save(str(out_dir / "traj_wrapped.pdb"))
    print(f"[analyze] wrote wrapped/centered solute trajectory "
          f"({wrapped.n_atoms} atoms, {wrapped.n_frames} frames) "
          f"-> {out_dir/'traj_wrapped.xtc'} (+ traj_wrapped.pdb topology)")


def analyze(traj_path: Path, topology: Path, out_dir: Path):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    traj = md.load(str(traj_path), top=str(topology))

    # The DCD is *unwrapped*: a solute that drifts across a periodic face paints
    # itself split across opposite box walls (a visualization artifact, not real
    # motion). image_molecules() reassembles whole molecules across the PBC; we
    # then keep only the solute (protein + docked ligand + cofactors, dropping
    # water/ions) and center it, and write a solute-only xtc + a matching
    # first-frame pdb. Load the .xtc in PyMOL against that traj_wrapped.pdb
    # (NOT the full complex.pdb -- a solute-only xtc needs the solute-only
    # topology or PyMOL/mdtraj reject the atom-count mismatch).
    try:
        _write_wrapped(traj, out_dir)
    except Exception as e:  # wrapping is a convenience; never block analysis on it
        print(f"[analyze] trajectory wrapping skipped ({str(e)[:80]})")

    protein_ca = traj.topology.select("protein and name CA")
    has_protein = len(protein_ca) > 0

    eng_csv = out_dir / "energy.csv"
    e = None
    if eng_csv.exists():
        # StateDataReporter header is `#"Step","Time (ps)","Potential Energy (kJ/mole)",...`
        e = pd.read_csv(eng_csv)
        e.columns = [c.lstrip("#").strip('"') for c in e.columns]

    if has_protein:
        png = _analyze_protein(traj, out_dir, protein_ca, e)
    else:
        png = _analyze_molecule(traj, out_dir, e)
    return png


def _analyze_protein(traj, out_dir, protein_ca, e):
    """Protein/ligand RMSD + per-residue CA RMSF (the original analysis)."""
    protein_heavy = traj.topology.select("protein and not element H")
    # docked ligand is stamped resname "LIG" in build_system (OpenFF otherwise
    # leaves small molecules as "UNK", which mdtraj misreads as an amino acid).
    ligand = traj.topology.select("resname LIG")

    # protein CA RMSD (optimal alignment on CA)
    rmsd_prot = (
        md.rmsd(traj, traj[0], atom_indices=protein_ca, ref_atom_indices=protein_ca)
        if len(protein_ca) else np.zeros(traj.n_frames)
    )

    # ligand RMSD in the protein-aligned frame (align on protein heavy, no
    # further re-alignment when measuring the ligand)
    if len(ligand) and len(protein_heavy):
        aligned = traj.superpose(traj[0], atom_indices=protein_heavy, ref_atom_indices=protein_heavy)
        lig0 = aligned.xyz[0, ligand]
        rmsd_lig = np.array([
            np.sqrt(((aligned.xyz[i, ligand] - lig0) ** 2).sum(axis=1).mean())
            for i in range(traj.n_frames)
        ])
    else:
        rmsd_lig = np.zeros(traj.n_frames)

    df = pd.DataFrame({
        "frame": np.arange(traj.n_frames),
        "time_ns": traj.time / 1000.0,
        "protein_ca_rmsd_nm": rmsd_prot,
        "ligand_rmsd_nm": rmsd_lig,
    })
    df.to_csv(out_dir / "rmsd.csv", index=False)

    # per-residue CA RMSF (after superposition on CA)
    if len(protein_ca):
        ca_traj = traj.superpose(traj[0], atom_indices=protein_ca, ref_atom_indices=protein_ca)
        rmsf = md.rmsf(ca_traj, ca_traj[0], atom_indices=protein_ca)
        residues = [a.residue for a in traj.topology.atoms if a.name == "CA"]
        pd.DataFrame({
            "residue": [f"{r.name}{r.resSeq}" for r in residues],
            "rmsf_nm": rmsf,
        }).to_csv(out_dir / "rmsf.csv", index=False)

    # plots
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].plot(df["time_ns"], df["protein_ca_rmsd_nm"], label="protein CA")
    ax[0].plot(df["time_ns"], df["ligand_rmsd_nm"], label="ligand")
    ax[0].set_xlabel("time (ns)")
    ax[0].set_ylabel("RMSD (nm)")
    ax[0].legend()
    ax[0].set_title("RMSD")
    _plot_energy(ax[1], e)
    fig.tight_layout()
    fig.savefig(out_dir / "analysis.png", dpi=120)
    plt.close(fig)

    print(f"[analyze] {traj.n_frames} frames. "
          f"protein CA RMSD {rmsd_prot[-1]:.3f} nm, ligand RMSD {rmsd_lig[-1]:.3f} nm "
          f"(vs frame 0). -> {out_dir/'analysis.png'}")
    return out_dir / "analysis.png"


def _analyze_molecule(traj, out_dir, e):
    """Single-molecule analysis: whole-solute RMSD, per-atom RMSF, radius of gyration.

    No protein CA -> the solute is one residue (stamped "MOL" by build_mol_system).
    RMSD/RMSF are over the solute heavy atoms (aligned on the same); Rg is over the
    whole solute (mass-weighted), a compactness/conformational proxy.
    """
    solute = np.array(_solute_selection(traj.topology))
    if not len(solute):
        raise ValueError("no solute atoms found for single-molecule analysis")
    solute_heavy = np.array(
        [i for i in solute if traj.topology.atom(i).element.symbol != "H"]
    )

    # whole-solute heavy-atom RMSD, optimally aligned on the heavy atoms
    if len(solute_heavy):
        rmsd_mol = md.rmsd(
            traj, traj[0], atom_indices=solute_heavy, ref_atom_indices=solute_heavy
        )
        # per-atom RMSF over solute heavy atoms (superposed on heavy atoms)
        mol_traj = traj.superpose(
            traj[0], atom_indices=solute_heavy, ref_atom_indices=solute_heavy
        )
        rmsf = md.rmsf(mol_traj, mol_traj[0], atom_indices=solute_heavy)
        atoms = [traj.topology.atom(i) for i in solute_heavy]
        pd.DataFrame({
            "atom": [f"{a.residue.name}{a.residue.resSeq}:{a.name}" for a in atoms],
            "element": [a.element.symbol for a in atoms],
            "rmsf_nm": rmsf,
        }).to_csv(out_dir / "rmsf.csv", index=False)
    else:
        rmsd_mol = np.zeros(traj.n_frames)

    # radius of gyration of the solute only (solvent must be excluded -- Rg over
    # the full box is meaningless). md.compute_rg has no atom_indices arg, so slice.
    rg = md.compute_rg(traj.atom_slice(solute))

    df = pd.DataFrame({
        "frame": np.arange(traj.n_frames),
        "time_ns": traj.time / 1000.0,
        "solute_rmsd_nm": rmsd_mol,
        "solute_rg_nm": rg,
    })
    df.to_csv(out_dir / "rmsd.csv", index=False)

    # plots: RMSD, Rg, Energy (3 panels for the single-molecule case)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    ax[0].plot(df["time_ns"], df["solute_rmsd_nm"])
    ax[0].set_xlabel("time (ns)"); ax[0].set_ylabel("RMSD (nm)"); ax[0].set_title("Solute RMSD")
    ax[1].plot(df["time_ns"], df["solute_rg_nm"])
    ax[1].set_xlabel("time (ns)"); ax[1].set_ylabel("Rg (nm)"); ax[1].set_title("Radius of gyration")
    _plot_energy(ax[2], e)
    fig.tight_layout()
    fig.savefig(out_dir / "analysis.png", dpi=120)
    plt.close(fig)

    print(f"[analyze] {traj.n_frames} frames (single molecule). "
          f"solute RMSD {rmsd_mol[-1]:.3f} nm, Rg {rg[-1]:.3f} nm (vs frame 0). "
          f"-> {out_dir/'analysis.png'}")
    return out_dir / "analysis.png"


def _plot_energy(ax, e):
    if e is None:
        ax.set_title("Energy (none)")
        return
    time_col = next((c for c in e.columns if c.lower().startswith("time")), e.columns[0])
    pe_col = next((c for c in e.columns if "potential" in c.lower()), None)
    if pe_col is not None:
        ax.plot(e[time_col], e[pe_col])
        ax.set_xlabel(time_col)
        ax.set_ylabel(pe_col)
        ax.set_title("Energy")