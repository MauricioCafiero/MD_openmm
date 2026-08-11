#!/usr/bin/env python
"""Regenerate the deliverables for a finished rotaxane WT-METAD run.

Run from ``rotaxanes/metad/`` after the metad run is done (i.e. after
``metad_run/HILLS``, ``metad_run/traj.dcd`` and -- if it was resumed --
``metad_run/traj_resume.dcd`` exist). Produces, in ``metad_run/``:

  * ``Fes.dat``         -- free-energy surface F(d) from the full HILLS history
                           (shells out to ``plumed sum_hills --mintozero``).
  * ``Fes_plot.png``    -- F(d) plotted in kcal/mol, with the two terminal wells
                           (d = +/-4.6 A), the central metastable well (d = 0)
                           and the shuttle barrier marked.
  * ``pymol_rotaxane.pdb`` + ``pymol_rotaxane.dcd`` -- solute-only (114 atoms)
                           trajectory with the rod centred at the origin and the
                           wheel placed at its minimum-image position. All frames
                           are kept and there is NO PBC snapping, because the
                           wheel never leaves half a box of the rod (see README
                           "gotchas" for why the two naive wraps are wrong).

Atom layout in the built ``complex.pdb`` (0-based, matches the PLUMED 1-based
indices in ``plumed.dat``):

  * rod   = atoms 0:58   (58 atoms, residue ROD)
  * wheel = atoms 58:114 (56 atoms, residue WHL)
  * CV axis = rod amide N's at 0-based 11, 28  (PLUMED 12, 29)
  * wheel-O centroid = 0-based 58,61,64,67,70,73,76,79  (PLUMED 59,62,...,80)

Usage:
  python analyze_metad.py                 # do all three steps
  python analyze_metad.py --no-fes        # skip sum_hills + plot (PyMOL only)
  python analyze_metad.py --no-pymol      # skip the trajectory wrap (FES only)
  python analyze_metad.py --topology ../outputs/complex.pdb --out-dir metad_run

Requires the ``openmm-md`` env (mdtraj, numpy, matplotlib, and the ``plumed``
binary on PATH). Activate it first:
  source ~/miniforge3/etc/profile.d/conda.sh && conda activate openmm-md
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import numpy as np

# atom layout (0-based)
ROD = slice(0, 58)
WHL = slice(58, 114)
N1, N2 = 11, 28                 # rod amide N's defining the CV axis
WHL_O = [58, 61, 64, 67, 70, 73, 76, 79]   # the 8 wheel oxygens (every 3rd atom)
KJ_TO_KCAL = 1.0 / 4.184


def make_fes(hills: Path, out: Path,
             grid_min: float = -12.77, grid_max: float = 12.77, grid_bin: int = 360) -> Path:
    """F(d) from the full HILLS history via plumed sum_hills (--mintozero).

    The grid spans the same range as ``plumed.dat``'s METAD GRID_MIN/MAX/BIN
    (±12.77 A, 360 bins) so the FES is reproducible run-to-run rather than left to
    sum_hills' hills-derived auto-boundaries (which tighten to the sampled range
    and shift the binning). Same physical F(d) either way; this one is fixed.
    """
    cmd = ["plumed", "sum_hills", "--hills", str(hills), "--mintozero",
           "--outfile", str(out),
           "--min", str(grid_min), "--max", str(grid_max), "--bin", str(grid_bin)]
    print("[analyze] $ " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"[analyze] Fes.dat -> {out}")
    return out


def plot_fes(fes: Path, out_png: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = np.loadtxt(fes, comments="#")
    cv, F = d[:, 0], d[:, 1]
    F = F - np.nanmin(F)               # mintozero already does this, but be safe
    kcal = F * KJ_TO_KCAL
    ok = ~np.isnan(kcal)

    # shuttle barrier = max-min over the threaded region |d| <= 6.5 A
    m = ok & (np.abs(cv) <= 6.5)
    barrier = float(np.nanmax(kcal[m]) - np.nanmin(kcal[m]))

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(cv[ok], kcal[ok], lw=2, color="#1f5fa8")
    ax.set_xlabel("shuttle coordinate  d  (Å)  [wheel-O centroid, projected on N···N axis]")
    ax.set_ylabel("F(d)  (kcal/mol)")
    ax.set_title("Rotaxane free-energy surface along the rod  (WT-MetaD)")
    ax.axhline(0, color="0.7", lw=0.8)
    for dd, c in [(-4.7, "#c0392b"), (0, "#27ae60"), (4.7, "#c0392b")]:
        i = int(np.argmin(np.abs(cv - dd)))
        ax.plot(cv[i], kcal[i], "o", color=c, ms=6)
        ax.annotate(f"{kcal[i]:.1f}", (cv[i], kcal[i]),
                    textcoords="offset points", xytext=(6, 6), fontsize=9)
    ax.text(0.02, 0.95, f"barrier ≈ {barrier:.1f} kcal/mol",
            transform=ax.transAxes, va="top", fontsize=10,
            bbox=dict(fc="white", ec="0.7", boxstyle="round"))
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    print(f"[analyze] Fes_plot.png -> {out_png}  (barrier ≈ {barrier:.1f} kcal/mol)")
    return out_png


def _make_molecules_whole(P: np.ndarray, L: np.ndarray) -> np.ndarray:
    """Reference every atom of each molecule to that molecule's first atom via
    minimum-image, so each molecule is contiguous (no split across a PBC face).

    Both the rod (~13 A) and the wheel (~6 A) are smaller than L/2 (~15.8 A), so
    referencing to the first atom is exact. ``P`` is (F, N, 3) in nm; ``L`` is
    (F, 3) in nm (per-frame, NPT box may drift).
    """
    W = P.copy()
    for m0, m1 in [(ROD.start, ROD.stop), (WHL.start, WHL.stop)]:
        ref = P[:, m0:m0 + 1, :]                      # (F, 1, 3)
        delta = W[:, m0:m1, :] - ref
        delta -= L[:, None, :] * np.round(delta / L[:, None, :])
        W[:, m0:m1, :] = ref + delta
    return W


def make_pymol(trajs, topology: Path, out_dir: Path) -> tuple[Path, Path]:
    import mdtraj as md

    t = trajs[0]
    for extra in trajs[1:]:
        t = t.join(extra)
    solute = t.atom_slice(range(ROD.start, WHL.stop))
    L = np.asarray(t.unitcell_lengths)              # (F, 3) nm, cubic ~3.17

    W = _make_molecules_whole(solute.xyz.copy(), L)
    rod = W[:, ROD, :].mean(1)                       # (F, 3) whole-molecule COM
    whl = W[:, WHL, :].mean(1)
    rel = whl - rod
    rel -= L * np.round(rel / L)                     # minimum-image (PLUMED-style)

    P = W - rod[:, None, :]                         # rod -> origin
    P[:, WHL, :] += (rel - (whl - rod))[:, None, :] # wheel -> min-image rel (no snap)
    solute.xyz = P

    pdb = out_dir / "pymol_rotaxane.pdb"
    dcd = out_dir / "pymol_rotaxane.dcd"
    solute[0].save_pdb(str(pdb))
    solute.save_dcd(str(dcd))

    # report: wheel-rod distance + frame jump (sanity that there are no L-snaps)
    dvec = solute.xyz[:, WHL, :].mean(1) - solute.xyz[:, ROD, :].mean(1)
    dist = np.linalg.norm(dvec, axis=1) * 10.0
    jmp = np.linalg.norm(np.diff(dvec, axis=0), axis=1) * 10.0
    print(f"[analyze] pymol_traj -> {pdb}, {dcd}  "
          f"({solute.n_frames} frames; wheel-rod {dist.min():.2f}..{dist.max():.2f} A; "
          f"frame jump max {jmp.max():.2f} A)")
    return pdb, dcd


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--topology", type=Path, default=Path("../outputs/complex.pdb"))
    ap.add_argument("--out-dir", type=Path, default=Path("metad_run"))
    ap.add_argument("--no-fes", action="store_true", help="skip sum_hills + plot")
    ap.add_argument("--no-pymol", action="store_true", help="skip the trajectory wrap")
    args = ap.parse_args()

    out = args.out_dir
    if not args.no_fes:
        make_fes(out / "HILLS", out / "Fes.dat")
        plot_fes(out / "Fes.dat", out / "Fes_plot.png")
    if not args.no_pymol:
        import mdtraj as md
        trajs = [md.load_dcd(str(out / "traj.dcd"), top=str(args.topology))]
        resume = out / "traj_resume.dcd"
        if resume.exists():
            trajs.append(md.load_dcd(str(resume), top=str(args.topology)))
        make_pymol(trajs, args.topology, out)


if __name__ == "__main__":
    main()