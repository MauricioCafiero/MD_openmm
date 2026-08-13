#!/usr/bin/env python
"""Extract representative MD snapshots at target shuttle-coordinate (d) values
from the finished rot1 metad run, for QM single-point rescoring.

Recomputes d directly from the solute-only pymol_rotaxane.{pdb,dcd} trajectory
(rod N1/N2 midpoint -> wheel-O centroid, projected on the N1->N2 axis -- same
definition PLUMED used in ../metad/make_plumed.py) rather than back-mapping to
COLVAR rows, since COLVAR's time axis resets across the resume segment (see
../metad/README.md's "time is in fs" gotcha) while the joined PyMOL trajectory
(pymol_rotaxane.dcd, produced by ../metad/analyze_metad.py) already has frames
in correct end-to-end order across both run legs.

For each target d, writes the single nearest frame as a plain XYZ
(frames/d{d:+.2f}.xyz) for scoring with score_tblite.py.

``--auto-grid`` builds the target list from the MetaD's own Fes.dat instead of
a manual list: dense spacing across the physically interesting span (both
terminal wells + the bumpy central plateau between them -- rot1's landscape
has ~7 minima and ~5 maxima in that span, not a simple twin-well, so "dense at
the wells" in practice means dense across the whole span they occupy), sparse
spacing in the tails beyond the terminal wells (nothing physically interesting
out there, and it's the least-converged part of the sampled range -- see
../metad/README.md's "barrier window" gotcha).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import mdtraj as md
import numpy as np
from scipy.signal import argrelextrema

# rot1 atom layout (0-based), from ../metad/analyze_metad.py's atom_layout()
# on ../../outputs/complex.pdb: ROD=0:88 WHL=88:144 N1,N2=27,46
# WHL_O=[88,91,94,97,100,103,106,109]. Duplicated here (not imported) to keep
# this folder self-contained; re-derive via atom_layout() if the build changes.
N1, N2 = 27, 46
WHL_O = [88, 91, 94, 97, 100, 103, 106, 109]

# the four wells + an approximate saddle from the final Fes_plot.png
# (-10.85/+11.20 terminal wells, -4.90/+4.20 intermediate features, ~+7 near
# the true barrier peak identified between the station wells). Used when
# --auto-grid is not passed.
DEFAULT_TARGETS = [-10.85, -4.90, 0.0, 4.20, 7.00, 11.20]


def build_auto_grid(fes_path: Path, colvar_path: Path,
                    dense_min: float, dense_max: float,
                    dense_spacing: float, sparse_spacing: float) -> list[float]:
    """Dense grid across [dense_min, dense_max] (default: spans both terminal
    wells + the bumpy central plateau), sparse grid in the tails out to the
    observed COLVAR range. Values from Fes.dat aren't used to place individual
    points (rot1's landscape is bumpy enough -- 7 minima + 5 maxima in the
    reliable range -- that per-feature clustering would just re-derive
    "dense everywhere in the middle" anyway); Fes.dat is only used to report
    what's in that span, for the printed summary."""
    fes = np.loadtxt(fes_path, comments="#")
    colvar = np.loadtxt(colvar_path, comments="#")
    obs_min, obs_max = float(colvar[:, 1].min()), float(colvar[:, 1].max())
    cv, F = fes[:, 0], fes[:, 1]
    m = (cv >= obs_min) & (cv <= obs_max) & ~np.isnan(F)
    cv_m, kcal = cv[m], (F[m] - F[m].min()) / 4.184
    k = 5
    smooth = np.convolve(kcal, np.ones(k) / k, mode="same")
    order = max(3, len(smooth) // 40)
    mins = sorted(round(float(cv_m[i]), 2) for i in argrelextrema(smooth, np.less_equal, order=order)[0])
    maxs = sorted(round(float(cv_m[i]), 2) for i in argrelextrema(smooth, np.greater_equal, order=order)[0])
    print(f"[auto-grid] observed range {obs_min:.2f}..{obs_max:.2f} A")
    print(f"[auto-grid] minima in dense span: {[x for x in mins if dense_min <= x <= dense_max]}")
    print(f"[auto-grid] maxima in dense span: {[x for x in maxs if dense_min <= x <= dense_max]}")

    dense = list(np.arange(dense_min, dense_max + dense_spacing / 2, dense_spacing))
    sparse_lo = list(np.arange(obs_min, dense_min - sparse_spacing / 2, sparse_spacing))
    sparse_hi = list(np.arange(dense_max + sparse_spacing, obs_max + sparse_spacing / 2, sparse_spacing))
    grid = sorted(round(float(x), 2) for x in sparse_lo + dense + sparse_hi)
    print(f"[auto-grid] {len(sparse_lo)} sparse (lo) + {len(dense)} dense + "
          f"{len(sparse_hi)} sparse (hi) = {len(grid)} target points")
    return grid


def compute_d(traj: md.Trajectory) -> np.ndarray:
    """d = (wheel-O centroid - midpoint(N1,N2)) . (N2-N1)/|N2-N1|, in Angstrom.
    Same definition as ../metad/make_plumed.py's PLUMED CV, computed directly
    from Cartesian coordinates -- no PBC handling needed here since
    pymol_rotaxane.dcd is already unwrapped/whole and solute-only."""
    xyz = traj.xyz * 10.0  # nm -> A
    n1, n2 = xyz[:, N1, :], xyz[:, N2, :]
    mid = (n1 + n2) / 2.0
    whl = xyz[:, WHL_O, :].mean(axis=1)
    axis = n2 - n1
    axis_hat = axis / np.linalg.norm(axis, axis=1, keepdims=True)
    return np.einsum("ij,ij->i", whl - mid, axis_hat)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--topology", type=Path,
                    default=Path("../metad/metad_run/pymol_rotaxane.pdb"))
    ap.add_argument("--traj", type=Path,
                    default=Path("../metad/metad_run/pymol_rotaxane.dcd"))
    ap.add_argument("--targets", type=float, nargs="+", default=DEFAULT_TARGETS,
                    help="target d values in A (default: the identified wells "
                         "+ approximate saddle from the final Fes_plot.png); "
                         "ignored if --auto-grid is passed")
    ap.add_argument("--auto-grid", action="store_true",
                    help="build the target list from Fes.dat instead of --targets: "
                         "dense across [--dense-min, --dense-max], sparse elsewhere")
    ap.add_argument("--fes", type=Path, default=Path("../metad/metad_run/Fes.dat"))
    ap.add_argument("--colvar", type=Path, default=Path("../metad/metad_run/COLVAR"))
    ap.add_argument("--dense-min", type=float, default=-12.0)
    ap.add_argument("--dense-max", type=float, default=12.0)
    ap.add_argument("--dense-spacing", type=float, default=0.5)
    ap.add_argument("--sparse-spacing", type=float, default=1.5)
    ap.add_argument("--out-dir", type=Path, default=Path("frames"))
    ap.add_argument("--n-frames", type=int, default=4,
                    help="MD frames to extract per target d, for later averaging "
                         "(default 4; 1 = old single-frame behavior)")
    ap.add_argument("--window", type=float, default=0.15,
                    help="A -- candidate frames must be within this of the target")
    ap.add_argument("--min-gap", type=int, default=500,
                    help="minimum frame-index separation between picked replicates, "
                         "so they sample different excursions rather than adjacent, "
                         "highly-correlated MD steps (500 frames = 1 ns apart here)")
    args = ap.parse_args()

    if args.auto_grid:
        targets = build_auto_grid(args.fes, args.colvar, args.dense_min, args.dense_max,
                                  args.dense_spacing, args.sparse_spacing)
    else:
        targets = args.targets

    t = md.load_dcd(str(args.traj), top=str(args.topology))
    d = compute_d(t)
    print(f"loaded {t.n_frames} frames; d range {d.min():+.2f} .. {d.max():+.2f} A")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    elements = [a.element.symbol for a in t.topology.atoms]
    for target in targets:
        cand = np.where(np.abs(d - target) < args.window)[0]
        cand = cand[np.argsort(np.abs(d[cand] - target))]  # closest-to-target first
        picked: list[int] = []
        for i in cand:
            if all(abs(int(i) - j) >= args.min_gap for j in picked):
                picked.append(int(i))
            if len(picked) == args.n_frames:
                break
        if not picked:  # window too tight for a sparse region -- fall back to nearest
            picked = [int(np.argmin(np.abs(d - target)))]

        for rep, i in enumerate(picked):
            xyz = t.xyz[i] * 10.0  # A
            suffix = f"_r{rep}" if args.n_frames > 1 else ""
            out = args.out_dir / f"d{target:+.2f}{suffix}.xyz"
            with open(out, "w") as f:
                f.write(f"{len(elements)}\n")
                f.write(f"rot1 snapshot, target d={target:+.2f} A, actual d={d[i]:+.2f} A, frame {i}\n")
                for el, (x, y, z) in zip(elements, xyz):
                    f.write(f"{el:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")
        print(f"  target {target:+6.2f} A -> {len(picked)} frame(s) {picked} "
              f"(d={[round(float(d[i]),2) for i in picked]})")


if __name__ == "__main__":
    main()
