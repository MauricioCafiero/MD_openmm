#!/usr/bin/env python
"""Snapshot the shuttle FES from an in-progress WT-METAD HILLS file, for a
periodic check-in while the run is still going.

Reports three things, each from ``plumed sum_hills``:
  * ALL HILLS    -- F(d) from every hill deposited so far.
  * FIRST HALF   -- F(d) from only the first half of hills (by deposition order).
  * SECOND HALF  -- F(d) from only the second half.

FIRST HALF vs SECOND HALF is a rough convergence diagnostic: if their well
positions/barrier already agree, the bias has likely flattened the surface and
more sampling is mostly just refining noise; if they disagree a lot, the run
needs more time. This mirrors the standard WT-MetaD block-average convergence
check, done on two blocks instead of many (cheap, one-shot per check-in).

Wells are found as local minima of the (lightly smoothed) F(d) restricted to
the actually-sampled cv range (from COLVAR min/max) -- unsampled tail bins sit
at the flat unbiased grid edge and would otherwise look like spurious wells.

Run from rotaxanes/metad/ while run_metad.py is still writing HILLS/COLVAR:
  python checkin_fes.py
"""
from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy.signal import argrelextrema

KJ_TO_KCAL = 1.0 / 4.184


def split_hills(hills_path: Path) -> tuple[list[str], list[str], list[str]]:
    lines = hills_path.read_text().splitlines(keepends=True)
    header = [l for l in lines if l.startswith("#")]
    data = [l for l in lines if not l.startswith("#") and l.strip()]
    half = len(data) // 2
    return header, data[:half], data[half:]


def sum_hills(text: str, grid_min: float, grid_max: float, grid_bin: int,
              tmpdir: Path, tag: str) -> Path | None:
    if not text.strip():
        return None
    hills_path = tmpdir / f"HILLS_{tag}"
    hills_path.write_text(text)
    out_path = tmpdir / f"Fes_{tag}.dat"
    cmd = ["plumed", "sum_hills", "--hills", str(hills_path), "--mintozero",
           "--outfile", str(out_path), "--min", str(grid_min), "--max", str(grid_max),
           "--bin", str(grid_bin)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not out_path.exists():
        return None
    return out_path


def analyze(fes_path: Path, obs_min: float, obs_max: float):
    """Wells (local minima) over the whole sampled range, but the barrier is
    restricted to the window BETWEEN the two station wells -- not the whole
    sampled range. The tails beyond each station (out toward obs_min/obs_max)
    are the least-visited, least-flattened parts of the surface and can carry
    a spuriously steep, unconverged rise that has nothing to do with the real
    shuttling transition state (caught in practice: a reported "barrier" that
    was actually the sampled-range edge, not the saddle between the wells)."""
    d = np.loadtxt(fes_path, comments="#")
    cv, F = d[:, 0], d[:, 1]
    m = (cv >= obs_min) & (cv <= obs_max) & ~np.isnan(F)
    if m.sum() < 5:
        return None
    cv_m = cv[m]
    kcal = (F[m] - F[m].min()) * KJ_TO_KCAL

    k = min(5, (len(kcal) // 2) * 2 - 1)
    smooth = np.convolve(kcal, np.ones(k) / k, mode="same") if k >= 3 else kcal
    order = max(3, len(smooth) // 40)
    idx = argrelextrema(smooth, np.less_equal, order=order)[0]
    wells = sorted({(round(float(cv_m[i]), 2), round(float(kcal[i]), 2)) for i in idx},
                   key=lambda w: w[1])[:4]

    if len(wells) < 2:
        return float(kcal.max()), wells  # only one well found; no saddle to bound yet

    # station A = global minimum; station B = lowest-energy well on the
    # opposite side of A (the "other" shuttle station, not a shoulder of A's
    # own well) -- barrier is the max F strictly between them.
    a_cv, a_k = wells[0]
    opposite = [w for w in wells if (w[0] > 0) != (a_cv > 0)]
    if not opposite:
        return float(kcal.max()), wells  # no opposite-side well yet; can't bound a saddle
    b_cv, b_k = min(opposite, key=lambda w: w[1])
    lo, hi = sorted((a_cv, b_cv))
    win = (cv_m >= lo) & (cv_m <= hi)
    barrier = float(kcal[win].max() - min(a_k, b_k)) if win.any() else float(kcal.max())
    return barrier, wells


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hills", type=Path, default=Path("metad_run/HILLS"))
    ap.add_argument("--colvar", type=Path, default=Path("metad_run/COLVAR"))
    ap.add_argument("--grid-min", type=float, default=-31.5)
    ap.add_argument("--grid-max", type=float, default=31.5)
    ap.add_argument("--grid-bin", type=int, default=360)
    args = ap.parse_args()

    if not args.hills.exists() or args.hills.stat().st_size == 0:
        print("[checkin] no HILLS yet")
        return

    header, first, second = split_hills(args.hills)
    n_hills = len(first) + len(second)
    if n_hills < 20:
        print(f"[checkin] only {n_hills} hills deposited so far -- too early for a meaningful FES")
        return

    colvar = np.loadtxt(args.colvar, comments="#")
    cv_col = colvar[:, 1]
    obs_min, obs_max = float(cv_col.min()), float(cv_col.max())

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        header_text = "".join(header)
        results = {}
        for tag, lines in [("all", first + second), ("first_half", first), ("second_half", second)]:
            fes = sum_hills(header_text + "".join(lines), args.grid_min, args.grid_max,
                             args.grid_bin, tmp, tag)
            results[tag] = analyze(fes, obs_min, obs_max) if fes else None

    print(f"[checkin] {n_hills} hills deposited; sampled cv range {obs_min:.2f}..{obs_max:.2f} A")
    for tag, label in [("all", "ALL HILLS"), ("first_half", "FIRST HALF"), ("second_half", "SECOND HALF")]:
        r = results[tag]
        if r is None:
            print(f"  {label}: n/a")
            continue
        barrier, wells = r
        wells_str = ", ".join(f"{c:+.2f}A({k:.1f}kcal)" for c, k in wells) or "none found"
        print(f"  {label:12s}: barrier={barrier:5.2f} kcal/mol   wells: {wells_str}")


if __name__ == "__main__":
    main()
