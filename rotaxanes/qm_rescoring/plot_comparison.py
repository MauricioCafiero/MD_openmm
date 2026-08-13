#!/usr/bin/env python
"""Plot the QM-refined (GFN2-xTB, ALPB water) points from relax_and_score.py
against the MD/GAFF2 Fes.dat curve, so the comparison (and the point-to-point
noise) is actually visible instead of just tabulated.

Usage:
  .venv/bin/python plot_comparison.py                       # reads relax_dense.log
  .venv/bin/python plot_comparison.py --log relax_dense.log --out qm_vs_md.png
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np


def parse_log(log_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Pull the (d target, E_rel mean kcal/mol, std) table rows out of a
    relax_and_score.py log -- the grouped-by-target summary table
    ('d target', 'n', 'E_rel mean...', 'std' columns)."""
    d, e, std = [], [], []
    in_table = False
    for line in log_path.read_text().splitlines():
        if line.strip().startswith("d target"):
            in_table = True
            continue
        if not in_table or not line.strip():
            continue
        parts = line.split()
        if len(parts) != 4:
            continue
        try:
            d.append(float(parts[0]))
            e.append(float(parts[2]))
            std.append(float(parts[3]))
        except ValueError:
            continue
    return np.array(d), np.array(e), np.array(std)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--log", type=Path, default=Path("relax_dense.log"))
    ap.add_argument("--fes", type=Path, default=Path("../metad/metad_run/Fes.dat"))
    ap.add_argument("--colvar", type=Path, default=Path("../metad/metad_run/COLVAR"))
    ap.add_argument("--out", type=Path, default=Path("qm_vs_md.png"))
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    qm_d, qm_e, qm_std = parse_log(args.log)
    print(f"parsed {len(qm_d)} QM points from {args.log} "
          f"(mean std across replicates: {qm_std.mean():.2f} kcal/mol)")

    fes = np.loadtxt(args.fes, comments="#")
    colvar = np.loadtxt(args.colvar, comments="#")
    obs_min, obs_max = float(colvar[:, 1].min()), float(colvar[:, 1].max())
    cv, F = fes[:, 0], fes[:, 1]
    m = (cv >= obs_min) & (cv <= obs_max) & ~np.isnan(F)
    md_d, md_kcal = cv[m], (F[m] - F[m].min()) / 4.184

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(md_d, md_kcal, lw=2, color="#1f5fa8", label="MD/GAFF2 (Fes.dat)")
    ax.errorbar(qm_d, qm_e, yerr=qm_std, fmt="o-", lw=1, ms=4, capsize=2,
               color="#c0392b", alpha=0.85,
               label="QM-refined (GFN2-xTB, ALPB water, frame-averaged)")
    ax.set_xlabel("shuttle coordinate  d  (Å)")
    ax.set_ylabel("F or E_rel  (kcal/mol, each curve zeroed to its own minimum)")
    ax.set_title("rot1 shuttle: MD/GAFF2 vs QM-refined single-snapshot rescoring")
    ax.legend(loc="upper center")
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(args.out, dpi=150)
    print(f"-> {args.out}")


if __name__ == "__main__":
    main()
