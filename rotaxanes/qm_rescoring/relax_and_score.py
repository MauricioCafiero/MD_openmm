#!/usr/bin/env python
"""Loose local relaxation + GFN2-xTB single point on the extracted MD frames.

The raw single points in score_tblite.py turned out to carry a lot of
instantaneous MD strain (d=-4.90 came out +43 kcal/mol above the minimum --
way more than the MD/GAFF2 surface's 7.2 kcal/mol, and more than the gas-phase
pipeline's own careful numbers), because a single dynamics frame isn't a
relaxed structure -- ring pucker, bond angles etc. are all mid-motion.

Fix: pin the minimum needed to keep the reaction coordinate meaningful, then
let a loose LBFGS pass shake out the local steric/strain noise:
  - ONE rod atom (the most -u-extreme rod atom on our existing N1->N2 CV axis)
    -- NOT both rod tips. Pinning both tips (~30 A apart on rot1's long rod)
    is exactly the over-constraint that strains the central region and
    artificially lifts the central well -- the concern raised about the
    Rotaxanes repo's own naive `run_scan`. One point removes that risk while
    still preventing gross rigid-body drift during a *loose*, few-step pass.
  - ONE wheel unit (3 heavy atoms of one -O-CH2-CH2- triplet) -- the same
    `pinu` trick displace_wheel.py uses: 3 non-collinear pins fix the wheel's
    6 rigid-body DOF (can't translate/rotate as a whole) while the rest of the
    56-atom ring (and all H's) is free to pucker and relieve threading strain.

Not a production-quality relaxed scan (no stopper-walking, no symmetric
mirroring, no careful anchor-mode selection) -- just enough to get each
MD-sampled station off its instantaneous strain before scoring, cheaply.

Frame averaging: a dense scan at 0.5 A spacing showed >10 kcal/mol swings
between ADJACENT grid points -- not physical at that resolution. Cause: each
point was still a single MD snapshot, so snapshot-to-snapshot conformational
noise (ring pucker, local rotamer state) dominated over the underlying trend.
extract_stations.py's --n-frames now pulls several MD frames per target d
(spread across the trajectory, not adjacent in time, so they sample different
excursions rather than one correlated stretch); this script groups files by
their d-target (the `d{target}_r{rep}.xyz` naming) and averages the relaxed
energies per group, reporting the std-dev so the residual noise is visible
rather than hidden.

Usage:
  .venv/bin/python relax_and_score.py                  # gas phase
  .venv/bin/python relax_and_score.py --solvent water   # ALPB water
"""
from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

N1, N2 = 27, 46
WHL_UNIT = [88, 89, 90]  # first -O-CH2-CH2- triplet (non-collinear, heavy)
EV_TO_KCAL = 23.060548


def compute_d(pos: np.ndarray) -> float:
    n1, n2 = pos[N1], pos[N2]
    mid = (n1 + n2) / 2.0
    whl = pos[[88, 91, 94, 97, 100, 103, 106, 109]].mean(axis=0)
    axis = n2 - n1
    axis_hat = axis / np.linalg.norm(axis)
    return float((whl - mid) @ axis_hat)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frames-dir", type=Path, default=Path("frames"))
    ap.add_argument("--method", default="GFN2-xTB")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--multiplicity", type=int, default=1)
    ap.add_argument("--solvent", default=None,
                    help="e.g. 'water' -> ALPB implicit solvation; default gas phase")
    ap.add_argument("--fmax", type=float, default=0.15, help="loose eV/A tolerance")
    ap.add_argument("--steps", type=int, default=100)
    args = ap.parse_args()

    from ase.io import read
    from ase.constraints import FixAtoms
    from ase.optimize import LBFGS
    from tblite.ase import TBLite

    frames = sorted(args.frames_dir.glob("d*.xyz"))
    if not frames:
        raise SystemExit(f"no d*.xyz frames in {args.frames_dir} -- run extract_stations.py first")

    kwargs = dict(method=args.method, charge=args.charge,
                 multiplicity=args.multiplicity, verbosity=0)
    if args.solvent:
        kwargs["solvation"] = ("alpb", args.solvent)

    results = []
    for f in frames:
        atoms = read(f)
        pos0 = atoms.get_positions()
        rod_anchor = int(np.argmin(
            (pos0[:88] - pos0[N1]) @ ((pos0[N2] - pos0[N1]) / np.linalg.norm(pos0[N2] - pos0[N1]))
        ))
        atoms.set_constraint(FixAtoms([rod_anchor] + WHL_UNIT))
        atoms.calc = TBLite(**kwargs)

        d_before = compute_d(pos0)
        opt = LBFGS(atoms, logfile=None)
        opt.run(fmax=args.fmax, steps=args.steps)
        e_ev = atoms.get_potential_energy()
        d_after = compute_d(atoms.get_positions())

        m = re.match(r"d([+-][\d.]+)(?:_r\d+)?\.xyz", f.name)
        d_target = float(m.group(1))
        results.append((d_target, e_ev, d_before, d_after, opt.get_number_of_steps(),
                        bool(opt.converged())))
        print(f"  {f.name}: rod_anchor={rod_anchor}  d {d_before:+.2f}->{d_after:+.2f} A  "
              f"E={e_ev:.6f} eV  ({opt.get_number_of_steps()} steps, "
              f"{'converged' if opt.converged() else 'step-capped'})")

    # group replicates by target d and average
    groups: dict[float, list[float]] = defaultdict(list)
    for d_target, e_ev, *_ in results:
        groups[d_target].append(e_ev)
    grouped = sorted((d_target, np.mean(es), np.std(es), len(es))
                     for d_target, es in groups.items())

    e_min = min(mean for _, mean, _, _ in grouped)
    tag = f"ALPB/{args.solvent}" if args.solvent else "gas-phase"
    print(f"\n{'d target':>9}  {'n':>3}  {'E_rel mean (kcal/mol, ' + tag + ')':>34}  {'std':>8}")
    for d_target, e_mean, e_std, n in grouped:
        rel = (e_mean - e_min) * EV_TO_KCAL
        std_kcal = e_std * EV_TO_KCAL
        print(f"{d_target:9.2f}  {n:3d}  {rel:34.2f}  {std_kcal:8.2f}")


if __name__ == "__main__":
    main()
