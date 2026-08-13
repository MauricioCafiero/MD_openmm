#!/usr/bin/env python
"""GFN2-xTB single-point rescoring of the frames extracted by
extract_stations.py, using the tblite ASE calculator -- same engine as the
Rotaxanes repo's gas-phase pipeline (code/vib_stations.py,
code/displace_wheel.py --engine tblite).

Purpose: test whether GAFF2 (the classical force field driving our
explicit-solvent MetaD) and GFN2-xTB (tight-binding QM) agree on the RELATIVE
energetics across the wells/saddle rot1's MetaD found -- particularly whether
the solution-phase well-depth reordering we saw (stopper-region wells deeper
than central wells, opposite of the gas-phase pipeline's finding) survives a
QM re-scoring of the same geometries, or is a GAFF2 artifact. See PLAN.md.

These are gas-phase single points of solution-phase-sampled geometries (no
implicit/explicit solvent here) -- a first-pass consistency check, not a
solvent-corrected QM free energy. PLAN.md's step 2/3 are the follow-ups.

Requires: ase, tblite (with the tblite.ase ASE calculator) -- NOT in the
openmm-md conda env; see PLAN.md for install notes. Does NOT need
fairchem-core / torch / HF_TOKEN (those are only for the UMA engine, unused
here, and per the Rotaxanes repo's CLAUDE.md cannot share a process with
tblite anyway -- both bundle their own libomp and segfault if combined).

Usage:
  python score_tblite.py                     # scores everything in frames/
  python score_tblite.py --frames-dir frames --method GFN2-xTB
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

EV_TO_KCAL = 23.060548


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frames-dir", type=Path, default=Path("frames"))
    ap.add_argument("--method", default="GFN2-xTB",
                    help="tblite method (default GFN2-xTB; also GFN1-xTB, GFN0-xTB, CEH)")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--multiplicity", type=int, default=1)
    args = ap.parse_args()

    from ase.io import read
    from tblite.ase import TBLite

    frames = sorted(args.frames_dir.glob("d*.xyz"))
    if not frames:
        raise SystemExit(f"no d*.xyz frames in {args.frames_dir} -- run extract_stations.py first")

    results = []
    for f in frames:
        atoms = read(f)
        atoms.calc = TBLite(method=args.method, charge=args.charge,
                            multiplicity=args.multiplicity, verbosity=0)
        e_ev = atoms.get_potential_energy()
        m = re.match(r"d([+-][\d.]+)\.xyz", f.name)
        d_target = float(m.group(1))
        results.append((d_target, e_ev))
        print(f"  {f.name}: E = {e_ev:.6f} eV")

    results.sort(key=lambda r: r[0])
    e_min = min(e for _, e in results)
    print(f"\n{'d (A)':>8}  {'E_rel (kcal/mol, GFN2-xTB, gas-phase)':>38}")
    for d_target, e_ev in results:
        rel = (e_ev - e_min) * EV_TO_KCAL
        print(f"{d_target:8.2f}  {rel:38.2f}")
    print(f"\nCompare against the MetaD/GAFF2 relative free energies in "
          f"../metad/metad_run/Fes.dat (wells: -10.85(0.0), +11.20(2.0), "
          f"-4.90(7.2), +4.20(7.8) kcal/mol) -- see PLAN.md for how to read this.")


if __name__ == "__main__":
    main()
