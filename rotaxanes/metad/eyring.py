#!/usr/bin/env python
"""Eyring/TST rate estimate from a free-energy barrier, matching the exact
convention the Rotaxanes repo uses in rot1_shuttles.md:

    k = (kB*T / h) * exp(-DeltaG_double_dagger / (R*T))

with the TST prefactor kB*T/h ~ 6.21e12 s^-1 at T = 300 K (their rounding;
computed exactly here from CODATA constants, which comes out to ~6.25e12 --
close enough that neither number changes the reported rate/timescale at the
precision these barriers are known to). tau = 1/k is the mean hop time.

These are estimates: a constrained-surface (their case) or classical-MD
(ours) Delta-G barrier plugged into a bare TST prefactor, no
recrossing/friction (Kramers) correction, no transmission coefficient --
same caveat the Rotaxanes repo states for its own numbers.

Usage:
  python eyring.py 12.59                  # single barrier, kcal/mol
  python eyring.py 12.59 9.2 14.2 10.7 9.51 14.25 --labels ...
"""
from __future__ import annotations

import argparse

KB = 1.380649e-23      # J/K
H = 6.62607015e-34     # J*s
R = 8.314462618        # J/(mol*K)
KCAL_TO_J = 4184.0


def rate(barrier_kcal: float, T: float = 300.0) -> tuple[float, float]:
    """(k in s^-1, tau in s) for a barrier in kcal/mol at temperature T (K)."""
    dG = barrier_kcal * KCAL_TO_J          # J/mol
    prefactor = KB * T / H                  # s^-1
    k = prefactor * pow(2.718281828459045, -dG / (R * T))
    return k, 1.0 / k


def fmt_time(tau_s: float) -> str:
    if tau_s < 1e-6:
        return f"{tau_s*1e9:.2f} ns"
    if tau_s < 1e-3:
        return f"{tau_s*1e6:.2f} us"
    if tau_s < 1:
        return f"{tau_s*1e3:.2f} ms"
    if tau_s < 60:
        return f"{tau_s:.2f} s"
    if tau_s < 3600:
        return f"{tau_s/60:.2f} min"
    if tau_s < 86400:
        return f"{tau_s/3600:.2f} h"
    return f"{tau_s/86400:.2f} d"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("barriers", type=float, nargs="+", help="barrier(s) in kcal/mol")
    ap.add_argument("--labels", nargs="*", default=None)
    ap.add_argument("--temp", type=float, default=300.0, help="K (default 300, matches the MD)")
    args = ap.parse_args()

    labels = args.labels or [f"{b:.2f} kcal/mol" for b in args.barriers]
    if len(labels) != len(args.barriers):
        labels = [f"{b:.2f} kcal/mol" for b in args.barriers]

    print(f"T = {args.temp} K   kB*T/h = {KB*args.temp/H:.3e} s^-1\n")
    print(f"{'label':38s} {'DeltaG (kcal/mol)':>18s} {'k (s^-1)':>12s} {'tau':>12s}")
    for label, b in zip(labels, args.barriers):
        k, tau = rate(b, args.temp)
        print(f"{label:38s} {b:18.2f} {k:12.3e} {fmt_time(tau):>12s}")


if __name__ == "__main__":
    main()
