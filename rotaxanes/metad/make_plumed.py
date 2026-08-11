#!/usr/bin/env python
"""Generate a PLUMED input for well-tempered metadynamics on the rotaxane
shuttle coordinate, with atom indices read from the built complex.pdb.

Shuttle CV (internal coordinate -- invariant to whole-molecule rotation / PBC):
    cv = ( W - midpoint(N1,N2) ) . (N2 - N1) / |N2 - N1|
where N1, N2 are the two rod amide nitrogens (the central diamide, ~7.5 A apart,
spanning the rod's long axis) and W is the centroid of the wheel's 8 crown-ether
oxygens. The midpoint origin makes cv symmetric about the rod centre, so the two
shuttle wells land at ~ +/- d. For rot2htpuma the range is ~ +/- 7.5 A.

The WT-METAD deposits Gaussian hills on cv to flatten the shuttling barrier;
PLUMED writes HILLS (the bias history) and COLVAR (cv vs time) for reweighting /
plotting. Lengths in Angstrom (matches the Rotaxanes repo's d-in-A convention).
"""
from __future__ import annotations

import argparse
from pathlib import Path

from openmm import app, unit
import numpy as np


def find_cv_atoms(pdb_path: Path) -> tuple[list[int], list[int]]:
    """1-based indices of the two rod N atoms and the wheel O atoms."""
    pdb = app.PDBFile(str(pdb_path))
    rod_n, wheel_o = [], []
    for a in pdb.topology.atoms():
        if a.residue.name == "ROD" and a.element and a.element.symbol == "N":
            rod_n.append(a.index + 1)
        elif a.residue.name == "WHL" and a.element and a.element.symbol == "O":
            wheel_o.append(a.index + 1)
    if len(rod_n) != 2:
        raise ValueError(f"expected 2 rod N atoms, found {len(rod_n)}: {rod_n}")
    if not wheel_o:
        raise ValueError("no wheel O atoms found (is there a WHL residue?)")
    return rod_n, wheel_o


def cv_range_ang(pdb_path: Path, rod_n: list[int], wheel_o: list[int]) -> float:
    """Half-range of the shuttle CV (A), from the built geometry: the wheel can
    travel from one rod tip to the other, ~ the rod length. We approximate the
    extent by 1.5x the N-N axis length + a margin (the stopper arms extend past
    the N's)."""
    pdb = app.PDBFile(str(pdb_path))
    pos = pdb.positions.value_in_unit(unit.angstrom)
    p = lambda i: np.array([float(pos[i - 1][k]) for k in range(3)])
    n1, n2 = rod_n
    d_nn = np.linalg.norm(p(n2) - p(n1))
    return float(d_nn * 1.5 + 1.5)  # A; grid half-width


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--topology", type=Path, default=Path("../outputs/complex.pdb"),
                   help="built complex.pdb (ROD + WHL residues).")
    p.add_argument("--out", type=Path, default=Path("plumed.dat"),
                   help="output PLUMED input file.")
    # WT-METAD parameters
    p.add_argument("--sigma", type=float, default=0.5, help="hill width in A (default 0.5)")
    p.add_argument("--height", type=float, default=2.5,
                   help="initial hill height in kJ/mol (default 2.5 ~ kBT)")
    p.add_argument("--pace", type=int, default=500, help="hill deposition stride in steps")
    p.add_argument("--biasfactor", type=float, default=15.0,
                   help="well-tempered biasfactor (default 15)")
    p.add_argument("--temp", type=float, default=300.0, help="temperature K (default 300)")
    p.add_argument("--grid-bin", type=int, default=360)
    return p.parse_args()


def main():
    args = parse_args()
    rod_n, wheel_o = find_cv_atoms(args.topology)
    half = cv_range_ang(args.topology, rod_n, wheel_o)
    n1, n2 = rod_n
    olist = ",".join(str(i) for i in wheel_o)

    plumed = f"""# === rot2htpuma shuttle -- well-tempered metadynamics ===
# CV: signed distance (A) of the wheel-O centroid from the midpoint of the two
# rod amide N's, projected onto the N->N axis. Internal coordinate -> invariant
# to whole-molecule rotation / PBC wrapping.
UNITS LENGTH=A TIME=fs ENERGY=kj/mol

# virtual atoms: rod-N midpoint + wheel-O centroid
MID: CENTER ATOMS={n1},{n2}
WHL: CENTER ATOMS={olist}

# (N2 - N1) axis vector and (WHL - MID) vector, as components
axis: DISTANCE ATOMS={n1},{n2} COMPONENTS
dvec: DISTANCE ATOMS=MID,WHL COMPONENTS

# cv = (dvec . axis) / |axis|   (6 args -> PLUMED needs explicit VAR= names)
cv: CUSTOM ARG=dvec.x,dvec.y,dvec.z,axis.x,axis.y,axis.z \
    VAR=dvx,dvy,dvz,ax,ay,az \
    FUNC=(dvx*ax+dvy*ay+dvz*az)/sqrt(ax*ax+ay*ay+az*az) PERIODIC=NO

# well-tempered metadynamics on the shuttle CV
METAD ...
  ARG=cv
  SIGMA={args.sigma}
  HEIGHT={args.height}
  PACE={args.pace}
  BIASFACTOR={args.biasfactor}
  TEMP={args.temp}
  GRID_MIN=-{half:.2f}
  GRID_MAX={half:.2f}
  GRID_BIN={args.grid_bin}
  FILE=HILLS
... METAD

# log the CV (and the bias) for reweighting / plotting
PRINT ARG=cv STRIDE=100 FILE=COLVAR
FLUSH STRIDE=100
"""
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(plumed)
    print(f"wrote {args.out}")
    print(f"  rod N atoms (1-based): {n1}, {n2}   | axis length ~ {half/1.5:.2f} A")
    print(f"  wheel O atoms ({len(wheel_o)}): {wheel_o}")
    print(f"  CV grid: -{half:.2f} .. +{half:.2f} A  ({args.grid_bin} bins)")
    print(f"  WT-METAD: SIGMA={args.sigma} A  HEIGHT={args.height} kJ/mol  "
          f"PACE={args.pace}  BIASFACTOR={args.biasfactor}  TEMP={args.temp} K")


if __name__ == "__main__":
    main()