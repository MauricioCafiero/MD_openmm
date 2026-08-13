#!/usr/bin/env python
"""Regenerate the deliverables for a finished rotaxane WT-METAD run.

Run from ``rotaxanes/metad/`` after the metad run is done (i.e. after
``metad_run/HILLS``, ``metad_run/traj.dcd`` and -- if it was resumed --
``metad_run/traj_resume.dcd`` exist). Produces, in ``metad_run/``:

  * ``Fes.dat``         -- free-energy surface F(d) from the full HILLS history
                           (shells out to ``plumed sum_hills --mintozero``).
  * ``Fes_plot.png``    -- F(d) plotted in kcal/mol, with auto-detected wells
                           (local minima within the sampled cv range) and the
                           shuttle barrier marked.
  * ``pymol_rotaxane.pdb`` + ``pymol_rotaxane.dcd`` -- solute-only trajectory
                           with the rod centred at the origin and the wheel
                           placed at its minimum-image position. All frames
                           are kept and there is NO PBC snapping, because the
                           wheel never leaves half a box of the rod (see README
                           "gotchas" for why the two naive wraps are wrong).

Atom layout (ROD/WHL atom-index ranges, the two rod amide N's, the wheel
oxygens) is derived from the built ``complex.pdb`` at runtime -- the same
approach ``make_plumed.py`` uses -- rather than hardcoded, so this works for
any rod/wheel pair (rot2htpuma's 58+56=114 atoms, rot1's 88+56=144, etc.)
without a molecule-specific copy of this script. Likewise the FES grid
(GRID_MIN/MAX/BIN) is read from the run's ``plumed.dat`` so the grid used here
always matches the grid the METAD bias was actually deposited on.

Usage:
  python analyze_metad.py                 # do all three steps
  python analyze_metad.py --no-fes        # skip sum_hills + plot (PyMOL only)
  python analyze_metad.py --no-pymol      # skip the trajectory wrap (FES only)
  python analyze_metad.py --topology ../../outputs/complex.pdb --plumed plumed.dat --out-dir metad_run

Requires the ``openmm-md`` env (openmm, mdtraj, numpy, scipy, matplotlib, and
the ``plumed`` binary on PATH). Activate it first:
  source ~/miniforge3/etc/profile.d/conda.sh && conda activate openmm-md
"""
from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

import numpy as np
from scipy.signal import argrelextrema

KJ_TO_KCAL = 1.0 / 4.184


def atom_layout(topology_path: Path) -> tuple[slice, slice, int, int, list[int]]:
    """(rod_slice, whl_slice, n1, n2, whl_o) -- all 0-based -- from the built
    complex.pdb: ROD/WHL atom-index ranges (assumed contiguous per residue, as
    build_rotaxane.py writes them), the two rod amide N's, and the wheel O's.
    Mirrors make_plumed.py's find_cv_atoms (which finds the same atoms in
    1-based PLUMED indices)."""
    from openmm import app

    pdb = app.PDBFile(str(topology_path))
    rod_idx, whl_idx, rod_n, whl_o = [], [], [], []
    for a in pdb.topology.atoms():
        if a.residue.name == "ROD":
            rod_idx.append(a.index)
            if a.element and a.element.symbol == "N":
                rod_n.append(a.index)
        elif a.residue.name == "WHL":
            whl_idx.append(a.index)
            if a.element and a.element.symbol == "O":
                whl_o.append(a.index)
    if len(rod_n) != 2:
        raise ValueError(f"expected 2 rod N atoms, found {len(rod_n)}: {rod_n}")
    if not whl_o:
        raise ValueError("no wheel O atoms found (is there a WHL residue?)")
    rod_slice = slice(min(rod_idx), max(rod_idx) + 1)
    whl_slice = slice(min(whl_idx), max(whl_idx) + 1)
    return rod_slice, whl_slice, rod_n[0], rod_n[1], whl_o


def parse_grid(plumed_path: Path) -> tuple[float, float, int]:
    """Read GRID_MIN/GRID_MAX/GRID_BIN out of the run's plumed.dat, so the FES
    here is always binned on the same grid the bias was deposited on."""
    text = plumed_path.read_text()
    gmin = float(re.search(r"GRID_MIN=(-?[\d.]+)", text).group(1))
    gmax = float(re.search(r"GRID_MAX=(-?[\d.]+)", text).group(1))
    gbin = int(re.search(r"GRID_BIN=(\d+)", text).group(1))
    return gmin, gmax, gbin


def observed_cv_range(colvar_path: Path) -> tuple[float, float]:
    d = np.loadtxt(colvar_path, comments="#")
    return float(d[:, 1].min()), float(d[:, 1].max())


def make_fes(hills: Path, out: Path, grid_min: float, grid_max: float, grid_bin: int) -> Path:
    """F(d) from the full HILLS history via plumed sum_hills (--mintozero)."""
    cmd = ["plumed", "sum_hills", "--hills", str(hills), "--mintozero",
           "--outfile", str(out),
           "--min", str(grid_min), "--max", str(grid_max), "--bin", str(grid_bin)]
    print("[analyze] $ " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"[analyze] Fes.dat -> {out}")
    return out


def find_wells(cv: np.ndarray, kcal: np.ndarray, max_wells: int = 4) -> list[tuple[float, float]]:
    """Local minima of the (lightly smoothed) F(d), same method as
    checkin_fes.py's live snapshots, so the final plot's wells are found the
    same way the hourly check-ins found them."""
    k = min(5, (len(kcal) // 2) * 2 - 1)
    smooth = np.convolve(kcal, np.ones(k) / k, mode="same") if k >= 3 else kcal
    order = max(3, len(smooth) // 40)
    idx = argrelextrema(smooth, np.less_equal, order=order)[0]
    wells = sorted({(round(float(cv[i]), 2), round(float(kcal[i]), 2)) for i in idx},
                   key=lambda w: w[1])
    return wells[:max_wells]


def station_barrier(cv: np.ndarray, kcal: np.ndarray, wells: list[tuple[float, float]]) -> float:
    """Barrier restricted to the window BETWEEN the two station wells (global
    min + the lowest-energy well on the opposite side of the CV axis), not the
    whole sampled range. The tails beyond each station, out toward the sampled
    range's edges, are the least-visited/least-flattened parts of the surface
    and can carry a spuriously steep, unconverged rise unrelated to the real
    shuttling transition state -- caught in practice on rot1: the naive
    max-min over the full range picked up a peak at the extreme edge of
    sampling, past the deep well, inflating the reported barrier by ~2x."""
    if len(wells) < 2:
        return float(kcal.max())
    a_cv, a_k = wells[0]
    opposite = [w for w in wells if (w[0] > 0) != (a_cv > 0)]
    if not opposite:
        return float(kcal.max())
    b_cv, b_k = min(opposite, key=lambda w: w[1])
    lo, hi = sorted((a_cv, b_cv))
    win = (cv >= lo) & (cv <= hi)
    return float(kcal[win].max() - min(a_k, b_k)) if win.any() else float(kcal.max())


def plot_fes(fes: Path, out_png: Path, obs_min: float, obs_max: float) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = np.loadtxt(fes, comments="#")
    cv, F = d[:, 0], d[:, 1]
    F = F - np.nanmin(F)               # mintozero already does this, but be safe
    kcal = F * KJ_TO_KCAL
    ok = ~np.isnan(kcal)

    # wells restricted to the actually-sampled cv range -- unsampled tail bins
    # sit at the flat unbiased grid edge and would otherwise look like
    # spurious deep wells (same reasoning as checkin_fes.py). Barrier is then
    # further restricted to between the two station wells (see station_barrier).
    m = ok & (cv >= obs_min) & (cv <= obs_max)
    wells = find_wells(cv[m], kcal[m])
    barrier = station_barrier(cv[m], kcal[m], wells)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(cv[ok], kcal[ok], lw=2, color="#1f5fa8")
    ax.set_xlabel("shuttle coordinate  d  (Å)  [wheel-O centroid, projected on N···N axis]")
    ax.set_ylabel("F(d)  (kcal/mol)")
    ax.set_title("Rotaxane free-energy surface along the rod  (WT-MetaD)")
    ax.axhline(0, color="0.7", lw=0.8)
    ax.axvspan(obs_min, obs_max, color="#1f5fa8", alpha=0.06, label="sampled range")
    for c, k in wells:
        ax.plot(c, k, "o", color="#c0392b", ms=6)
        ax.annotate(f"{k:.1f}", (c, k), textcoords="offset points", xytext=(6, 6), fontsize=9)
    ax.text(0.02, 0.95, f"barrier ≈ {barrier:.1f} kcal/mol",
            transform=ax.transAxes, va="top", fontsize=10,
            bbox=dict(fc="white", ec="0.7", boxstyle="round"))
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    wells_str = ", ".join(f"{c:+.2f}A({k:.1f}kcal)" for c, k in wells)
    print(f"[analyze] Fes_plot.png -> {out_png}  (barrier ≈ {barrier:.1f} kcal/mol; wells: {wells_str})")
    return out_png


def _bond_adjacency(topology, m0: int, m1: int) -> dict[int, list[int]]:
    """Adjacency list (local 0-based indices within [m0, m1)) from the mdtraj
    topology's bonds, restricted to bonds internal to that atom range."""
    import collections
    adj: dict[int, list[int]] = collections.defaultdict(list)
    for b in topology.bonds:
        i, j = b[0].index, b[1].index
        if m0 <= i < m1 and m0 <= j < m1:
            adj[i - m0].append(j - m0)
            adj[j - m0].append(i - m0)
    return adj


def _unwrap_by_bonds(P: np.ndarray, L: np.ndarray, adj: dict[int, list[int]]) -> np.ndarray:
    """Make one molecule whole by walking its bond graph (BFS from atom 0) and
    unwrapping each atom relative to its already-unwrapped bonded parent via
    minimum-image. Bond lengths are always ~1-2 A, far below L/2 regardless of
    the molecule's total extent -- unlike referencing every atom to a single
    fixed atom, which silently breaks for an elongated molecule whose overall
    span exceeds L/2 (rot1's rod: 30.4 A end-to-end vs L/2 = 21.65 A here).
    ``P`` is (F, n, 3) in nm for just this molecule's atoms; ``L`` is (F, 3) nm."""
    import collections
    F, n, _ = P.shape
    W = P.copy()
    visited = np.zeros(n, dtype=bool)
    visited[0] = True
    order = [0]
    parent = {}
    q = collections.deque([0])
    while q:
        u = q.popleft()
        for v in adj.get(u, ()):
            if not visited[v]:
                visited[v] = True
                parent[v] = u
                order.append(v)
                q.append(v)
    if not visited.all():
        missing = np.where(~visited)[0].tolist()
        raise ValueError(f"atoms not reachable via bonds from atom 0: {missing} "
                          f"(disconnected fragment or missing CONECT records)")
    for v in order[1:]:
        u = parent[v]
        delta = W[:, v, :] - W[:, u, :]
        delta -= L * np.round(delta / L)
        W[:, v, :] = W[:, u, :] + delta
    return W


def _make_molecules_whole(P: np.ndarray, L: np.ndarray, rod: slice, whl: slice,
                          topology) -> np.ndarray:
    """Make both the rod and the wheel whole (no split across a PBC face) via
    bond-graph unwrapping. ``P`` is (F, N, 3) in nm; ``L`` is (F, 3) in nm."""
    W = P.copy()
    for m0, m1 in [(rod.start, rod.stop), (whl.start, whl.stop)]:
        adj = _bond_adjacency(topology, m0, m1)
        W[:, m0:m1, :] = _unwrap_by_bonds(P[:, m0:m1, :], L, {k: v for k, v in adj.items()})
    return W


def make_pymol(trajs, topology: Path, out_dir: Path, rod: slice, whl: slice) -> tuple[Path, Path]:
    import mdtraj as md

    t = trajs[0]
    for extra in trajs[1:]:
        t = t.join(extra)
    solute = t.atom_slice(range(rod.start, whl.stop))
    L = np.asarray(t.unitcell_lengths)              # (F, 3) nm

    # rod/whl slices are relative to the full topology; re-zero them to the
    # solute-only slice's atom indices (solute starts at rod.start).
    rod0 = slice(0, rod.stop - rod.start)
    whl0 = slice(rod0.stop, rod0.stop + (whl.stop - whl.start))

    W = _make_molecules_whole(solute.xyz.copy(), L, rod0, whl0, solute.topology)
    rodc = W[:, rod0, :].mean(1)                     # (F, 3) whole-molecule COM
    whlc = W[:, whl0, :].mean(1)
    rel = whlc - rodc
    rel -= L * np.round(rel / L)                     # minimum-image (PLUMED-style)

    P = W - rodc[:, None, :]                         # rod -> origin
    P[:, whl0, :] += (rel - (whlc - rodc))[:, None, :]  # wheel -> min-image rel (no snap)
    solute.xyz = P

    pdb = out_dir / "pymol_rotaxane.pdb"
    dcd = out_dir / "pymol_rotaxane.dcd"
    solute[0].save_pdb(str(pdb))
    solute.save_dcd(str(dcd))

    # report: wheel-rod distance + frame jump (sanity that there are no L-snaps)
    dvec = solute.xyz[:, whl0, :].mean(1) - solute.xyz[:, rod0, :].mean(1)
    dist = np.linalg.norm(dvec, axis=1) * 10.0
    jmp = np.linalg.norm(np.diff(dvec, axis=0), axis=1) * 10.0
    print(f"[analyze] pymol_traj -> {pdb}, {dcd}  "
          f"({solute.n_frames} frames; wheel-rod {dist.min():.2f}..{dist.max():.2f} A; "
          f"frame jump max {jmp.max():.2f} A)")
    return pdb, dcd


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--topology", type=Path, default=Path("../../outputs/complex.pdb"))
    ap.add_argument("--plumed", type=Path, default=Path("plumed.dat"))
    ap.add_argument("--out-dir", type=Path, default=Path("metad_run"))
    ap.add_argument("--no-fes", action="store_true", help="skip sum_hills + plot")
    ap.add_argument("--no-pymol", action="store_true", help="skip the trajectory wrap")
    args = ap.parse_args()

    out = args.out_dir
    rod, whl, n1, n2, whl_o = atom_layout(args.topology)
    print(f"[analyze] atom layout from {args.topology}: "
          f"ROD={rod.start}:{rod.stop}  WHL={whl.start}:{whl.stop}  "
          f"N1,N2={n1},{n2}  WHL_O={whl_o}")

    if not args.no_fes:
        gmin, gmax, gbin = parse_grid(args.plumed)
        obs_min, obs_max = observed_cv_range(out / "COLVAR")
        print(f"[analyze] grid from {args.plumed}: {gmin}..{gmax} ({gbin} bins); "
              f"sampled cv range {obs_min:.2f}..{obs_max:.2f} A")
        make_fes(out / "HILLS", out / "Fes.dat", gmin, gmax, gbin)
        plot_fes(out / "Fes.dat", out / "Fes_plot.png", obs_min, obs_max)
    if not args.no_pymol:
        import re
        import mdtraj as md
        trajs = [md.load_dcd(str(out / "traj.dcd"), top=str(args.topology))]
        # join every resume leg in order (traj_resume.dcd, traj_resume2.dcd, ...)
        # -- a run may have been extended more than once (see
        # run_metad_resume.py's latest_leg()), so don't stop after the first.
        def leg_idx(p):
            m = re.match(r"traj_resume(\d*)\.dcd", p.name)
            return int(m.group(1)) if m.group(1) else 1
        for resume in sorted(out.glob("traj_resume*.dcd"), key=leg_idx):
            trajs.append(md.load_dcd(str(resume), top=str(args.topology)))
        make_pymol(trajs, args.topology, out, rod, whl)


if __name__ == "__main__":
    main()
