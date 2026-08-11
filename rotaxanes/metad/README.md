# Well-tempered metadynamics of the rotaxane shuttle (PLUMED + OpenMM)

Drives the wheel along the rod with WT-MetaD to reconstruct the shuttling
free-energy profile F(d) — the solution-phase analog of the gas-phase GFN2-xTB
scan in the Rotaxanes repo (`displace_wheel.py` + `vib_stations.py`), but here
in explicit TIP3P solvent at 300 K, sampling the *real* shuttling coordinate
with the wheel free to rotate/pucker.

## Collective variable

`d` = signed distance (Å) of the wheel-O centroid from the midpoint of the two
rod amide nitrogens, projected onto the N→N axis:

```
d = ( W_centroid − (N1+N2)/2 ) · (N2 − N1) / |N2 − N1|
```

This is an **internal coordinate** (two rod atoms + the wheel oxygens) — invariant
to whole-molecule rotation and PBC wrapping, so it survives the tumbling and
image-wrapping that split a Cartesian "wheel x-coordinate" CV. `make_plumed.py`
reads `outputs/complex.pdb`, finds the rod `N` atoms and the wheel `O` atoms, and
emits `plumed.dat` with the right 1-based PLUMED indices.

## Files

| file | role |
|---|---|
| `make_plumed.py` | topology → `plumed.dat` (CV + WT-METAD + COLVAR logging) |
| `run_metad.py` | OpenMM driver: unbiased NPT equil → biased WT-METAD production |
| `run_metad.sh` | detached launcher (caffeinate + start_new_session) |
| `run_metad_resume.py` | warm-restart after a crash: last DCD frame + `RESTART` (HILLS appended) |
| `run_resume.sh` | detached launcher for the resume |
| `analyze_metad.py` | regenerate the deliverables (`Fes.dat`, `Fes_plot.png`, `pymol_rotaxane.*`) |
| `plumed.dat` | generated PLUMED input (regenerate after rebuilding the topology) |

## One-time setup: install the OpenMM PLUMED plugin

Already done (Aug 2026) — `openmm-plumed 2.1` is in the `openmm-md` env. Kept here
for reference / a fresh-machine rebuild:

```sh
source ~/miniforge3/etc/profile.d/conda.sh
conda activate openmm-md
conda install -c conda-forge openmm-plumed      # py311 build (2.1); downgrades openmm 8.5 -> 8.4
python -c "from openmmplumed import PlumedForce; print('plumed plugin ok')"
```

## Run

```sh
# 1. (re)build the topology -> plumed.dat (run once after omd build-multimol)
python make_plumed.py --topology ../outputs/complex.pdb --out plumed.dat

# 2. launch the WT-METAD run detached (it is a multi-hour job)
./run_metad.sh
```

`HILLS` (bias history) and `COLVAR` (d vs time) land in `metad_run/`. The
reweighted free energy F(d) is reconstructed from the HILLS file — e.g. with
`plumed sum_hills --hills metad_run/HILLS --mintozero --outfile Fes.dat`.

## Defaults / tuning

`run_metad.py` defaults: 5000-step unbiased NPT equilibration, then 5 M steps
(10 ns) biased production, auto platform (OpenCL single on this Mac). WT-METAD
defaults (in `plumed.dat`): `SIGMA=0.5 Å`, `HEIGHT=2.5 kJ/mol` (~kBT),
`PACE=500` steps, `BIASFACTOR=15`, `TEMP=300 K`.

To converge a ~14 kcal/mol shuttling barrier expect **10s of ns** of biased
sampling (the GFN2 scan found ΔG‡ ≈ 14.4 kcal on rot2htp). For a first look,
shorten with `--steps 2500000` (5 ns) and inspect `COLVAR` for barrier
crossings; extend if it hasn't crossed. To bias harder/faster: raise `--height`
or lower `--biasfactor`; to resolve the wells better: lower `--sigma`.

## Notes / caveats

- The CV uses the rod's two amide N atoms as the axis. They sit in the central
  diamide (~7.5 Å apart) but define the rod long-axis *direction* correctly; the
  midpoint origin keeps `d` symmetric about the rod centre.
- `PlumedForce` adds the bias as a regular force term; it runs on OpenCL single
  precision on this Mac. The bias force is computed in double inside PLUMED, so
  single-precision OpenMM integration is fine.
- The unbiased equilibration phase (no `PlumedForce`) settles the box before any
  hill is deposited, so the metad doesn't fight the solvent relaxation.
- `HILLS`/`COLVAR` paths are rewritten to absolute inside `metad_run/` by
  `run_metad.py`, so a run launched from any cwd lands them with the trajectory.

---

## Run record — Aug 2026 (25M-step rot2htpuma shuttle)

The first full WT-METAD run on `rot2htpuma` (rod + 24-crown-8 wheel, 114 solute
atoms + ~1008 waters, padding 1.2 nm). Target 25M steps (50 ns at dt = 2 fs); it
did not finish in one go — the OpenCL platform silently died twice (see
gotchas). The final F(d) is reconstructed from the 23.63M hills actually
deposited (94.5% of the run), which is effectively converged.

### Process & commands

```sh
source ~/miniforge3/etc/profile.d/conda.sh && conda activate openmm-md
cd rotaxanes/metad

# 0. topology already built via `omd build-multimol` (../outputs/complex.pdb,
#    system.xml). Build the PLUMED input from it:
python make_plumed.py --topology ../outputs/complex.pdb --out plumed.dat

# 1. launch the 25M-step WT-METAD run detached (multi-hour job):
./run_metad.sh                # -> metad_run/{traj.dcd, energy.csv, HILLS, COLVAR}

# 2. the original run died at 11.944M steps (47.8%) with no checkpoint. Resume:
./run_resume.sh               # warm restart from last traj.dcd frame + RESTART
                              # -> metad_run/{traj_resume.dcd, energy_resume.csv,
                              #                checkpoint_resume.bin}; HILLS appended

# 3. regenerate the deliverables (FES + PyMOL trajectory):
python analyze_metad.py       # -> Fes.dat, Fes_plot.png, pymol_rotaxane.{pdb,dcd}
```

### Results

| quantity | value |
|---|---|
| shuttle barrier ΔF‡ | **10.7 kcal/mol** (max−min of F(d) over \|d\| ≤ 6.5 Å) |
| terminal wells | **d ≈ ±4.6 Å**, degenerate: F = 0.0 and 0.5 kcal/mol |
| central metastable well | **d = 0**, ~4.9 kcal/mol above the terminals |
| CV range sampled | d = −7.86 … +7.86 Å (the wheel shuttles between the two stations, just past each amide N) |
| rod rigidity | N···N axis 7.30 ± 0.19 Å (no bending) |
| wheel escape | **none** — 3D min-image wheel–rod distance 0.04–7.48 Å, never near L/2 = 15.84 Å |

Physical picture: a symmetric double-well along the rod (two degenerate shuttle
stations at ±4.6 Å) with a shallow central metastable well at the rod midpoint
and a ~10.7 kcal/mol barrier between them — the solution-phase analog of the
gas-phase GFN2-xTB scan (ΔG‡ ≈ 14.4 kcal), lowered by solvent as expected.

### Where to find the deliverables (all in `metad_run/`, gitignored)

| file | what |
|---|---|
| `metad_run/Fes.dat` | the free-energy surface F(d): col 0 = d (Å), col 1 = F (kJ/mol, min = 0), col 2 = derivative |
| `metad_run/Fes_plot.png` | F(d) plotted in kcal/mol with the wells + barrier marked — **the "potential surface across the rod"** |
| `metad_run/pymol_rotaxane.pdb` | solute-only topology (114 atoms, first frame) for PyMOL |
| `metad_run/pymol_rotaxane.dcd` | solute-only trajectory, 23,631 frames, rod centred, wheel at its minimum-image position (no PBC snapping) |
| `metad_run/HILLS` | full bias history (47,264 hills; original + resume appended) — the source for Fes.dat |
| `metad_run/COLVAR` | d vs time (col 0 = time fs, col 1 = d Å) — the CV trace |
| `metad_run/traj.dcd` + `traj_resume.dcd` | raw wrapped solute+solvent trajectories (original + resume segments) |
| `metad_run/energy.csv` + `energy_resume.csv` | thermodynamic energy / T / volume vs time |
| `metad_run/energy_plot.png` | thermodynamic plot (E/T/V vs time) — NOT the FES; the FES is `Fes_plot.png` |

Load the structures in PyMOL:

```
load metad_run/pymol_rotaxane.pdb
load_traj metad_run/pymol_rotaxane.dcd
```

### Gotchas (each cost real time on this run)

- **OpenMM OpenCL silently dies on multi-hour runs on Apple Silicon.** Both the
  original (47.8%) and the detached resume (94.5%) crashed with no Python
  traceback, no OCL error, no `pmset` sleep event, and (for the resume) no
  harness reap (it was launched detached via `caffeinate -w $$` +
  `Popen(start_new_session=True)`). The platform is deprecated/fragile here.
  Mitigation: run metad on **CPU** (~2.9 ms/step vs 0.9 OpenCL, but stable), or
  keep OpenCL and checkpoint frequently (`run_metad_resume.py` writes
  `checkpoint_resume.bin` every 50k steps → a `loadCheckpoint` restart loses
  almost nothing). See `memory/rotaxane-md-plumed-metad.md`.
- **The wheel does NOT escape the rod — do not add walls.** `plumed.dat` uses
  plain `DISTANCE … COMPONENTS` (no `NOPBC`), so the CV is **minimum-image**
  (bounded by L/2 ≈ 15.8 Å) and only ever reaches ±7.86 Å; the measured 3D
  min-image wheel–rod distance never exceeds 7.48 Å. The "wheel escapes to
  ~48 Å, add LOWER/UPPER_WALLS" idea was a mis-unwrapping artifact, not physics.
- **The CV is minimum-image, so the FES grid `±12.77` is wider than needed.**
  The wheel only reaches ±7.86 Å; `GRID_MIN/MAX = ±12.77, GRID_BIN = 360` is
  fine (no escape to clamp), just wider than the sampled range. `sum_hills`
  auto-detects tighter boundaries; `analyze_metad.py` pins the grid to ±12.77 /
  360 for run-to-run reproducibility.
- **Two naive ways to "wrap" the PyMOL trajectory are both wrong:**
  1. Per-frame minimum-image on the naive COM *snaps* — the wrapped DCD splits
     molecules across PBC faces, so the COM is garbage mid-split, and min-image
     then teleports the wheel to the opposite side when \|rel\| crosses L/2. This
     is the "wheel constantly disappearing to the other side" failure.
  2. Per-atom continuity unwrapping *fakes a drift* — it accumulates the bulk
     diffusion of the rod and the wheel as if they separate, manufacturing a
     spurious ±60 Å wander.
  The **correct wrap** (in `analyze_metad.py`) is: make each molecule whole
  (reference every atom to its molecule's first atom via minimum-image — exact,
  since both molecules < L/2), then minimum-image the wheel–rod vector. Smooth
  because the wheel never crosses half a box; frame jump drops 62.99 Å → 1.96 Å
  (real 2-ps motion, no teleport); all frames kept.
- **COLVAR / energy.csv time is in fs, not steps.** dt = 2 fs → steps = time_fs / 2.
  The resume segment's step counter restarts at 0, so to get the absolute step
  add the original-run length (11,944,000) to the resume steps. `analyze_metad.py`
  joins `traj.dcd` + `traj_resume.dcd` directly so this bookkeeping is avoided.
- **`energy_plot.png` is not the FES.** It is total/potential/kinetic energy,
  temperature and box volume vs time — useful for sanity-checking the
  thermostat/barostat, but it shows no wells. The free-energy surface along the
  rod is `Fes_plot.png` (from `Fes.dat`).
- **`UNITS LENGTH=A` in `plumed.dat`** — PLUMED writes the CV in Å (COLVAR col 1
  is Å, not nm). Easy to misread as nm and conclude the wheel wanders ±78 Å.