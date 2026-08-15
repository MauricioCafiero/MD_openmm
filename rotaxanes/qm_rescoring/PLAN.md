# QM rescoring of the rot1 MetaD wells — plan

## Why

**Note:** the table below is the 10M-step snapshot this plan was originally
written against. The MetaD run was subsequently extended a third time (+5M,
15M steps total) and the barrier settled higher (~14.2 kcal/mol, see
`../metad/README.md`'s run record) with a bigger well-depth gap (5.1 kcal/mol
vs the 2.0 below) — the qualitative "terminal wells favored over central"
finding is unchanged, but if re-running any of this folder's scripts, pull
fresh target `d` values from the current `Fes.dat` rather than reusing the
numbers here verbatim.

The rot1 WT-MetaD run (`../metad/`, 10M steps at the time, corrected
`station_barrier()` metric) found a solution-phase free-energy surface with
wells at:

| d (Å) | role | F above min (kcal/mol) |
|------:|------|---:|
| −10.85 | deepest well (global min) | 0.0 |
| +11.20 | other terminal well | 2.0 |
| −4.90 | intermediate feature | 7.2 |
| +4.20 | intermediate feature | 7.8 |

The [mauriciocafiero/Rotaxanes](https://github.com/mauriciocafiero/Rotaxanes)
repo's **gas-phase** GFN2-xTB/UMA scan + constrained partial-Hessian
vibrational free-energy pipeline (`rot1_shuttles.md`) found the *same rod*
has four analogous stations, but with the **opposite relative depth**:

| d (Å), their frame | role | ΔG (kcal/mol) |
|------:|------|---:|
| ±4.80 | central wells | **0.0 (global min)** |
| ±9.6 | stopper wells | **+6.2** |

i.e. in gas phase the *central* wells are deepest and the *stopper* wells are
~6 kcal above; in our explicit-solvent MD the *stopper*-region wells are
deepest and the *central* wells sit ~7-8 kcal above (at the 10M-step
snapshot; the gap grew further by 15M) — solvent appears to **flip which pair
of stations is thermodynamically favored**, not just shift numbers around.
That's the headline result worth verifying before trusting it, and it's
gotten more pronounced with more sampling, not less.

Caveat: their `d` is a rigid-body wheel-translation scan coordinate; ours is a
dynamic wheel-O-centroid-projected-on-N···N-axis PLUMED CV. The topology match
(two central-ish wells, two stopper-ish wells, roughly symmetric) is solid;
the exact Å values are not guaranteed to be on identical origins/scales, so
treat the qualitative reordering as the finding, not the precise numbers.

## What this checks

Is the well-depth reordering a **real solvent effect**, or a **GAFF2
artifact**? GAFF2 (generic small-molecule force field, what drives the
explicit-solvent MetaD's non-bonded terms) may simply be getting the
wheel-rod non-covalent interactions wrong at one or both station types
(CH-π, amide H-bonding, CF₃ contacts) — reordering the wells for the wrong
reason. A higher-level QM single point on the *same* MD-sampled geometries is
a fast, direct way to test this.

## Plan (in order of effort)

### 1. Single-point QM rescoring (this folder, ready to run)

Pull representative frames from `../metad/metad_run/pymol_rotaxane.dcd` at
each well/saddle `d` value, single-point them with GFN2-xTB (tblite, the same
engine `vib_stations.py` uses), and compare relative energies to what the
MetaD/GAFF2 surface says at the same stations.

- `extract_stations.py` — recomputes `d` directly from the (already
  correctly-wrapped, solute-only) PyMOL trajectory using the same CV
  definition PLUMED used, and writes the nearest frame to each target `d` as
  a plain XYZ. Deliberately does NOT try to back-map to `COLVAR` rows, since
  `COLVAR`'s time axis resets across the resume segment (see
  `../metad/README.md`'s "time is in fs" gotcha) while the joined PyMOL
  trajectory's frame order is already correct end-to-end.
- `score_tblite.py` — GFN2-xTB single point (via `tblite`'s ASE calculator,
  charge=0, multiplicity=1 — closed-shell, matching the Rotaxanes repo's own
  convention) on each extracted frame, tabulates relative energies.

This is a **gas-phase single point of a solution-phase-sampled geometry** —
no implicit or explicit solvent in this script. It answers "does GAFF2 and
GFN2-xTB agree on the relative electronic energetics of these exact
geometries," not "what is the true solvent-corrected barrier." Good first
pass because it's cheap (minutes) and directly testable.

**Expected outcomes and what they'd mean:**
- If GFN2-xTB *also* ranks the stopper-like stations below the central-like
  ones (even without solvent) → the reordering isn't purely a solvation
  effect; GAFF2 was pointing at something real that the gas-phase pipeline's
  entropy correction (see below) may have been masking or that the vibrational
  correction shifted.
- If GFN2-xTB ranks them the *other* way (central deeper, matching the
  gas-phase pipeline) → the reordering is solvent-driven, and GAFF2's
  qualitative picture (deep stopper wells) needs the solvent context to trust;
  worth then checking with implicit solvent (below) before fully believing it.

### 2. Implicit-solvent xtb (ALPB) — cheap middle rung

`tblite`'s ASE wrapper has native implicit solvation (verified against the
installed 0.7.0 source, `tblite/ase.py`): pass `solvation=("alpb", "water")`
(a tuple, not separate kwargs) to `TBLite(...)` -- it forwards to
`Calculator.add("alpb-solvation", "water")` under the hood. Also supports
`"gbsa"` the same way. Re-running the extracted frames (or, better, a short
xtb-implicit-solvent re-scan analogous to `displace_wheel.py --engine
tblite`) with ALPB water would test whether the reordering shows up with
*implicit* polarity alone, or whether it specifically requires explicit
water/H-bonding (only visible in our TIP3P MD) to appear. Much cheaper than
the full explicit-solvent MetaD extension (minutes vs. ~20h), so a good
sanity check to run before committing to more MD.

### 3. Vibrational free-energy correction on our stations

Reuse `vib_stations.py`'s constrained partial-Hessian approach (rigid wheel +
rod-tip anchors fixed, so the reaction coordinate can't drift — the
established fix for the imaginary-mode problem a naive Hessian would hit)
applied to structures pulled from our MD wells, layering a QM-quality
ZPE/entropy correction on top of classical, solvent-inclusive configurational
sampling. Most effort (needs adapting the constraint set to our topology and
atom indices), best physics: MD supplies configurational + solvent entropy,
QM supplies accurate intramolecular vibrations/electronics.

## Package requirements

Not installed in the `openmm-md` conda env, and deliberately kept **out** of
it — this analysis is fully separate from the OpenMM/PLUMED MD (it only reads
already-finished trajectory output), and `openmm-md` has a delicately-pinned
solve (`openmm-plumed` pins `openmm` to 8.4; see `../metad/README.md`) not
worth risking by adding unrelated packages into that env's dependency graph.

Isolated in a **`uv`-managed `.venv`** in this folder (matches the Rotaxanes
repo's own convention: `.venv`, Python 3.12, `uv pip install`) rather than a
new conda env — separate, disposable, and doesn't touch the conda solver at
all:

```sh
cd rotaxanes/qm_rescoring
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python ase "tblite[ase]" mdtraj
```

Done — `.venv/` now has `ase==3.29.0`, `tblite==0.7.0` (with its
`tblite.ase.TBLite` ASE calculator), and `mdtraj==1.11.1` (so both scripts run
fully standalone from this folder, no conda env needed at all). Verified all
three import cleanly.

**This one env covers all three plan steps, not just step 1** — verified by
reading the installed packages' source, not assumed: step 2's ALPB solvation
is a `tblite` built-in (`tblite/ase.py`, forwards to
`Calculator.add("alpb-solvation", ...)`, no extra package), and step 3's
constrained partial-Hessian needs only `ase.vibrations.Vibrations` +
`ase.constraints.FixAtoms` + `ase.thermochemistry.HarmonicThermo`, all
standard `ase` modules already installed, driving the same `TBLite`
calculator. No further installs needed for any of the three.

**Not needed**: `fairchem-core` / `torch` / an `HF_TOKEN` — those are only for
the UMA MLIP engine (step 2/3 could optionally add a UMA cross-check later,
but it's heavier — brings torch, needs a HuggingFace token, and per the
Rotaxanes repo's own `CLAUDE.md`, **cannot share a process with tblite** —
both bundle their own `libomp` and segfault if loaded together, so a UMA path
would need its own subprocess/venv, same as their `bench_uma_tblite.py
--driver` pattern). If ever added, it'd get its own separate `.venv` for the
same reason tblite gets one here.

## Status (as of the rot1 20M-step, converged MetaD result)

Env is ready and everything below has actually been run at least once. Summary
of what worked, what didn't, and what's queued next:

1. **`score_tblite.py`** (raw single-point on 6 hand-picked stations) — ran,
   but the numbers were badly inflated (d=-4.90 came out +43 kcal/mol above
   the minimum, vs MD's 7.2) because a raw MD snapshot isn't a relaxed
   structure. **Superseded by (2).**
2. **`relax_and_score.py`** (loose LBFGS relax — one rod atom + one wheel
   `pinu`-style unit pinned, NOT both rod tips — then single point) — fixed
   the magnitude problem (same station: 43 -> 14 kcal/mol) while barely
   moving `d` (max 0.09 Å drift), and confirmed the terminal-wells-favored
   finding independently, both gas-phase and with `--solvent water` (ALPB).
   ALPB pushed the central region further above the terminal wells (same
   direction as the MD result), on the same solvent-derived geometries in
   both cases — see the file's docstring for the full caveat about that not
   being a clean gas-vs-solution test.
3. **Dense auto-grid scan** (`extract_stations.py --auto-grid`, 51 points at
   0.5 Å across the physically interesting -12..+12 Å span, `relax_and_score.py
   --solvent water`) — ran (~39 min), and **revealed the single-snapshot
   approach doesn't survive densification**: >10 kcal/mol swings between
   *adjacent* 0.5 Å grid points (visible in `qm_vs_md.png`, sent to the user).
   Each point was still one MD snapshot's post-relax energy, so
   snapshot-to-snapshot conformational noise (ring pucker, rotamer state)
   dominates over the underlying trend at that resolution.
4. **Frame averaging + dense scan — RUN, against the converged 20M-step MD
   result.** `extract_stations.py --auto-grid --n-frames 4` (default
   `--min-gap 500`) pulled 205 frames: 52 targets, most with 4 decorrelated
   MD replicates each (one sparse-region target only got 1 -- too few
   candidates within the window there). Launched detached via
   `run_dense_scan.sh` + `systemd-inhibit` (same pattern as the MD runs; full
   scope, ~2.6h estimated) and finished in ~1h55m (52 min faster than
   estimated). `relax_and_score.py` grouped replicates by target `d` and
   reported mean + std; `plot_comparison.py` plotted the mean with error bars
   against `Fes.dat` (`qm_vs_md.png`, sent to the user).

   **Result: averaging fixed the point-to-point noise.** Std devs across
   replicates are ~1-5 kcal/mol (real, visible residual noise) instead of the
   >40 kcal/mol adjacent-point swings the single-frame dense scan showed. The
   qualitative picture holds: low points cluster near both terminal regions
   (e.g. d=-11.00 at 0.0, d=+12.00 at 1.83, d=+9.50 at 4.26 kcal/mol) with the
   middle broadly higher (~8-15 kcal/mol) and bumpy -- consistent with the
   MD's own bumpy central plateau (rot1_shuttles.md's "central plateau of
   soft-mode rattling" description), not a clean two-well curve.

   **Gotcha hit:** `run_dense_scan.sh` didn't set `PYTHONUNBUFFERED=1`
   (unlike the MD launch scripts, which do), so Python block-buffered stdout
   when redirected to a file -- the log showed 0 lines for the first ~30 min
   despite the process genuinely computing (279% CPU, confirmed via `ps`).
   Worked around by estimating progress from wall-clock/CPU time instead of
   log line count for that run; fixed the script for next time (added the
   env var, matching the MD scripts' convention) -- doesn't retroactively fix
   the run that already happened, but the next one will log live.

5. **Eyring/TST rates** (`../metad/eyring.py`, new) — applied
   `k = (kB*T/h)*exp(-DeltaG/RT)` to the MD barriers at each stage and to the
   gas-phase pipeline's reference values, validated against
   `rot1_shuttles.md`'s own published numbers (matches within rounding). The
   converged 20M-step barrier (12.59 kcal/mol) implies tau ~ 237.5 us --
   between the gas-phase Shuttle A (central saddle, 1.36 us) and Shuttle B
   (stopper saddle, 3.85 ms) timescales, and ~24x slower than rot2's
   converged tau (9.97 us).

6. **Vibrational free-energy correction** (original plan's step 3) — not
   started. Now that a frame-averaged, lower-noise electronic curve exists,
   the natural next move is `ase.vibrations.Vibrations` + `FixAtoms` (same
   one rod atom + one wheel unit already used for the loose relax --
   consistent with what `vib_stations.py`'s "constrained partial Hessian
   avoids the reaction-coordinate imaginary mode" logic needs) +
   `HarmonicThermo`, to add ZPE/entropy on top of the electronic energy at
   whichever stations end up mattering most (the two terminal wells + the
   barrier region).

### Commands (already run once against the 20M-step result; re-run if the MD changes again)

```sh
cd rotaxanes/qm_rescoring
rm -f frames/*.xyz                            # clear any prior extraction
.venv/bin/python extract_stations.py --auto-grid --n-frames 4   # or a coarser scope
./run_dense_scan.sh                            # detached + systemd-inhibit, ~2h
.venv/bin/python plot_comparison.py            # -> qm_vs_md.png, mean + error bars
python ../metad/eyring.py <barrier_kcal>       # rate/timescale for any barrier
```
