# Notes for an automated rotaxane WT-METAD run

This folder is a self-contained rotaxane build + WT-MetaD + analysis pipeline.
The full process, commands, results, and gotchas are in `metad/README.md` — read
that first. This file is the short operational guardrail list for running it
end-to-end (the things that are easy to get wrong and not obvious from the code).

## The task
Run a rotaxane WT-METAD shuttle simulation on a rot — `rot1` (inputs already in
this folder) or a new rot — and regenerate the FES + PyMOL deliverables.

## To run a new rot
1. Write `rotaxanes/<name>.txt` with two lines: `rod: <SMILES>` and `wheel: <SMILES>`.
2. Provide `rotaxanes/<name>_displaced_pinu.xyz`: threaded, scan-relaxed starting
   coords (coords only, no bonds). The atom count MUST equal the SMILES atom
   count (rod heavy+H + wheel heavy+H) — `build_rotaxane.py` asserts this.
3. The rod SMILES must contain exactly 2 amide N's and the wheel SMILES exactly
   8 crown-ether O's — `make_plumed.py` asserts both (the generic CV is residue
   ROD N's + residue WHL O's, by element).
4. Follow `metad/README.md` step 0 → 1 → 2.

## Operational guardrails (the part that bites)
- **Multi-hour job — launch DETACHED, not as a tracked background task.** Use
  `caffeinate -w $$ &` near the top of the run script, then launch via
  `subprocess.Popen(..., start_new_session=True)` from a normal (non-background)
  shell call. Do NOT use the Bash tool's `run_in_background: true` — the harness
  reaps tracked tasks ~30–50 min in. See the user's global CLAUDE.md for the
  exact pattern. Verify with `ps -ax -o pid,ppid,sess,command` (expect PPID 1).
- **Prefer CPU over OpenCL for metad.** On Apple Silicon, OpenCL is deprecated
  and silently crashes on multi-hour runs (no python traceback, no sleep event,
  not a harness reap — the rot2htpuma run died twice this way mid-run on OpenCL).
  CPU is ~2.9 ms/step vs ~0.9 OpenCL but stable. If you do use OpenCL, the
  `CheckpointReporter` in `run_metad_resume.py` (every 50k steps) makes a death
  recoverable by `Simulation.loadCheckpoint`.
- **Resume after a death:** `run_metad_resume.py` (+ `run_resume.sh`) does a warm
  restart from the last `traj.dcd` frame and puts PLUMED in RESTART mode so
  HILLS/COLVAR are appended (not truncated). New reporters (`traj_resume.dcd`,
  `energy_resume.csv`) keep the originals intact.
- **Do NOT add LOWER/UPPER_WALLS.** The wheel never deshuttles: the CV is plain
  `DISTANCE ... COMPONENTS` with PBC on (minimum-image), so it is bounded by
  L/2; COLVAR only ever reaches ±~7.9 Å. The "wheel escapes to 48 Å" story was a
  mis-unwrapping artifact, not physics.
- **Deliverables come from `analyze_metad.py`**, in `metad_run/`:
  `Fes.dat` + `Fes_plot.png` (the free-energy surface F(d), kcal/mol) and
  `pymol_rotaxane.pdb` + `pymol_rotaxane.dcd` (solute-only, rod-centred,
  make-whole + min-image wrap — all frames kept, no PBC snap). The
  `energy_plot.png` is thermodynamic (E/T/V vs time), NOT the FES.

## Env
`conda env create -f environment.yml` (at repo root) is metad-ready as of commit
789dd38 — it includes `openmm-plumed` (which brings the `plumed` CLI and pins
openmm to 8.4). Then `pip install -e .` for the `omd` console script.