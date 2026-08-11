#!/bin/zsh
# First solvated rotaxane MD run: rot2htpuma (rod + 24-crown-8 wheel), 1 ns / 500k
# steps, 300 K, NPT, on whichever platform auto-selects (OpenCL single on this Mac).
# Detached + caffeinate-held so it survives the ~30-50 min background-task reap.
#
# Launch with the Popen(start_new_session=True) one-liner in run_rot2htpuma.sh --
# do NOT run this via run_in_background; it must be untracked (PPID 1, own session).
source ~/miniforge3/etc/profile.d/conda.sh
conda activate openmm-md
export DYLD_FALLBACK_LIBRARY_PATH="/opt/homebrew/opt/libomp/lib:/lib:/usr/lib${DYLD_FALLBACK_LIBRARY_PATH:+:$DYLD_FALLBACK_LIBRARY_PATH}"
export OMP_NUM_THREADS=4
cd /Users/cafierom/python_mac/openmm/rotaxanes
caffeinate -w $$ &          # dies with this script; holds off macOS idle sleep
set -e

# --- stage 1: assemble the 2-fragment SDF (bonds from SMILES, coords from the
#     displaced XYZ already in the repo root) ---
if [[ ! -f rot2htpuma.sdf ]]; then
  echo "=== build SDF $(date) ==="
  python build_rotaxane.py --smiles rot2htpuma.txt \
        --from-xyz ../rot2htpuma_displaced_pinu.xyz --out rot2htpuma.sdf
fi

# --- stage 2: solvate (skip if already built; padding 1.2 nm balances PBC
#     wrapping safety vs cost -- 3144 particles, ~1008 waters) ---
if [[ ! -f outputs/system.xml ]]; then
  echo "=== build-multimol $(date) ==="
  omd build-multimol --sdf rot2htpuma.sdf --out-dir outputs \
                    --padding 1.2 --resnames ROD WHL
fi

# --- stage 3: 1 ns production at 300 K, NPT, auto platform ---
echo "=== RUN start $(date) ==="
omd run --system outputs/system.xml --topology outputs/complex.pdb \
        --steps 500000 --platform auto --out-dir outputs

# --- stage 4: analyze (RMSD/RMSF + energy plots) ---
echo "=== ANALYZE start $(date) ==="
omd analyze --traj outputs/traj.dcd --topology outputs/complex.pdb --out-dir outputs
echo "=== DONE $(date) ==="