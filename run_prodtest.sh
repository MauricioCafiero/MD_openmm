#!/bin/zsh
# Detached + caffeinate-held production test (run -> analyze) for the SULT1A3
# dimer + PAP + ligand system (~1 ns / 500k steps). Mirrors the Rotaxanes launch
# pattern: caffeinate -w $$ dies with this script, holding off macOS idle sleep
# for the whole run; launched detached via a
# subprocess.Popen(..., start_new_session=True) one-liner so it survives the
# ~45-min background-task reap (macOS has no setsid binary; start_new_session
# calls setsid() in the child -> new session, reparented to launchd).
source ~/miniforge3/etc/profile.d/conda.sh
conda activate openmm-md
export DYLD_FALLBACK_LIBRARY_PATH="/opt/homebrew/opt/libomp/lib:/lib:/usr/lib${DYLD_FALLBACK_LIBRARY_PATH:+:$DYLD_FALLBACK_LIBRARY_PATH}"
export OMP_NUM_THREADS=4
cd /Users/cafierom/python_mac/openmm
caffeinate -w $$ &          # hold off idle sleep for the whole script
set -e
echo "=== RUN start $(date) ==="
omd run --system outputs/prodtest/system.xml --topology outputs/prodtest/complex.pdb \
        --steps 500000 --out-dir outputs/prodtest
echo "=== ANALYZE start $(date) ==="
omd analyze --traj outputs/prodtest/traj.dcd --topology outputs/prodtest/complex.pdb \
           --out-dir outputs/prodtest
echo "=== DONE $(date) ==="