#!/bin/zsh
# Detached + caffeinate-held RESUME of the WT-METAD run (died at step ~11.944M).
# Warm restart from the last traj.dcd frame + PLUMED RESTART (appends to HILLS).
# Multi-hour job -> launch via the Popen(start_new_session=True) one-liner below,
# NOT via run_in_background.
source ~/miniforge3/etc/profile.d/conda.sh
conda activate openmm-md
export DYLD_FALLBACK_LIBRARY_PATH="/opt/homebrew/opt/libomp/lib:/lib:/usr/lib${DYLD_FALLBACK_LIBRARY_PATH:+:$DYLD_FALLBACK_LIBRARY_PATH}"
export OMP_NUM_THREADS=4
export PYTHONUNBUFFERED=1
cd /Users/cafierom/python_mac/openmm/rotaxanes/metad
caffeinate -w $$ &
set -e

echo "=== METAD RESUME start $(date) ==="
python run_metad_resume.py \
    --system ../outputs/system.xml --topology ../outputs/complex.pdb \
    --plumed plumed.dat --out-dir metad_run --platform auto
echo "=== DONE $(date) ==="
echo
echo "Reconstruct F(d) with:"
echo "  plumed sum_hills --hills metad_run/HILLS --mintozero --outfile Fes.dat"