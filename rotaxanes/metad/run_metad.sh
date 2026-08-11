#!/bin/zsh
# Detached + caffeinate-held WT-METAD run of the rotaxane shuttle.
# Multi-hour job -> launch via the Popen(start_new_session=True) one-liner
# (see the launch pattern below), NOT via run_in_background.
#
# Pre-req: conda install -c conda-forge openmm-plumed  (into the openmm-md env)
source ~/miniforge3/etc/profile.d/conda.sh
conda activate openmm-md
export DYLD_FALLBACK_LIBRARY_PATH="/opt/homebrew/opt/libomp/lib:/lib:/usr/lib${DYLD_FALLBACK_LIBRARY_PATH:+:$DYLD_FALLBACK_LIBRARY_PATH}"
export OMP_NUM_THREADS=4
export PYTHONUNBUFFERED=1            # flush prints so the log is live
cd /Users/cafierom/python_mac/openmm/rotaxanes/metad
caffeinate -w $$ &
set -e

# (re)generate plumed.dat from the current topology (cheap; idempotent)
python make_plumed.py --topology ../outputs/complex.pdb --out plumed.dat

echo "=== METAD start $(date) ==="
python run_metad.py \
    --system ../outputs/system.xml --topology ../outputs/complex.pdb \
    --plumed plumed.dat --out-dir metad_run \
    --steps "${1:-5000000}" --platform auto
echo "=== DONE $(date) ==="
echo
echo "Reconstruct the free-energy profile F(d) with:"
echo "  plumed sum_hills --hills metad_run/HILLS --mintozero --outfile Fes.dat"