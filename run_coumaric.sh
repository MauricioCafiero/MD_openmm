#!/bin/zsh
# Single-molecule production run of m-coumaric acid with a 300 -> 400 K linear
# temperature ramp over 500k steps (1 ns), NVT. Detached (caffeinate -w $$ +
# start_new_session). ~2.4h on CPU for this 2604-particle system.
source ~/miniforge3/etc/profile.d/conda.sh
conda activate openmm-md
cd /Users/cafierom/python_mac/openmm
export OMP_NUM_THREADS=4
caffeinate -w $$ &          # dies with this script; holds off macOS idle sleep
echo "=== RUN start $(date) ==="
omd run --system outputs/coumaric/system.xml --topology outputs/coumaric/complex.pdb \
        --steps 500000 --ramp-end 400 --pressure 0 --platform CPU --out-dir outputs/coumaric
omd analyze --traj outputs/coumaric/traj.dcd --topology outputs/coumaric/complex.pdb --out-dir outputs/coumaric
echo "=== DONE $(date) ==="