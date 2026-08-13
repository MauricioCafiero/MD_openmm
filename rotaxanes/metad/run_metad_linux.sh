#!/bin/bash
# Detached + systemd-inhibit-held WT-METAD run of the rot1 rotaxane shuttle.
# Linux port of run_metad.sh (which is macOS/zsh/caffeinate-specific). This
# machine has no GPU -- only Reference/CPU OpenMM platforms exist -- so we run
# on CPU explicitly (also the CLAUDE.md-preferred platform for metad stability
# regardless of GPU availability).
#
# Multi-hour job -> launched via the Popen(start_new_session=True) one-liner
# (see rotaxanes/CLAUDE.md), NOT via the Bash tool's run_in_background.
source ~/miniforge3/etc/profile.d/conda.sh
conda activate openmm-md
export PYTHONUNBUFFERED=1        # flush prints so the log is live
export OPENMM_CPU_THREADS=8      # this machine's core count
cd /home/cafierom/python_linux/MD_openmm/rotaxanes/metad
set -e

echo "=== METAD start $(date) ==="
systemd-inhibit --what=sleep:idle --who="rot1-metad" \
    --why="rot1 WT-METAD run (multi-hour, must not suspend)" --mode=block \
    python run_metad.py \
        --system ../../outputs/system.xml --topology ../../outputs/complex.pdb \
        --plumed plumed.dat --out-dir metad_run \
        --steps "${1:-25000000}" --platform CPU
echo "=== DONE $(date) ==="
echo
echo "Reconstruct the free-energy profile F(d) with:"
echo "  plumed sum_hills --hills metad_run/HILLS --mintozero --outfile Fes.dat"
