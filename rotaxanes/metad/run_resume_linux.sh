#!/bin/bash
# Detached + systemd-inhibit-held RESUME of the rot1 WT-METAD run: warm
# restart from the last traj.dcd frame + PLUMED RESTART (appends to HILLS).
# Linux port of run_resume.sh (macOS/zsh/caffeinate-specific). See
# run_metad_linux.sh for the platform/inhibitor rationale.
#
# Multi-hour job -> launched via the Popen(start_new_session=True) one-liner
# (see rotaxanes/CLAUDE.md), NOT via the Bash tool's run_in_background.
source ~/miniforge3/etc/profile.d/conda.sh
conda activate openmm-md
export PYTHONUNBUFFERED=1
export OPENMM_CPU_THREADS=8
cd /home/cafierom/python_linux/MD_openmm/rotaxanes/metad
set -e

echo "=== METAD RESUME start $(date) ==="
systemd-inhibit --what=sleep:idle --who="rot1-metad-resume" \
    --why="rot1 WT-METAD resume (multi-hour, must not suspend)" --mode=block \
    python run_metad_resume.py \
        --system ../../outputs/system.xml --topology ../../outputs/complex.pdb \
        --plumed plumed.dat --out-dir metad_run \
        --add-steps "${1:-5000000}" --platform CPU
echo "=== DONE $(date) ==="
echo
echo "Reconstruct F(d) with:"
echo "  plumed sum_hills --hills metad_run/HILLS --mintozero --outfile Fes.dat"
