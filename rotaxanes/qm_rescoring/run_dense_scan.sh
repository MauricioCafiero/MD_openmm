#!/bin/bash
# Detached + systemd-inhibit-held frame-averaged QM rescoring run.
# ~2.6h for the full 205-frame scope -- same auto-suspend risk as the MD
# runs (CPU-busy != "not idle" for this machine's 1h timeout), so launched
# the same way: systemd-inhibit + Popen(start_new_session=True), NOT
# run_in_background.
cd /home/cafierom/python_linux/MD_openmm/rotaxanes/qm_rescoring
export PYTHONUNBUFFERED=1        # flush prints so the log is live, not buffered
set -e

echo "=== DENSE SCAN start $(date) ==="
systemd-inhibit --what=sleep:idle --who="rot1-qm-dense-scan" \
    --why="rot1 QM frame-averaged rescoring (multi-hour, must not suspend)" --mode=block \
    .venv/bin/python relax_and_score.py --solvent water
echo "=== DONE $(date) ==="
