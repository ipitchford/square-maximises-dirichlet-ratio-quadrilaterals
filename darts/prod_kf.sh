#!/bin/bash
# Production: spike/collapse face, C in KF (high-xi zone), FEM spike charts at k >= 5.  Resumable; one process.
cd /Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
PY=/Users/admin/.venvs/analysis-py313/bin/python
nice -n 10 $PY -u run_spike_all.py KF COLLAR certs2/kf_collar.jsonl 0 1 >> certs2/kf_collar.out 2>&1
nice -n 10 $PY -u run_spike_all.py KF SPIKE certs2/kf_spike.jsonl 0 1 >> certs2/kf_spike.out 2>&1
