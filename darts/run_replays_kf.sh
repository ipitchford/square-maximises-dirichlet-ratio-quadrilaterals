#!/bin/bash
# fail-closed replays of the KF runs (3 shards for the spike log, 1 for the collar log); exit codes recorded per shard
cd /Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
PY=/Users/admin/.venvs/analysis-py313/bin/python
for i in 0 1 2; do
  ( nice -n 10 $PY -u replay2.py certs2/kf_spike.jsonl KF_SPIKE manifests/kf_spike.json --shard $i/3 --out replays/kf_spike.$i.jsonl > replays/kf_spike.$i.log 2>&1; echo "exit $?" >> replays/kf_spike.$i.log ) &
done
( nice -n 10 $PY -u replay2.py certs2/kf_collar.jsonl KF_COLLAR manifests/kf_collar.json --out replays/kf_collar.jsonl > replays/kf_collar.log 2>&1; echo "exit $?" >> replays/kf_collar.log ) &
wait
