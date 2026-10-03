#!/bin/bash
# Replays of the first-session and KF logs (replay2.py, full re-certification), one process, in the exact form used by
# replay_all.sh (replays/<tag>.0.jsonl, shard 0/1), so that replay_all.sh can reuse them while log and verifier are unchanged.
cd "$(dirname "$0")"
PY=/Users/admin/.venvs/analysis-py313/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
while read log fam man tag; do
  rm -f replays/$tag.[0-9]*.jsonl replays/$tag.[0-9]*.log
  nice -n 10 $PY replay2.py $log $fam $man --shard 0/1 --out replays/$tag.0.jsonl > replays/$tag.0.log 2>&1
  echo "$tag exit $? $(grep -h '"leaves_checked"' replays/$tag.0.log | cut -c1-200)"
done << 'LIST'
certs/crit_collar.jsonl COLLAR manifests/crit_collar.json crit_collar
certs/crit_spike_X1.jsonl SPIKE_X1 manifests/crit_spike_X1.json crit_spike_X1
certs/main_cbox.jsonl MAIN_CBOX manifests/main_cbox.json main_cbox
certs2/kf_collar.jsonl KF_COLLAR manifests/kf_collar.json kf_collar
certs/shallow_crit.jsonl BOX manifests/shallow_crit.json shallow_crit
certs2/kf_spike.jsonl KF_SPIKE manifests/kf_spike.json kf_spike
LIST
echo LEGACY_DONE
