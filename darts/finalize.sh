#!/bin/bash
# After the pool has finished (every worker printed DONE): stop the scaler, merge every family's logs and run the full fail-closed
# replay (replay_all.sh) with SHARDS processes; then print the runs table (summarize.py).  Usage: ./finalize.sh [SHARDS]
cd "$(dirname "$0")"
PY=/Users/admin/.venvs/analysis-py313/bin/python
SH=${1:-7}
pkill -f pool_scaler.sh; pkill -f pool_keepalive.sh
if pgrep -f "pool.py [0-9]" > /dev/null; then echo "pool workers still running"; exit 1; fi
for F in PP SP_col SP_spk PF FP FF SF_col SF_spk; do $PY merge_pool.py $F || { echo "merge failed: $F"; exit 1; }; done
SHARDS=$SH ./replay_all.sh > replays/replay_all.out 2>&1; st=$?
tail -40 replays/replay_all.out
$PY summarize.py > replays/runs_table.md
cat replays/runs_table.md
exit $st
