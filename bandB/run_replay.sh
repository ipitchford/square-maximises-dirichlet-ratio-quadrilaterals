#!/bin/bash
# run_replay.sh — fail-closed replay (replay.jl) of the re-certified logs, nice 10.  Author: Anonymous.
# exit status = replay.jl's: 0 PASSED, 1 FAILED, 3 test mode.
set -euo pipefail
cd "$(dirname "$0")"
T=${T:-3}
echo "replay: $(date) threads $T"
code=0
nice -n 10 /opt/homebrew/bin/julia -t "$T" replay.jl runs/recert_chart1.jsonl runs/recert_chart2.jsonl > replay.log 2>&1 || code=$?
echo "replay exit $code $(date)"
exit "$code"
