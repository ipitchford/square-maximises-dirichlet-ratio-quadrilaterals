#!/bin/bash
# run_recert.sh — re-certify both chart covers with the current (fixed) cert.jl (recert.jl).  Author: Anonymous.
# New logs runs/recert_chart{1,2}.jsonl; the original logs runs/cover_chart{1,2}.jsonl are only read.
set -euo pipefail
cd "$(dirname "$0")"
T=${T:-3}
for chart in 2 1; do
  echo "recert chart $chart: $(date) threads $T"
  if ! nice -n 10 /opt/homebrew/bin/julia -t "$T" recert.jl "$chart" "runs/cover_chart$chart.jsonl" "runs/recert_chart$chart.jsonl" > "runs/recert_chart$chart.out" 2>&1; then
    echo "recert chart $chart FAILED $(date)"; exit 1
  fi
  echo "recert chart $chart exit 0 $(date)"
done
