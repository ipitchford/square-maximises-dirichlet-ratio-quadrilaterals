#!/bin/bash
# number of LOUD FAILURE lines written by the pool workers since their last restart marker
cd "$(dirname "$0")"
n=0
for f in pool_runs/w*.err; do
  c=$(awk '/^=== restart/{c=0} /LOUD FAILURE/{c++} END{print c+0}' "$f"); n=$((n + c))
done
echo $n
