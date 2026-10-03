#!/bin/bash
# Launch pool workers FIRST..LAST (pool.py) at nice 10, one CPU each; output appended to pool_runs/w<ID>.out / .err.
# All workers get the same family list (same pseudo-random root order).  Usage: ./pool_launch.sh FIRST LAST
cd "$(dirname "$0")"
PY=/Users/admin/.venvs/analysis-py313/bin/python
mkdir -p pool_runs claims
for ((i=$1; i<=$2; i++)); do
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 nohup nice -n 10 $PY -u pool.py $i PP SP_spk SP_col PF FP \
    >> pool_runs/w$i.out 2>> pool_runs/w$i.err &
  echo "worker $i pid $!"
done
