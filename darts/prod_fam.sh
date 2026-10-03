#!/bin/bash
# production: one cover5 family shard (resumable); usage: prod_fam.sh FAMILY SHARD NSHARDS
cd /Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
F=$1; I=$2; N=$3; L=$(echo $F | tr A-Z a-z)
nice -n 10 /Users/admin/.venvs/analysis-py313/bin/python -u run_fam.py $F certs3/${L}.$I.jsonl $I $N >> certs3/${L}.$I.out 2>> certs3/${L}.$I.err
echo "EXIT $?" >> certs3/${L}.$I.out
