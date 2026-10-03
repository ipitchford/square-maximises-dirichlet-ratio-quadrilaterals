#!/bin/bash
# run a sequence of production jobs in one process slot; arguments are "FAM SHARD NSHARDS" triples
cd /Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts
while [ $# -ge 3 ]; do ./prod_fam.sh $1 $2 $3; shift 3; done
