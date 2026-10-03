#!/bin/bash
# run several families sequentially in one process slot: prod_seq.sh FAM1 FAM2 ...
cd /Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts
for F in "$@"; do ./prod_fam.sh $F 0 1; done
