#!/bin/bash
# Independent replay of the local theorem: 0 < |p|_inf <= 0.04  =>  xi(Q_p) < 5/2.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$HOME/.venvs/analysis-py313/bin/python}
mkdir -p replay
$PY lemma_firstorder.py > replay/lemma_firstorder.log
$PY lemma_residual.py   > replay/lemma_residual.log
grep -q "VERIFIED" replay/lemma_firstorder.log && grep -q "VERIFIED" replay/lemma_residual.log
$PY assemble.py 8 10 12 8/100 replay/coeffs.pkl > replay/assemble.log
$PY schur_mtm.py replay/coeffs.pkl replay/schur.pkl 36,38,40,42,43,44,45,46,46.5,47,47.5,48,48.5,49 > replay/schur.log
$PY cover2.py replay/schur.pkl 4 16 8 > replay/cover.log
$PY - <<'PYEOF'
import json, re
line = [l for l in open("replay/cover.log") if l.startswith("SUMMARY")][0]
S = json.loads(line[len("SUMMARY "):])
assert S["failures"] == 0, S
from fractions import Fraction as Fr; assert Fr(S["min_smax"]) >= Fr(1,25), S
print(f"LOCAL THEOREM REPLAYED: {S['certified_boxes']} direction boxes, 0 failures, cube half-width {S['min_smax']} (exact), {S['seconds']} s")
PYEOF
