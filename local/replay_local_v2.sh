#!/bin/bash
# Independent replay of the local theorem: 0 < |p|_inf <= 1/20  =>  xi(Q_p) < 5/2.
# Fail-closed: the code is checked against MANIFEST.sha256 before anything runs, every stage asserts its own conclusions, and the
# run ends with a receipt binding the SHA-256 of the code, every output and the summary (replay_v2/RECEIPT.json).
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$HOME/.venvs/analysis-py313/bin/python}
CODE="engine2.py assemble.py localcert.py mtm.py schur_mtm.py cover2.py lemma_firstorder.py lemma_residual.py dump_local_data.py replay_local_v2.sh"
for f in $CODE; do grep -q "  $f\$" MANIFEST.sha256 || { echo "MANIFEST lacks $f"; exit 1; }; done
grep -E "  ($(echo $CODE | tr ' ' '|'))\$" MANIFEST.sha256 | shasum -a 256 -c - > /dev/null || { echo "code differs from MANIFEST.sha256"; exit 1; }
mkdir -p replay_v2
$PY lemma_firstorder.py > replay_v2/lemma_firstorder.log
$PY lemma_residual.py   > replay_v2/lemma_residual.log
grep -q "VERIFIED" replay_v2/lemma_firstorder.log && grep -q "VERIFIED" replay_v2/lemma_residual.log
$PY assemble.py 10 12 14 1/10 replay_v2/coeffs.pkl > replay_v2/assemble.log
$PY schur_mtm.py replay_v2/coeffs.pkl replay_v2/schur.pkl 36,37,38,39,40,41,42,43,44,45,46,46.5,47,47.5,48,48.5,49 > replay_v2/schur.log
$PY cover2.py replay_v2/schur.pkl 4 16 10 > replay_v2/cover.log
$PY dump_local_data.py replay_v2/coeffs.pkl replay_v2/schur.pkl replay_v2 > replay_v2/dump.log
$PY - <<'PYEOF'
import json, hashlib, os, sys, platform
line = [l for l in open("replay_v2/cover.log") if l.startswith("SUMMARY")][0]
S = json.loads(line[len("SUMMARY "):])
assert S["failures"] == 0, S
from fractions import Fraction as Fr; assert Fr(S["min_smax"]) >= Fr(1,20), S
h = lambda f: hashlib.sha256(open(f, "rb").read()).hexdigest()
code = "engine2.py assemble.py localcert.py mtm.py schur_mtm.py cover2.py lemma_firstorder.py lemma_residual.py dump_local_data.py replay_local_v2.sh".split()
outs = sorted(f for f in os.listdir("replay_v2") if f != "RECEIPT.json" and os.path.isfile(os.path.join("replay_v2", f)))
import flint
rec = dict(claim="0 < |p|_inf <= 1/20 implies xi(Q_p) < 5/2 (Theorem 2)", summary=S,
           code={f: h(f) for f in code}, outputs={f: h(os.path.join("replay_v2", f)) for f in outs},
           python=sys.version, python_flint=getattr(flint, "__version__", None), platform=platform.platform())
json.dump(rec, open("replay_v2/RECEIPT.json", "w"), indent=1, sort_keys=True)
print(f"LOCAL THEOREM REPLAYED: {S['certified_boxes']} direction boxes, 0 failures, cube half-width {S['min_smax']} (exact), {S['seconds']} s")
print("RECEIPT replay_v2/RECEIPT.json sha256", h("replay_v2/RECEIPT.json"))
PYEOF
