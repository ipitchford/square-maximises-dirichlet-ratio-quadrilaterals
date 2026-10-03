"""Re-check boxes NOT verified by xsample.py, with xcert.certify_adaptive (range subdivision). Usage: xsample2.py in.jsonl out.jsonl"""
import json, sys, time
from fractions import Fraction as Fr
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from xcert import certify_adaptive
src, out = sys.argv[1], sys.argv[2]
R = [json.loads(l) for l in open(src)]
todo = [r for r in R if not r["verified"]] if "--all" not in sys.argv else R
print(f"{len(todo)} boxes to check", flush=True); nv = 0; t0 = time.time()
with open(out, "w") as f:
    for i, r in enumerate(todo):
        t = time.time(); ok, st, beta = certify_adaptive([Fr(v) for v in r["lo"]], [Fr(v) for v in r["hi"]])
        nv += ok
        f.write(json.dumps({"path": r["path"], "verified": ok, "subboxes": st["leaves"], "beta": beta, "agent_kappa": r["agent_kappa"], "sec": round(time.time()-t, 1)}) + "\n"); f.flush()
        print(f"[{i+1}/{len(todo)}] {r['path'][:40]} ok={ok} subboxes={st['leaves']} {time.time()-t:.0f}s (verified so far {nv})", flush=True)
print(f"DONE verified {nv}/{len(todo)}, {time.time()-t0:.0f}s", flush=True)
