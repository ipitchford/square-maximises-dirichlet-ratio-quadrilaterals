"""Independent re-certification of sampled transport boxes (verify/xcert.py), combined strategy:
(1) one certificate for the box with adaptive range subdivision; (2) if that fails, bisect the BOX (fresh test functions,
fresh beta) and recurse.  A box is VERIFIED iff every leaf of this recursion is certified.  Usage: xsample3.py in.jsonl out.jsonl"""
import json, sys, time
from fractions import Fraction as Fr
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from xcert import certify_adaptive
class Budget(Exception): pass
import xcover as _XC     # independent exact exclusion predicates (E2 Lemma A, E4 non-convex, N no convex point)
def excluded(lo, hi):
    lo_f = [float(v) for v in lo]; hi_f = [float(v) for v in hi]
    if any(Fr(a) != b for a, b in zip(lo_f, lo)) or any(Fr(a) != b for a, b in zip(hi_f, hi)):
        lo_f, hi_f = lo, hi                     # predicates accept exact Fractions as well
    return _XC.pred("N", lo_f, hi_f) or _XC.pred("E4", lo_f, hi_f) or _XC.pred("E2", lo_f, hi_f)
CERT_BUDGET = 40          # work units (certificate setups + range sub-boxes) per sampled box
def combined(lo, hi, depth, maxdepth=7, stats=None):
    if excluded(lo, hi): stats["excluded"] = stats.get("excluded", 0) + 1; return True
    if stats["certs"] + stats["subboxes"] >= CERT_BUDGET: raise Budget()      # budget in measured work units (~21 s each under load)
    ok, st, beta = certify_adaptive(lo, hi, maxdepth=8, maxleaves=32)
    stats["certs"] += 1; stats["subboxes"] += st["leaves"]
    if ok: return True
    if depth >= maxdepth: return False
    w = [h-l for l, h in zip(lo, hi)]; d = w.index(max(w)); mid = (lo[d]+hi[d])/2
    h1 = list(hi); h1[d] = mid; l2 = list(lo); l2[d] = mid
    return combined(lo, h1, depth+1, maxdepth, stats) and combined(l2, hi, depth+1, maxdepth, stats)
if __name__ == "__main__":
    src, out = sys.argv[1], sys.argv[2]
    import hashlib
    _h = lambda fn: hashlib.sha256(open(fn, "rb").read()).hexdigest()
    _here = __file__.rsplit("/", 1)[0]
    print("PROVENANCE", json.dumps({"input": _h(src), "xsample3.py": _h(__file__), "xcert.py": _h(_here + "/xcert.py")}), flush=True)
    R = [json.loads(l) for l in open(src)]
    todo = [r for r in R if not r["verified"]] if "--all" not in sys.argv else R
    print(f"{len(todo)} boxes to check", flush=True); nv = 0; t0 = time.time()
    with open(out, "w") as f:
        for i, r in enumerate(todo):
            t = time.time(); stats = {"certs": 0, "subboxes": 0}
            try:
                ok = combined([Fr(v) for v in r["lo"]], [Fr(v) for v in r["hi"]], 0, stats=stats); stats["outcome"] = "verified" if ok else "failed"
            except Budget:
                ok = False; stats["outcome"] = "budget"
            nv += ok
            f.write(json.dumps({"path": r["path"], "verified": ok, **stats, "agent_kappa": r["agent_kappa"], "sec": round(time.time()-t, 1)}) + "\n"); f.flush()
            print(f"[{i+1}/{len(todo)}] {r['path'][:40]} {stats['outcome']} certs={stats['certs']} subboxes={stats['subboxes']} {time.time()-t:.0f}s (verified so far {nv})", flush=True)
    print(f"DONE verified {nv}/{len(todo)}, {time.time()-t0:.0f}s", flush=True)
