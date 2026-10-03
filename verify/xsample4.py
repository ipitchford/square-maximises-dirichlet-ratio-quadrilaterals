"""xsample4: xsample3 with the per-box work budget from $XCERT_BUDGET and the boxes processed in parallel ($XCERT_WORKERS).
Independent re-certification of sampled transport boxes (verify/xcert.py), combined strategy:
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
CERT_BUDGET = int(__import__("os").environ.get("XCERT_BUDGET", "40"))   # work units per sampled box (default 40 as in run 6)
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
def one(r):
    t = time.time(); stats = {"certs": 0, "subboxes": 0}
    try:
        ok = combined([Fr(v) for v in r["lo"]], [Fr(v) for v in r["hi"]], 0, stats=stats); stats["outcome"] = "verified" if ok else "failed"
    except Budget:
        ok = False; stats["outcome"] = "budget"
    return {"path": r["path"], "verified": ok, **stats, "budget": CERT_BUDGET, "agent_kappa": r["agent_kappa"], "sec": round(time.time()-t, 1)}
if __name__ == "__main__":
    import hashlib, os, multiprocessing as mp
    src, out = sys.argv[1], sys.argv[2]
    _h = lambda fn: hashlib.sha256(open(fn, "rb").read()).hexdigest()
    _here = __file__.rsplit("/", 1)[0]
    print("PROVENANCE", json.dumps({"input": _h(src), "xsample4.py": _h(__file__), "xcert.py": _h(_here + "/xcert.py"), "xcover.py": _h(_here + "/xcover.py"), "budget": CERT_BUDGET}), flush=True)
    todo = [json.loads(l) for l in open(src)]
    W = int(os.environ.get("XCERT_WORKERS", "1")); t0 = time.time(); nv = 0
    print(f"{len(todo)} boxes, {W} workers, budget {CERT_BUDGET}", flush=True)
    with open(out, "w") as f, mp.get_context("spawn").Pool(W) as pool:
        for i, res in enumerate(pool.imap_unordered(one, todo)):
            nv += res["verified"]; f.write(json.dumps(res) + "\n"); f.flush()
            print(f"[{i+1}/{len(todo)}] {res['path'][:44]} {res['outcome']} certs={res['certs']} subboxes={res['subboxes']} {res['sec']:.0f}s (verified so far {nv})", flush=True)
    print(f"DONE verified {nv}/{len(todo)}, {time.time()-t0:.0f}s", flush=True)
