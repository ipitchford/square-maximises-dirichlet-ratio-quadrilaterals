"""Independent re-certification of a stratified random sample of transport-certified boxes with xcert.certify.
A box counts as VERIFIED if xcert certifies it directly or all children of a bisection (depth <= MAXD) certify.
Usage: xsample.py <jsonl> <n_sample> <out.jsonl> [seed]"""
import json, sys, random, math, time
from fractions import Fraction as Fr
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from xcert import certify
MAXD = 4
def exact(v): return Fr(v)            # Float64 box ends are dyadic: Fraction(float) is exact
def try_box(lo, hi, depth):
    ok, res, beta = certify(lo, hi)
    if ok: return True, 1, [res]
    if depth >= MAXD: return False, 1, [res]
    w = [h-l for l, h in zip(lo, hi)]; d = w.index(max(w)); mid = (lo[d]+hi[d])/2
    hi1 = list(hi); hi1[d] = mid; lo2 = list(lo); lo2[d] = mid
    a = try_box(lo, hi1, depth+1); b = try_box(lo2, hi, depth+1)
    return a[0] and b[0], a[1]+b[1], a[2]+b[2]
if __name__ == "__main__":
    src, n, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]; seed = int(sys.argv[4]) if len(sys.argv) > 4 else 7
    recs = [r for r in map(json.loads, open(src)) if r.get("status") == "certified"]
    strata = {}
    for r in recs:
        vol = math.prod(h-l for l, h in zip(r["lo"], r["hi"]))
        key = (r["path"].split(".")[0], int(math.floor(math.log10(vol))))
        strata.setdefault(key, []).append(r)
    rng = random.Random(seed); keys = sorted(strata)
    picks = []
    while len(picks) < n and any(strata[k] for k in keys):
        for k in keys:
            if strata[k] and len(picks) < n:
                picks.append(strata[k].pop(rng.randrange(len(strata[k]))))
    print(f"{len(recs)} certified records, {len(keys)} strata, sampling {len(picks)}", flush=True)
    nv = nf = 0; t0 = time.time()
    with open(out, "w") as f:
        for i, r in enumerate(picks):
            lo = [exact(v) for v in r["lo"]]; hi = [exact(v) for v in r["hi"]]
            t = time.time(); ok, nboxes, _ = try_box(lo, hi, 0)
            nv += ok; nf += (not ok)
            f.write(json.dumps({"path": r["path"], "lo": r["lo"], "hi": r["hi"], "verified": ok, "xboxes": nboxes,
                                "agent_mode": r.get("mode"), "agent_kappa": r.get("kappa_cert"), "sec": round(time.time()-t, 1)}) + "\n"); f.flush()
            if (i+1) % 10 == 0: print(f"[{i+1}/{len(picks)}] verified {nv}, not verified {nf}, {time.time()-t0:.0f}s", flush=True)
    print(f"DONE verified {nv}, not verified {nf}, {time.time()-t0:.0f}s", flush=True)
