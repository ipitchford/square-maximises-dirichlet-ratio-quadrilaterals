"""Independent Monte Carlo coverage check of the dart decomposition against the final logs (receipt with provenance hashes)."""
import sys, random, time, json, hashlib, collections, math
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import xdart as X
h = lambda fn: hashlib.sha256(open(fn, "rb").read()).hexdigest()
prov = {k: h(X.DARTS + rel) for k, rel in X.LOGS.items()}; prov["xdart.py"] = h(__file__.rsplit("/", 1)[0] + "/xdart.py"); prov["xdart_run.py"] = h(__file__)
print("PROVENANCE", json.dumps(prov), flush=True)
ix = X.load_indices(); print("CERT leaves:", {k: len(v.p0) for k, v in ix.items()}, flush=True)
rng = random.Random(20261002); res = collections.Counter(); unc = []; reasons = collections.Counter(); t0 = time.time()
def run(name, gen, n):
    for _ in range(n):
        d = gen(); r = X.covered(d, ix)
        res[name + (":covered" if r else ":UNCOVERED")] += 1
        reasons[(r or "NONE").split(" ")[0] if not (r or "").startswith("lemma") else " ".join(r.split()[:2])] += 1
        if r is None: unc.append((name, d))
run("BD-uniform(logr -6..5)", lambda: X.sample(rng, (-6.0, 5.0)), 50000)
run("BD-wide(logr -8..8)", lambda: X.sample(rng, (-8.0, 8.0)), 20000)
run("hull-uniform", lambda: X.sample_hull(rng), 50000)
run("huge-radius(logr -8..40)", lambda: X.sample(rng, (-8.0, 40.0)), 20000)
for name, d in X.boundary_samples(rng, 2000):
    r = X.covered(d, ix); res[name + (":covered" if r else ":UNCOVERED")] += 1
    if r is None: unc.append((name, d))
print("RESULTS", json.dumps(dict(sorted(res.items()))), flush=True)
print("covered-by (first reason):", dict(reasons.most_common(20)), flush=True)
print(f"UNCOVERED total: {len(unc)}  ({time.time()-t0:.0f}s)", flush=True)
for name, d in unc[:30]:
    bA, rA, bC, rC = d
    print("  UNCOVERED %-22s bA %.4f deg rA %.5g bC %.4f deg rC %.5g t=%.5f a_x=%.5f thetaB=%.4f deg beta=%.4f"
          % (name, math.degrees(bA), rA, math.degrees(bC), rC, math.tan(bA), rA*math.cos(bA), math.degrees(bA+bC), X.beta(*d)), flush=True)
print("XDART", "PASS" if not unc else "FAIL")
sys.exit(0 if not unc else 1)
