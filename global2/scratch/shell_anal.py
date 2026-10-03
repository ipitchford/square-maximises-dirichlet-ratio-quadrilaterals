import json, sys, collections, math
f = sys.argv[1]
passes = []
for ln in open(f):
    if 'pilot_pass' not in ln and '"certified"' not in ln: continue
    r = json.loads(ln); lo, hi = r['lo'], r['hi']
    vol = 1.0
    for i in range(4): vol *= hi[i]-lo[i]
    pinf = max(max(abs(lo[i]), abs(hi[i])) for i in range(4))
    passes.append((pinf, vol, r.get('kappa_f64', r.get('kappa_cert'))))
edges = [0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25]
def shellvol(r):  # volume of F ∩ {|p|inf <= r} = r^3 * 2r
    return 2*r**4
prev = 0.049
cum_n = 0; cum_v = 0
for e in edges:
    sel = [p for p in passes if prev < p[0] <= e + 1e-12]
    n = len(sel); v = sum(p[1] for p in sel)
    sv = shellvol(e) - shellvol(prev)
    cum_n += n; cum_v += v
    print(f"|p|inf in ({prev:.3f},{e:.3f}]: boxes {n:6d}  covered vol {v:.3e} of shell {sv:.3e} ({100*v/sv:5.1f}%)  cum boxes {cum_n}")
    prev = e
print("total passes", len(passes), "max kappa", max(p[2] for p in passes))
