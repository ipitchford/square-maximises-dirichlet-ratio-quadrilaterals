"""Statistics of a cover log (JSONL from cover.jl): counts by status, certificate details, smallest margin."""
import json, sys, collections, statistics
f = sys.argv[1]
st = collections.Counter(); which = collections.Counter(); mode = collections.Counter(); ns = collections.Counter()
kap = []; ts = []; vols = []; worst = None; tF = []
for ln in open(f):
    r = json.loads(ln); s = r['status']; st[s] += 1
    if s == 'certified':
        which[r['which']] += 1; mode[r['mode']] += 1; ns[r['n']] += 1
        k = r['kappa_cert']; kap.append(k); ts.append(r['t']); tF.append(r.get('tF', 0.0))
        v = 1.0
        for i in range(4): v *= r['hi'][i] - r['lo'][i]
        vols.append(v)
        if worst is None or k > worst['kappa_cert']: worst = r
print("status counts:", dict(st))
if kap:
    print("certified: %d; test used: %s; mode: %s; degree n: %s" % (len(kap), dict(which), dict(mode), dict(ns)))
    print("certified bound on xi over a box: max %.6f (smallest margin 5/2 - max = %.6f), median %.4f, min %.4f" % (max(kap), 2.5 - max(kap), statistics.median(kap), min(kap)))
    print("box volume: median %.3e, min %.3e, max %.3e; total certified volume %.4e" % (statistics.median(vols), min(vols), max(vols), sum(vols)))
    print("time per certified box (1 thread): mean %.3f s (Float64 part %.3f s); sum %.0f s" % (statistics.mean(ts), statistics.mean(tF), sum(ts)))
    print("worst box:", worst['path'], "lo", worst['lo'], "hi", worst['hi'], "beta", worst['beta'], "which", worst['which'])
