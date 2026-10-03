"""Rigorous BD-chart cover of an axis box in (ax, ay, cx, cy), split into roots of relative size ~REL.
Discards (rigorous): not a dart on the whole box; Prop 5.1 (beta <= 0.0339 on the whole box).
Usage: run_box.py LOG ax0 ax1 ay0 ay1 cx0 cx1 cy0 cy1 NROOT_PER_AXIS"""
import sys, json, time, numpy as np, itertools
sys.path.insert(0, '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts')
from cover2 import Cover2
from charts import bd_problem
from regions import fD_bounds, beta_max
log = sys.argv[1]; lim = [float(x) for x in sys.argv[2:10]]; nr = int(sys.argv[10])
def disc(p, r):
    lo, hi = fD_bounds(p, r)
    if lo >= 0: return 'not a dart (f_D >= 0 on box)'
    if beta_max(p, r) <= 0.0339: return 'Prop 5.1'
    return None
cv = Cover2(log, bd_problem, discard_fn=disc, kmin=4, kmax=6, maxdepth=32,
            split_scale_fn=lambda p: np.array([abs(p[1]), abs(p[1]), abs(p[3]), abs(p[3])]))
edges = [np.linspace(lim[2*i], lim[2*i+1], nr + 1) for i in range(4)]
t0 = time.time()
for idx in itertools.product(range(nr), repeat=4):
    p0 = [(edges[i][idx[i]] + edges[i][idx[i]+1])/2 for i in range(4)]
    r = [(edges[i][1] - edges[i][0])/2 for i in range(4)]
    rid = 'B_' + '_'.join(map(str, idx))
    ok = cv.process_root(rid, p0, r)
    print(json.dumps(dict(root=rid, ok=ok, stats=cv.stats, wall=round(time.time()-t0, 1))), flush=True)
print('DONE', json.dumps(cv.stats), round(time.time()-t0, 1))
