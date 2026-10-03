"""Rigorous cover of the critical spike neighbourhood (where xi -> 7/3):
   C in CBOX = [0.40, 0.60] x [-0.966, -0.766] (CBD near equilateral), A a spike/sliver:
     COLLAR : tan b_A <= 0.02 (every A on such a ray, any length)           chart SPIKE (t = 0.02 fixed), params (cx, cy)
     SPIKE  : 0.02 <= tan b_A <= 0.30 and a_x >= 1.5 (any length)            chart SPIKE2 (X0 = 1.5), params (t, cx, cy)
Usage: run_crit.py {COLLAR|SPIKE} LOG"""
import sys, numpy as np, json, time
sys.path.insert(0, '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts')
from cover2 import Cover2
from charts import spike_problem, spike2_problem
from regions import X0, T1, T2
CX = (0.40, 0.60); CY = (-0.966, -0.766)
which, log = sys.argv[1], sys.argv[2]
nC = 2   # C roots per axis
cxs = np.linspace(CX[0], CX[1], nC + 1); cys = np.linspace(CY[0], CY[1], nC + 1)
roots = []
if which == 'COLLAR':
    cv = Cover2(log, spike_problem, kmin=5, kmax=6, maxdepth=22, split_scale_fn=lambda p: np.array([1.0, 1.0, abs(p[3]), abs(p[3])]))
    for i in range(nC):
        for j in range(nC):
            roots.append((f'COL_{i}{j}', [T1, 0.0, (cxs[i]+cxs[i+1])/2, (cys[j]+cys[j+1])/2], [0.0, 0.0, (cxs[1]-cxs[0])/2, (cys[1]-cys[0])/2]))
else:
    cv = Cover2(log, lambda k: spike2_problem(k, X0), kmin=5, kmax=6, maxdepth=26,
                split_scale_fn=lambda p: np.array([abs(p[0]), 1.0, abs(p[3]), abs(p[3])]))
    # t slabs: geometric, ratio 1.25
    ts = [T1]
    while ts[-1] < T2: ts.append(min(ts[-1]*1.25, T2))
    for a in range(len(ts) - 1):
        for i in range(nC):
            for j in range(nC):
                roots.append((f'SPK_{a}_{i}{j}', [(ts[a]+ts[a+1])/2, 0.0, (cxs[i]+cxs[i+1])/2, (cys[j]+cys[j+1])/2],
                              [(ts[a+1]-ts[a])/2, 0.0, (cxs[1]-cxs[0])/2, (cys[1]-cys[0])/2]))
t0 = time.time()
for rid, p0, r in roots:
    ok = cv.process_root(rid, p0, r)
    print(json.dumps(dict(root=rid, ok=ok, stats=cv.stats, wall=round(time.time()-t0, 1))), flush=True)
print('DONE', json.dumps(cv.stats), round(time.time()-t0, 1))
