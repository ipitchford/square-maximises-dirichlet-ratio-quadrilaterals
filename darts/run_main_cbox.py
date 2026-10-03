"""Rigorous BD-chart cover: C in CBOX = [0.40,0.60] x [-0.966,-0.766] and A with tan b_A >= 0.3, a_y <= AYMAX.
Boxes are discarded only by rigorous predicates: not a dart on the whole box; Prop 5.1; tan b_A <= 0.3 on the whole box
(covered by the COLLAR/SPIKE certificates for C in CBOX, run_crit.py); a_y > AYMAX on the whole box (outside this run).
Usage: run_main_cbox.py LOG AYMAX"""
import sys, json, time, numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts')
from cover2 import Cover2
from charts import bd_problem
from regions import fD_bounds, beta_max, tanb_bounds, in_cbox
log = sys.argv[1]; AYMAX = float(sys.argv[2])
CX = (0.40, 0.60); CY = (-0.966, -0.766)
def disc(p, r):
    lo, hi = fD_bounds(p, r)
    if lo >= 0: return 'not a dart (f_D >= 0 on box)'
    if beta_max(p, r) <= 0.0339: return 'Prop 5.1'
    if p[0] - r[0] > 0:
        tA = tanb_bounds(p[0] - r[0], p[0] + r[0], p[1] - r[1], p[1] + r[1])
        if tA[1] <= 0.3 and in_cbox(p, r): return 'tan b_A <= 0.3 (COLLAR/SPIKE certificates, C in CBOX)'
    if p[1] - r[1] >= AYMAX: return f'a_y >= {AYMAX!r} (outside this run)'
    return None
cv = Cover2(log, bd_problem, discard_fn=disc, kmin=4, kmax=6, maxdepth=30,
            split_scale_fn=lambda p: np.array([abs(p[1]), abs(p[1]), abs(p[3]), abs(p[3])]))
# A roots: a_x in [1, AYMAX/0.3 + 0.1], a_y in [0.25, AYMAX]; grid with cells ~ 10% of |A|
roots = []
ays = [0.25]
while ays[-1] < AYMAX: ays.append(min(ays[-1]*1.2, AYMAX))
for a in range(len(ays) - 1):
    y0, y1 = ays[a], ays[a+1]
    xmax = y1/0.3 + 0.05
    xs = [1.0]
    while xs[-1] < xmax: xs.append(xs[-1] + 0.2*max(1.0, y1))
    for b in range(len(xs) - 1):
        roots.append((f'M_{a}_{b}', [(xs[b]+xs[b+1])/2, (y0+y1)/2, (CX[0]+CX[1])/2, (CY[0]+CY[1])/2],
                      [(xs[b+1]-xs[b])/2, (y1-y0)/2, (CX[1]-CX[0])/2, (CY[1]-CY[0])/2]))
t0 = time.time()
for rid, p0, r in roots:
    ok = cv.process_root(rid, p0, r)
    print(json.dumps(dict(root=rid, ok=ok, stats=cv.stats, wall=round(time.time()-t0, 1))), flush=True)
print('DONE', json.dumps(cv.stats), round(time.time()-t0, 1))
