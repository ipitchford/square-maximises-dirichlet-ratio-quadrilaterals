"""Compare the Float64 certificate pipeline (validation_points.csv, written by validate_points.jl) with
P2 FEM (../explore/quadfem.py, k=6).  Checks: L_LG <= FEM lambda1 (soft: FEM is an upper bound), FEM lambda1 <= theta1,
theta2 >= FEM lambda2 (both RR and conforming FEM are upper bounds; theta2 >= true lambda2 is what matters)."""
import sys, csv, numpy as np
sys.path.insert(0, '../explore')
from quadfem import eigs
rows = list(csv.DictReader(open('validation_points.csv')))
out = open('validation_fem.log', 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
P(f"{'tag':8s} {'n':>2s} {'maxang':>6s} {'w/D':>5s} {'rmin':>6s} {'L_LG':>10s} {'FEM l1':>10s} {'theta1':>10s} {'FEM l2':>10s} {'theta2':>10s} {'beta':>8s} {'xi_bound':>8s} {'FEM xi':>7s} ok")
cache = {}
nbad = 0; npts = 0
for r in rows:
    V = [(float(r[f'v{j}x']), float(r[f'v{j}y'])) for j in range(1,5)]
    key = tuple(V)
    if key not in cache:
        cache[key] = eigs(V, k=6, nev=3)
    l = cache[key]
    L = float(r['LG_pt']); t1 = float(r['theta1']); t2 = float(r['theta2'])
    ok_l = (np.isnan(L) or L <= l[0] * (1 + 1e-7))
    ok_t1 = l[0] <= t1 * (1 + 1e-7)
    ok_t2 = t2 >= l[1] * (1 - 1e-7)
    ok = ok_l and ok_t1 and ok_t2
    nbad += (not ok); npts += 1
    xb = t2 / L if L == L and L > 0 else float('nan')
    P(f"{r['tag']:8s} {int(r['n']):2d} {float(r['maxangle']):6.1f} {float(r['wd']):5.3f} {float(r['ratio_min']):6.3f} {L:10.5f} {l[0]:10.5f} {t1:10.5f} {l[1]:10.5f} {t2:10.5f} {float(r['beta_pt']):8.3f} {xb:8.5f} {l[1]/l[0]:7.5f} {'ok' if ok else 'VIOLATION l1 %s t1 %s t2 %s' % (ok_l, ok_t1, ok_t2)}")
P(f"rows {npts}, distinct points {len(cache)}, violations {nbad}")
