"""FEM validation (non-rigorous reference) of certified leaves of the SPIKE / COLLAR / MAIN covers.
For each sampled certified leaf, sample darts Q covered by the leaf and check   L1_low(leaf) <= lambda1_FEM(Q)   (the certified
lower bound must lie below an accurate FEM value) and report xi_FEM(Q).
  MAIN  : p = (ax, ay, cx, cy) uniform in the box.
  SPIKE : t uniform in the box, C uniform in the box, a_x = X0 * exp(U[0, log 4])  (A = a_x (1, t)).
  COLLAR: t uniform in (0, T1], C uniform in the box, a_x in {0.6, 1.2, X0 * exp(U[0, log 4])}.
Usage: validate_fem2.py LOG REGION NLEAVES"""
import sys, json, numpy as np, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts')
from polyfem import poly_eigs
from regions import X0, T1
B = np.array([0., 0.]); D = np.array([1., 0.])
def fem(A, C, k=5):
    return poly_eigs([(D, [B, np.asarray(A, float)]), (D, [np.asarray(C, float), B])], k=k, kappa=2.5, deg=3)
log, region, nl = sys.argv[1], sys.argv[2], int(sys.argv[3])
rng = np.random.default_rng(1)
leaves = [json.loads(l) for l in open(log) if '"status": "cert"' in l]
pick = [leaves[i] for i in rng.choice(len(leaves), size=min(nl, len(leaves)), replace=False)]
worst = 1.0; maxxi = 0.0; viol = 0; n = 0
for L in pick:
    p0 = np.array(L['p0']); r = np.array(L['r'])
    for _ in range(2):
        u = p0 + r*rng.uniform(-1, 1, 4)
        if region == 'MAIN': A = (u[0], u[1]); C = (u[2], u[3])
        elif region == 'SPIKE':
            ax = X0*np.exp(rng.uniform(0, np.log(4))); A = (ax, ax*u[0]); C = (u[2], u[3])
        else:
            t = rng.uniform(0.2, 1.0)*T1; ax = [0.6, 1.2, X0*np.exp(rng.uniform(0, np.log(4)))][rng.integers(3)]
            A = (ax, ax*t); C = (u[2], u[3])
        w = fem(A, C); n += 1
        slack = (w[0] - L['L1'][0])/w[0]; worst = min(worst, slack); maxxi = max(maxxi, w[1]/w[0])
        if slack < 0: viol += 1; print('VIOLATION', L['p0'], L['r'], A, C, w[:2], L['L1'])
print(json.dumps(dict(log=log, region=region, darts_checked=n, violations=viol, min_slack=round(worst, 5), max_xi_fem=round(maxxi, 5))))
