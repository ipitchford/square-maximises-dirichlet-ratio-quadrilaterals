"""Validation (non-rigorous reference): for certified boxes in a JSONL log, sample points p in the box, compute accurate FEM
eigenvalues (P3, graded mesh, fem_ref/polyfem) and check   L1_low(box) <= lambda1_FEM(p)   and   xi_FEM(p) < 5/2.
Also reports the gap (lambda1_FEM - L1_low)/lambda1_FEM (how much the certificate gives away)."""
import sys, json, numpy as np, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts')
from polyfem import poly_eigs
B = np.array([0., 0.]); D = np.array([1., 0.])
def fem(p, k=5):
    A = np.array([p[0], p[1]]); C = np.array([p[2], p[3]])
    w = poly_eigs([(D, [B, A]), (D, [C, B])], k=k, kappa=2.5, deg=3); return w
rng = np.random.default_rng(0)
viol = 0; n = 0; worst = 1.0; worst_xi = 0
for path in sys.argv[1:]:
    for line in open(path):
        d = json.loads(line)
        if d.get('status') != 'cert': continue
        p0 = np.array(d['p0']); r = np.array(d['r'])
        for _ in range(2):
            t = rng.uniform(-1, 1, 4); p = p0 + r*t
            w = fem(p); n += 1
            gap = (w[0] - d['L1'][0])/w[0]; worst = min(worst, gap); worst_xi = max(worst_xi, w[1]/w[0])
            if d['L1'][0] > w[0]: viol += 1; print('VIOLATION', path, p, d['L1'], w[:2])
print(f"checked {n} points: violations {viol}; min relative slack (lam1_FEM - L1_low)/lam1_FEM = {worst:.4f}; max xi_FEM = {worst_xi:.4f}")
