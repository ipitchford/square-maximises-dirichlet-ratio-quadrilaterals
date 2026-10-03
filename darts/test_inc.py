"""Numerical validation of the INC certificate machinery (inc.py, LEMMAS Lemma 11.4) -- NOT part of the proof.
 (1) d_range encloses the angle at D over random polar boxes (dense sampling);
 (2) for random PP boxes, sampled darts of the box contain the inner dart and are contained in the outer dart (polygon test);
 (3) for INC leaves of the pool logs (sample), sampled darts of the box: E <= lambda_1^{P2}(Q) (P2 is an upper bound of lambda_1,
     so a violation would expose a wrong E) and Q- subset Q subset Q+.
Usage: test_inc.py [N_LEAVES]"""
import sys, os, glob, json, math, numpy as np, warnings
warnings.filterwarnings('ignore')
HERE = '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts'
sys.path.insert(0, HERE)
import inc
from p2fix import _gc
from p2up import solve_p2
rng = np.random.default_rng(20261001)

def in_tri(P, X, tol=1e-12):
    """P in the closed triangle B=(0,0), D=(1,0), X (X in the upper half plane)"""
    return P[1] >= -tol and X[0]*P[1] - X[1]*P[0] <= tol and (X[0] - 1)*P[1] - X[1]*(P[0] - 1) >= -tol

def dart_in(Ain, Cin, A, C):
    """dart(Ain, Cin) subset dart(A, C): the apexes of the smaller triangles lie in the larger ones"""
    return in_tri(Ain, A) and in_tri((Cin[0], -Cin[1]), (C[0], -C[1]))

bad = 0; n1 = 0
for _ in range(400):
    b = rng.uniform(0.02, 3.0); l = rng.uniform(-2.5, 3.0); rb = rng.uniform(1e-4, 0.1); rl = rng.uniform(1e-4, 0.3)
    if b - rb <= 0 or b + rb >= math.pi: continue
    dr = inc.d_range(b - rb, b + rb, l - rl, l + rl)
    if dr is None: continue
    B, L = np.meshgrid(np.linspace(b - rb, b + rb, 201), np.linspace(l - rl, l + rl, 41)); R = np.exp(L)
    d = np.arctan2(R*np.sin(B), 1 - R*np.cos(B)); n1 += 1
    if d.min() < dr[0] - 1e-12 or d.max() > dr[1] + 1e-12: bad += 1
print(f'(1) d_range: {n1} boxes, {bad} enclosure violations')

bad = 0; n2 = 0
for _ in range(300):
    p0 = np.array([rng.uniform(0.05, 1.7), rng.uniform(-2.5, 2.0), rng.uniform(0.05, 1.7), rng.uniform(-2.5, 2.0)])
    r = np.array([rng.uniform(1e-3, 0.05), rng.uniform(1e-3, 0.2), rng.uniform(1e-3, 0.05), rng.uniform(1e-3, 0.2)])
    dd = inc.inc_darts(p0, r)
    if dd is None: continue
    Ao, Co, Ai, Ci = dd; n2 += 1
    for _ in range(40):
        q = p0 + r*rng.uniform(-1, 1, 4)
        A = (math.exp(q[1])*math.cos(q[0]), math.exp(q[1])*math.sin(q[0])); C = (math.exp(q[3])*math.cos(q[2]), -math.exp(q[3])*math.sin(q[2]))
        if not (dart_in(Ai, Ci, A, C) and dart_in(A, C, Ao, Co)): bad += 1
print(f'(2) inclusion: {n2} boxes x 40 darts, {bad} violations')

nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 40
leaves = []
for f in sorted(glob.glob(os.path.join(HERE, 'certs3', 'pp.*.jsonl'))):
    for line in open(f):
        if '"INC"' in line:
            d = json.loads(line)
            if d.get('method') == 'INC': leaves.append(d)
bad = 0; n3 = 0; worst = float('inf')
gc = _gc('bd', 5)
for d in [leaves[i] for i in rng.choice(len(leaves), min(nmax, len(leaves)), replace=False)] if leaves else []:
    p0 = np.array(d['p0']); r = np.array(d['r'])
    for _ in range(3):
        q = p0 + r*rng.uniform(-1, 1, 4)
        A = (math.exp(q[1])*math.cos(q[0]), math.exp(q[1])*math.sin(q[0])); C = (math.exp(q[3])*math.cos(q[2]), -math.exp(q[3])*math.sin(q[2]))
        l1 = float(solve_p2(gc, np.array([A[0], A[1], C[0], C[1]]), np.zeros(4), order=0)['lamp2'][0])
        ok = l1 >= d['E'] and dart_in(d['Ain'], d['Cin'], A, C) and dart_in(A, C, d['Aout'], d['Cout'])
        worst = min(worst, l1/d['E']); n3 += 1; bad += (not ok)
print(f'(3) INC leaves: {len(leaves)} in the logs, {n3} sampled darts, {bad} violations, min lambda1^P2/E = {worst:.4f}')
sys.exit(1 if bad else 0)
