"""FEM validation (non-rigorous reference) of cover5 certificates: for sampled certified leaves, sample darts of the family box
and compare with accurate FEM (P3, graded, polyfem):  BO2 leaves: E <= lambda_1(Q) and lambda_2(Q) <= lam2p2 (P2 bound of the
inner dart, which bounds lambda_2(Q) by domain monotonicity);  INC leaves: the same with E = L1 of the outer dart;  FEM leaves: L1_low <= lambda_1(Q).  Reports violations, the minimal
relative slack and the largest FEM xi.  Far arms are sampled at rho in [RF, 3 RF] (FEM of a finite dart).
Usage: validate_fem5.py LOG FAMILY NLEAVES [SEED]"""
import sys, json, math, numpy as np, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts')
import fam5
from polyfem import poly_eigs
B = np.array([0., 0.]); D = np.array([1., 0.])
def fem(A, C, k=5):
    return poly_eigs([(D, [B, np.asarray(A, float)]), (D, [np.asarray(C, float), B])], k=k, kappa=2.5, deg=3)
def sample(fam, p0, r, rng):
    u = p0 + r*rng.uniform(-1, 1, 4)
    name = fam.name
    if name.startswith('SP') or name.startswith('SF'):
        t = max(u[0], 1e-4); amin = 1.0 if name.endswith('spk') else 0.3
        a = amin*math.exp(rng.uniform(0, math.log(6)))
        A = (a, a*t)
    elif getattr(fam, 'Afar', False):
        rho = fam5.RF*math.exp(rng.uniform(0, math.log(3))); A = (rho*math.cos(u[0]), rho*math.sin(u[0]))
    else:
        rho = math.exp(u[1]); A = (rho*math.cos(u[0]), rho*math.sin(u[0]))
    if getattr(fam, 'Cfar', False) or name.startswith('SF'):
        rho = fam5.RF*math.exp(rng.uniform(0, math.log(3))); C = (rho*math.cos(u[2]), -rho*math.sin(u[2]))
    else:
        rho = math.exp(u[3]); C = (rho*math.cos(u[2]), -rho*math.sin(u[2]))
    return A, C
if __name__ == '__main__':
    log, name, nl = sys.argv[1], sys.argv[2], int(sys.argv[3]); seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    fam = fam5.FAMILIES[name](); rng = np.random.default_rng(seed)
    leaves = [json.loads(l) for l in open(log) if '"status": "cert"' in l]
    pick = [leaves[i] for i in rng.choice(len(leaves), size=min(nl, len(leaves)), replace=False)]
    viol = 0; n = 0; worst1 = 1.0; worst2 = 1.0; maxxi = 0.0
    for L in pick:
        p0 = np.array(L['p0']); r = np.array(L['r'])
        for _ in range(2):
            A, C = sample(fam, p0, r, rng)
            fD = (1 - C[0])*A[1] + C[1]*(A[0] - 1)
            if fD >= 0: continue                       # validate on darts only
            try: w = fem(A, C)
            except Exception as e: continue
            n += 1; maxxi = max(maxxi, w[1]/w[0])
            lo = L['E'] if L['method'] in ('BO2', 'INC') else L['L1'][0]
            s1 = (w[0] - lo)/w[0]; worst1 = min(worst1, s1)
            if s1 < 0: viol += 1; print('VIOLATION lambda1', L['root'], L['path'], A, C, w[:2], lo)
            if L['method'] in ('BO2', 'INC'):
                s2 = (L['lam2p2'] - w[1])/w[1]; worst2 = min(worst2, s2)
                if s2 < -1e-3: viol += 1; print('VIOLATION lambda2', L['root'], L['path'], A, C, w[:2], L['lam2p2'])
    print(json.dumps(dict(log=log, family=name, darts=n, violations=viol, min_slack_lambda1=round(worst1, 4),
                          min_slack_lambda2=round(worst2, 4), max_xi_fem=round(maxxi, 4))))
