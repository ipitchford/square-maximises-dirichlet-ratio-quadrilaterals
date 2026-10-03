"""Stage 3b: Schur complement of the corrected 3x3 pencil as multivariate Taylor models (before any direction substitution).
Exact structural facts used (lemma_firstorder.py, lemma_residual.py; they hold for the corrected modes because the
corrections vanish at p=0 and only couple a base mode to modes orthogonal to it at p=0):
  (L0) K(0)=diag(l1,l2,l2), M(0)=I, LL11(0)=l1^2      -> zero constants of off-diagonal K, M; eps2(0)=0; A_cc(0)=0
  (L1) d rho(0) = 0 and d tr A_cc(0) = 0               -> zero degree-1 of rho and of tr S
  (L2) d eps2(0) = 0                                    -> zero degree-1 of eps2
Each zeroed coefficient is asserted to contain 0 in its computed enclosure (consistency check)."""
import sys, time, pickle
from fractions import Fraction as Fr
from flint import arb
from localcert import load, corrections, build_pencil, truncate_to, fr2arb
from mtm import MTM, ZERO4, deg

BASE = [(1,1),(2,1),(1,2)]

def ser(t): return ({g: (v.mid().str(60), v.rad().str(10)) for g, v in t.c.items()}, (t.r.mid().str(60), t.r.rad().str(10)), t.m)
def deser(s, rho):
    cd, r, m = s
    ball = lambda x: arb(x[0]) + arb(x[1]) * arb("[+/- 1]")
    return MTM({g: ball(v) for g, v in cd.items()}, ball(r), m, rho)

def build(fn, betas, out=None, verbose=True):
    t0 = time.time()
    D = load(fn); m, mLL, rho0 = D['m'], D['mLL'], D['rho0']; rho = fr2arb(rho0)
    C = corrections(D); Kt, Mt, LLt = build_pencil(D, C)
    if verbose: print(f"[{time.time()-t0:.0f}s] pencil built (N={D['N']}, m={m}, mLL={mLL}, rho0={float(rho0)})", flush=True)
    K = [[MTM(*Kt[(BASE[i], BASE[j])], m, rho) for j in range(3)] for i in range(3)]
    M = [[MTM(*Mt[(BASE[i], BASE[j])], m, rho) for j in range(3)] for i in range(3)]
    cd, r, _ = truncate_to(*LLt[(BASE[0], BASE[0])], mLL, m, rho0); LL = MTM(cd, r, m, rho)
    for i in range(3):
        for j in range(3):
            if i != j: K[i][j].zero_degree(0); M[i][j].zero_degree(0)
    # M > 0 on the ball (Gershgorin)
    for i in range(3):
        rest = MTM({g: v for g, v in M[i][i].c.items() if g != ZERO4}, M[i][i].r, m, rho)
        low = M[i][i].c[ZERO4] - rest.sup() - sum((M[i][j].sup() for j in range(3) if j != i), arb(0))
        assert low > 0, f"Gershgorin M row {i}: {low}"
    Minv = M[0][0].inv(); rq = K[0][0] * Minv; rq.zero_degree(1)
    eps2 = LL * Minv - rq * rq; eps2.zero_degree(0); eps2.zero_degree(1)
    if verbose: print(f"[{time.time()-t0:.0f}s] rho, eps2 done; eps2 rem {eps2.r.upper()}, sup {eps2.sup().upper()}", flush=True)
    res = {}
    for beta in betas:
        Lt = rq - eps2 * (MTM.const(fr2arb(beta), m, rho) - rq).inv()
        A = [[K[i][j] - (Lt * M[i][j]).scale(arb(5)/2) for j in range(3)] for i in range(3)]
        for i in range(3):
            for j in range(3):
                if (i, j) != (0, 0) and ZERO4 in A[i][j].c: A[i][j].zero_degree(0)
        A00inv = A[0][0].inv()    # asserts A00 keeps the sign of -3 pi^2 on the ball
        S = {}
        for (i, j) in [(0,0), (0,1), (1,1)]:
            S[(i,j)] = A[1+i][1+j] - A[1+i][0] * A[0][1+j] * A00inv
        trS = S[(0,0)] + S[(1,1)]
        assert all(deg(g) >= 1 for g in trS.c); trS.zero_degree(1)
        detS = S[(0,0)] * S[(1,1)] - S[(0,1)] * S[(0,1)]
        assert all(deg(g) >= 2 for g in detS.c)
        a00 = A[0][0].c[ZERO4]; a00sup = MTM({g: v for g, v in A[0][0].c.items() if g != ZERO4}, A[0][0].r, m, rho).sup()
        assert a00 + a00sup < 0, "A00 not certified negative"
        res[str(beta)] = dict(detS=ser(detS), trS=ser(trS))
        if verbose:
            q = lambda t, g: float(t.c.get(g, arb(0)).mid())
            print(f"[{time.time()-t0:.0f}s] beta={float(beta)}: detS r={float(detS.r.upper()):.3g} trS r={float(trS.r.upper()):.3g}; "
                  f"trS[b^2]={q(trS,(0,2,0,0)):.4f} trS[d^2]={q(trS,(0,0,0,2)):.3f} detS[a^2]={q(detS,(2,0,0,0)):.2f} detS[d^2]={q(detS,(0,0,0,2)):.2f}", flush=True)
    meta = dict(fn=fn, N=D['N'], m=m, mLL=mLL, rho0=str(rho0))
    if out: pickle.dump(dict(meta=meta, res=res), open(out, 'wb'))
    return meta, res

if __name__ == "__main__":
    fn = sys.argv[1]; out = sys.argv[2]; betas = [Fr(x) for x in sys.argv[3].split(',')]
    build(fn, betas, out)
