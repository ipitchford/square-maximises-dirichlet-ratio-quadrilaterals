"""Rigorous P2 Rayleigh-Ritz upper bound for lambda_2 of ONE fixed dart Q_in = int(tri(B,D,A) u tri(B,D,C) u BD) (or of the single
triangle CBD when A is None), via the gcert P2 machinery at a point (zero box radius).  Domain monotonicity: if Q_in is contained
in every dart of a family (bo2.Arm.inner_vertex), then lambda_2(Q) <= lambda_2(Q_in) for every dart Q of the family.
p2_cert(A, C, E, k) certifies lambda_2(Q_in) < (5/2) E (Schur test of K2 - (5/2) E M2 < 0)."""
import numpy as np
from gcert import GCert, Problem
from charts import bd_problem, c_only_chart, bd_vfun
from p2up import solve_p2, p2_test

_GC = {}
def _gc(kind, k):
    key = (kind, k)
    if key not in _GC:
        if kind == 'bd': _GC[key] = GCert(bd_problem(k))
        else:
            ch = c_only_chart(k); _GC[key] = GCert(Problem(ch, ch, bd_vfun, name=f'Conly_k{k}'))
    return _GC[key]

def p2_cert(Ain, Cin, E, k=4):
    """Ain: (x, y) with y > 0 or None; Cin: (x, y) with y < 0.  Returns dict(ok, lam2p2, k, schur/why)."""
    if Ain is None:
        gc = _gc('c', k); p0 = np.array([1.0, 1.0, Cin[0], Cin[1]])   # RIGOR: exact: float64 arrays of Python floats (no rounding)
    else:
        gc = _gc('bd', k); p0 = np.array([Ain[0], Ain[1], Cin[0], Cin[1]])   # RIGOR: exact: float64 arrays of Python floats (no rounding)
    r = np.zeros(4)
    sol = solve_p2(gc, p0, r, order=0)
    lam2 = float(sol['lamp2'][1])   # RIGOR: heuristic: float pre-check only (it can only make the certificate fail)
    if not lam2 < 2.5*E: return dict(ok=False, why='float lam2 >= 2.5 E', lam2p2=lam2, k=k)
    res = p2_test(gc, p0, r, sol, E)
    res.update(lam2p2=lam2, k=k)
    return res
