"""Rigorous Rayleigh-Ritz UPPER bound for lambda_2 over a box (conforming P2, the same machinery as step (4) of
gcert.GCert.certify, without the Crouzeix-Raviart / Temple part):  the P2 chart of a gcert.Problem (Dirichlet everywhere; a
domain contained in every Q of the box), two float test vectors transported to first order in the box parameters, 2x2 stiffness
and mass forms as order-3 Taylor models with rigorous remainders (exact rational reference matrices, outward rounding).
p2_test(gc, p0, r, sol, Lconst, tails) certifies  K2(p) - (5/2) min(Lconst, tails) M2(p) < 0  for all p in the box, i.e.
lambda_2(Q(p)) <= lambda_max(K2, M2) < (5/2) min(Lconst, tails) for all p in the box (Schur complement test, as in gcert)."""
import numpy as np
from flint import arb
from nanguard import fin, afin, amax, fup
from gcert import GCert, pvars, poly_to_tm, block_quad, vecquad2, mono_poly
from tm import TM
from schur import schur_upper

def solve_p2(gc, p0, r, order=1):   # RIGOR-FUNCTION: heuristic: float eigen-pairs / test vectors only
    """float P2 eigen-pairs and (order 1) first-order transported vectors (gc.solve without the CR part); order 0 for points"""
    p0 = np.asarray(p0, float); r = np.asarray(r, float)
    fg_p2 = gc._float_geo(gc.pb.p2c, p0, r)
    TKp = gc._taylor_mats(gc.pb.p2, fg_p2, 'p2')
    wp, rp = gc._vectors(TKp, [0, 1], 3, order)
    return dict(p0=p0, r=r, yvecs=[rp[0][0], rp[1][0]], lamp2=wp)

def p2_forms(gc, p0, r, sol):
    pb = gc.pb
    Z = pvars(np.asarray(p0, float), np.asarray(r, float))   # RIGOR: exact: float64 arrays of Python floats (no rounding)
    geoP = pb.geometry(pb.p2c, Z)
    GP = []
    for (N11, N12, N22, d, adR, E) in geoP:
        dT = poly_to_tm(d); lo, hi = dT.range()
        if lo > 0: sgn = 1
        elif hi < 0: sgn = -1
        else: raise ValueError('P2 piece degenerates on the box')
        dpos = d*sgn
        invd = poly_to_tm(dpos).recip()
        GP.append(dict(N=[N11, N12, N22], invd=invd, J=dpos*(1/adR)))
    ys = sol['yvecs']
    K2 = [[None]*2 for _ in range(2)]; M2 = [[None]*2 for _ in range(2)]
    for i in range(2):
        for j in range(2):
            ki = list(ys[i].keys()); kj = list(ys[j].keys())
            Yi = np.array([ys[i][k] for k in ki]).T; Yj = np.array([ys[j][k] for k in kj]).T   # RIGOR: exact: the float test vectors (free choice)
            mi = [mono_poly(k) for k in ki]; mj = [mono_poly(k) for k in kj]
            kt = TM.const(0); mt = TM.const(0)
            for q, g in enumerate(GP):
                bK = [block_quad(Yi, pb.p2['K'][q][kl], Yj, pb.p2['absK'][q][kl], pb.p2['mrow'][q][kl]) for kl in range(3)]
                kt = kt + g['invd']*poly_to_tm(vecquad2(bK, g['N'], mi, mj))
                bMq = block_quad(Yi, pb.p2['M'][q], Yj, pb.p2['absM'][q], pb.p2['mrowM'][q])
                mt = mt + poly_to_tm(vecquad2([bMq], [g['J']], mi, mj))
            K2[i][j] = kt; M2[i][j] = mt
    return K2, M2

def p2_test(gc, p0, r, sol, Lconst, tails=()):
    """True/False + diagnostics: K2 - (5/2) s M2 < 0 on the box for s = Lconst and every tail constant"""
    K2, M2 = p2_forms(gc, p0, r, sol)
    detM = M2[0][0]*M2[1][1] - M2[0][1]*M2[1][0]
    if not (detM.lower() > 0 and M2[0][0].lower() > 0): return dict(ok=False, why='M2 not PD')
    worsts = []
    for s_ in [arb(Lconst)] + [arb(c) for c in tails]:
        N = [[K2[i][j] - M2[i][j]*(s_*arb(2.5)) for j in range(2)] for i in range(2)]
        ok_, schur, why = schur_upper(N)
        if not ok_: return dict(ok=False, why=why, **({'schur': fup(schur)} if schur is not None else {}))
        worsts.append(schur)
    return dict(ok=True, schur=fup(amax(*worsts)))

def lam2_upper(gc, p0, r, sol, lo=0.0, hi=None, iters=40):   # RIGOR-FUNCTION: info: diagnostics only (not used by any certificate)
    """certified upper bound X of lambda_2 over the box (bisection on the Schur test K2 - X M2 < 0) -- diagnostics only"""
    K2, M2 = p2_forms(gc, p0, r, sol)
    hi = hi or 4*float(sol['lamp2'][1]) + 1
    def ok(X):
        N = [[K2[i][j] - M2[i][j]*arb(X) for j in range(2)] for i in range(2)]
        return schur_upper(N)[0]
    if not ok(hi): return None
    lo = float(sol['lamp2'][1])*0.5
    for _ in range(iters):
        mid = 0.5*(lo + hi)
        if ok(mid): hi = mid
        else: lo = mid
    return hi
