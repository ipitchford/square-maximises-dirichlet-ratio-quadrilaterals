"""Inclusion certificate INC (LEMMAS Lemma 11.4) for a polar x polar box of darts (family PP, BD chart).

For a box  b_A in [bA1, bA2], ln rho_A in [lA1, lA2], b_C in [bC1, bC2], ln rho_C in [lC1, lC2]:
  * d(X) = angle at D = (1, 0) of the triangle BDX.  For X = rho (cos b, sin b), 0 < b < pi:
        d = atan2(rho sin b, 1 - rho cos b),  dd/drho = sin b / |X - D|^2 > 0,  dd/db = rho (cos b - rho) / |X - D|^2,
    so d increases in rho, and in b it is unimodal (increasing while cos b > rho, then decreasing; decreasing if rho >= 1).
    Hence  max d over the box = max over b in [b1, b2] of d(b, rho2)  (attained at b1, b2 or at b* = arccos rho2), and
           min d over the box = min(d(b1, rho1), d(b2, rho1)).  Both are evaluated in arb (d_range).
  * OUTER apex X+ with angle at B >= b2 and angle at D >= max d:  every arm vertex X of the box lies in the closed triangle
    BDX+ (a triangle on BD is the intersection of its angular sectors at B and D), so tri(B,D,X) is contained in tri(B,D,X+).
    INNER apex X- with angle at B <= b1 and angle at D <= min d:  X- lies in tri(B,D,X), so tri(B,D,X-) is contained in it.
  * Q- := dart(A-, C-)  is contained in every dart Q of the box, which is contained in  Q+ := dart(A+, C+);  hence
        lambda_2(Q) <= lambda_2(Q-) <= P2 Rayleigh-Ritz bound   and   lambda_1(Q) >= lambda_1(Q+) >= L1 (gcert at the point Q+).
    p2fix.p2_cert(A-, C-, E, k) with E = L1 certifies lambda_2(Q-) < (5/2) E, hence xi(Q) < 5/2 on the whole box.
Apex coordinates are plain floats (stored in the log); the replay re-verifies the angle inequalities in arb and recomputes
L1 and the P2 test from scratch."""
import math, numpy as np
from flint import arb
from nanguard import afin, fdown, fup, amax, amin, exact_box
from p2fix import _gc, p2_cert

RPT = 2.0**-30      # half-width of the tiny Cartesian box around Q+ used for the L1 certificate (contains the point Q+)
RMAX = 8.0          # apex radius cap (the BD-chart FEM / Liu bound is useless for very long triangles), same as fam5.FEM_RMAX

def _ival(lo, hi): return arb(lo).union(arb(hi))

def d_arb(b, rho):
    """angle at D of rho (cos b, sin b) (arb in, arb out; b in (0, pi))"""
    return arb.atan2(rho*b.sin(), 1 - rho*b.cos())

def d_range(b1, b2, l1, l2):
    """rigorous floats (dlo, dhi): dlo <= d <= dhi for every X = e^l (cos b, sin b), b in [b1, b2], l in [l1, l2]
    (0 < b1 <= b2 < pi).  Uses the monotonicity / unimodality stated in the module docstring."""
    if not (0 < b1 <= b2 < math.pi): return None
    B1, B2 = arb(b1), arb(b2)
    rho1 = arb(l1).exp(); rho2 = arb(l2).exp()
    cand = [d_arb(B1, rho2), d_arb(B2, rho2)]
    if rho2.lower() < 1:                                    # interior peak at b* = arccos(rho2) possible
        bstar = rho2.acos() if rho2.upper() < 1 else None
        if bstar is None or not bstar.is_finite():
            cand.append(d_arb(_ival(b1, b2), rho2))         # fallback: enclosure over the whole b-interval (a superset)
        elif bstar.upper() >= b1 and bstar.lower() <= b2:
            cand.append(d_arb(bstar, rho2))
    lows = [d_arb(B1, rho1), d_arb(B2, rho1)]
    if not all(c.is_finite() for c in cand + lows): return None
    return fdown(amin(*lows)), fup(amax(*cand))              # directed: lower end down, upper end up

def _apex(beta, delta):
    L = math.sin(delta)/math.sin(beta + delta)
    return (L*math.cos(beta), L*math.sin(beta))

def angB(X): return arb.atan2(arb(X[1]), arb(X[0]))
def angD(X): return arb.atan2(arb(X[1]), 1 - arb(X[0]))

def check_outer(X, bmax, dmax):
    """rigorous: X upper half-plane, angle at B >= bmax and angle at D >= dmax (so tri(B,D,X) contains all arm vertices)"""
    return bool(X[1] > 0 and angB(X).lower() >= bmax and angD(X).lower() >= dmax and (angB(X) + angD(X)).upper() < arb.pi().lower())

def check_inner(X, bmin, dmin):
    """rigorous: X upper half-plane, angle at B <= bmin and angle at D <= dmin (so X lies in every tri(B,D,arm vertex))"""
    return bool(X[1] > 0 and angB(X).upper() <= bmin and angD(X).upper() <= dmin)

def outer_apex(bmax, dmax):
    if not (bmax + dmax < math.pi - 1e-9): return None
    for pad in (1e-13, 1e-12, 1e-11, 1e-10, 1e-9):
        X = _apex(bmax + pad, dmax + pad)
        if check_outer(X, bmax, dmax): return X
    return None

def inner_apex(bmin, dmin):
    for pad in (1e-13, 1e-12, 1e-11, 1e-10, 1e-9):
        if bmin - pad <= 0 or dmin - pad <= 0: return None
        X = _apex(bmin - pad, dmin - pad)
        if check_inner(X, bmin, dmin): return X
    return None

def arm_ranges(p0, r):
    """((bA1, bA2, dA_lo, dA_hi), (bC1, bC2, dC_lo, dC_hi)) for a PP box (p0, r) = (b_A, l_A, b_C, l_C) +- r, or None"""
    exact_box(p0, r)
    out = []
    for ib, il in ((0, 1), (2, 3)):
        b1, b2 = p0[ib] - r[ib], p0[ib] + r[ib]
        dr = d_range(b1, b2, p0[il] - r[il], p0[il] + r[il])
        if dr is None: return None
        out.append((b1, b2, dr[0], dr[1]))
    return out

def inc_darts(p0, r):
    """(Aout, Cout, Ain, Cin) as float pairs in the BD chart (C with y < 0), or None"""
    ar = arm_ranges(p0, r)
    if ar is None: return None
    (bA1, bA2, dA1, dA2), (bC1, bC2, dC1, dC2) = ar
    Ao = outer_apex(bA2, dA2); Co = outer_apex(bC2, dC2); Ai = inner_apex(bA1, dA1); Ci = inner_apex(bC1, dC1)
    if None in (Ao, Co, Ai, Ci): return None
    if max(math.hypot(*Ao), math.hypot(*Co)) > RMAX: return None   # RIGOR: heuristic: routing (very long apexes are not tried)
    return Ao, (Co[0], -Co[1]), Ai, (Ci[0], -Ci[1])

def verify_inclusion(p0, r, Aout, Cout, Ain, Cin):
    """rigorous re-check (replay): the recorded apexes satisfy the outer / inner angle conditions for the box"""
    ar = arm_ranges(p0, r)
    if ar is None: return False
    (bA1, bA2, dA1, dA2), (bC1, bC2, dC1, dC2) = ar
    mC = lambda X: (X[0], -X[1])
    return (check_outer(Aout, bA2, dA2) and check_outer(mC(Cout), bC2, dC2) and check_inner(Ain, bA1, dA1)
            and check_inner(mC(Cin), bC1, dC1) and max(math.hypot(*Aout), math.hypot(*Cout)) <= RMAX)   # RIGOR: heuristic: consistency check only; soundness of INC does not need the radius cap

def l1_point(Aout, Cout, k, sol=None):
    """rigorous lower bound E <= lambda_1(dart(Aout, Cout)) (gcert at a point: CR + Temple + Liu), with sigma2 and diagnostics"""
    from gsigma import sigma_box
    gc = _gc('bd', k)
    pc = np.array([Aout[0], Aout[1], Cout[0], Cout[1]], float); rc = np.full(4, RPT)   # a tiny box around the point (the   # RIGOR: exact: float64 arrays of Python floats (no rounding)
    if sol is None: sol = gc.solve(pc, rc)                                               # Taylor-model code needs r > 0)
    s, info = sigma_box(gc, sol, pc, rc)
    if s is None: return None, dict(why='sigma: ' + str(info.get('why')), lam1h=float(sol['lamh'][0]))   # RIGOR: info: logged
    res = gc.certify(pc, rc, s, sol=sol)
    if 'L1' not in res: return None, dict(why='certify: ' + str(res.get('why')), lam1h=float(sol['lamh'][0]))   # RIGOR: info: logged
    E = res['L1'][0]                                         # tm_rng: the exact lower end of the L1 range rounded DOWN
    return E, dict(sigma=s, frac=info.get('frac'), lam1h=float(sol['lamh'][0]), h=res.get('h'), temple=res.get('temple'))   # RIGOR: info: logged

LIU = 0.1893
def h_est(Ao, Co, k):   # RIGOR-FUNCTION: heuristic: pre-check estimate
    """rough estimate of the CR mesh size h of the BD chart at level k (pre-check only; calibrated on certificates)"""
    D = (1.0, 0.0)
    diam = max(1.0, math.hypot(*Ao), math.hypot(*Co), math.hypot(Ao[0] - D[0], Ao[1]), math.hypot(Co[0] - D[0], Co[1]))
    return 0.09*diam*2.0**(4 - k)

def inc_float(p0, r, k=4, kp=4):   # RIGOR-FUNCTION: heuristic: float pre-check
    """float pre-check at CR level k (P2 level kp): (darts, sol at Q+, (lam1h, lam2p2), Liu-corrected ratio) or (None, ...)"""
    from p2up import solve_p2
    dd = inc_darts(p0, r)
    if dd is None: return None, None, None, float('inf')
    Ao, Co, Ai, Ci = dd
    gc = _gc('bd', k)
    sol = gc.solve(np.array([Ao[0], Ao[1], Co[0], Co[1]]), np.full(4, RPT))
    l1 = float(sol['lamh'][0])
    l2 = float(solve_p2(_gc('bd', kp), np.array([Ai[0], Ai[1], Ci[0], Ci[1]]), np.zeros(4), order=0)['lamp2'][1])
    liu = l1/(1 + (LIU*h_est(Ao, Co, k))**2*l1)
    return dd, sol, (l1, l2), l2/liu

def inc_cert(p0, r, k0=4, kmax=6, pre=2.45):   # RIGOR-FUNCTION: heuristic: search driver; every decision is made by verify_inclusion, l1_point and p2_cert
    """INC certificate for the PP box (p0, r), with an adaptive CR mesh level for the outer dart.  Start at level k0; while the
    (Liu-corrected) float ratio is in [pre, esc) try the next level (esc = 2.9, or 3.5 for small boxes, max r <= 0.008), and
    continue only while each step lowers the ratio by more than 15% (strongly mesh-limited; such a level is inherited by the
    children as 'ik').  All checks rigorous except the float pre-checks.  Returns dict(ok, ...); on failure 'ratio' (best
    float ratio) and 'ik' guide the bisection."""
    esc = 3.5 if float(np.max(r)) <= 0.008 else 2.9
    best = float('inf'); ik = k0; prev = None
    for k in range(k0, kmax + 1):
        dd, sol, lams, ratio = inc_float(p0, r, k)
        if dd is None: return dict(ok=False, why='inc: no outer/inner darts', ratio=float('inf'), ik=k0)
        if prev is not None:
            if ratio < 0.85*prev: ik = k                        # strongly mesh-limited: children start at this level
            elif not ratio < pre: break                         # a finer mesh does not pay for this box: split instead
        best = min(best, ratio)
        if ratio < pre:
            Ao, Co, Ai, Ci = dd
            if not verify_inclusion(p0, r, Ao, Co, Ai, Ci): return dict(ok=False, why='inc: inclusion check', ratio=ratio, ik=ik)
            E, info = l1_point(Ao, Co, k, sol=sol)
            if E is not None:
                for kp in (4, 5, 6):
                    res = p2_cert(Ai, Ci, E, kp)
                    if res.get('ok'):
                        return dict(ok=True, method='INC', k=k, kp=kp, Aout=list(Ao), Cout=list(Co), Ain=list(Ai), Cin=list(Ci),
                                    E=E, sigma=info['sigma'], lam1h=info['lam1h'], h=info['h'], temple=info['temple'],
                                    schur=res['schur'], lam2p2=res['lam2p2'], ratio=ratio, ik=ik)
                    if res.get('lam2p2', float('inf')) > 2.5*E/0.85**(6 - kp): break
        if not ratio < esc: break
        prev = ratio
    return dict(ok=False, why='inc: float ratio %.3f (ik=%d)' % (best, ik), ratio=best, ik=ik)
