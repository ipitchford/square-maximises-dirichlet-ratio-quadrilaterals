"""Concrete charts and problems.
 BD   : pieces A (B,D,A), C (B,D,C), both graded toward D; params (ax, ay, cx, cy).
 SPIKE: Lemma 9.1 -- pieces H (B,D,E), E=(1,t), and C; params (t, -, cx, cy); CR: Neumann on DE; P2: chart C alone (CBD).
        (lambda_2(Q) <= lambda_2(P) <= lambda_2(CBD); using CBD alone keeps the P2 chart t-independent)
"""
from fractions import Fraction as Fr
from poly import Poly
from chart import Chart
from gcert import Problem

UNIT = ((Fr(0), Fr(0)), (Fr(1), Fr(0)), (Fr(0), Fr(1)))
def P(v): return v if isinstance(v, Poly) else Poly.const(v)

# ---------------------------------------------------------------------------------------------------------------
# Rigorous tail constants (directed: each returned float is <= the true constant).
def spike_tail(p0, r):
    """pi^2 / t_hi^2 with t_hi = p0_0 + r_0 the upper end of the represented t-range (exact in ball arithmetic), rounded down"""
    from flint import arb
    from nanguard import fdown
    th = arb(p0[0]) + arb(r[0])
    if not (th.lower() > 0): raise ValueError('spike tail: t_hi <= 0')
    return fdown(arb.pi()*arb.pi()/(th*th))
def cut_tail(c, h, k0):
    """lower bound (float, rounded down) of pi^2 / max |Q1Q2(u)|^2 over u in [c - h, c + h] (exact ends), where
    |Q1Q2(u)|^2 = (1 - k0 u^2)^2 + (k0 u)^2 = g(u^2),  g(s) = (1 - k0 s)^2 + k0^2 s  (k0 a positive rational).
    g is a convex quadratic in s (g'' = 2 k0^2 > 0), so its maximum over the s-range [s_lo, s_hi] of u^2 is max(g(s_lo), g(s_hi));
    s_hi = max(lo^2, hi^2), s_lo = 0 if 0 in [lo, hi] else min(lo^2, hi^2).  Everything in ball arithmetic, then rounded up / down."""
    from flint import arb
    from nanguard import fdown, amax, amin, afin
    k0 = Fr(k0)
    if not k0 > 0: raise ValueError('cut tail: k0 <= 0')
    kk = arb(k0.numerator)/arb(k0.denominator)
    lo = arb(c) - arb(h); hi = arb(c) + arb(h)
    s_hi = amax(lo*lo, hi*hi)
    s_lo = arb(0) if (lo.lower() <= 0 and hi.upper() >= 0) else amin(lo*lo, hi*hi)
    def g(sv):
        w = 1 - kk*sv
        return afin(w*w + kk*kk*sv)
    m = amax(g(s_lo), g(s_hi))
    if not (m.lower() > 0): raise ValueError('cut tail: |Q1Q2| may vanish')
    return fdown(arb.pi()*arb.pi()/m)

def bd_chart(k, fresh=False):
    return Chart([dict(name='A', ref=UNIT, labels=('B', 'D', 'A'), grade=1),
                  dict(name='C', ref=UNIT, labels=('B', 'D', 'C'), grade=1)], k, fresh=fresh)
def bd_vfun(Z):
    return {'B': (P(0), P(0)), 'D': (P(1), P(0)), 'A': (P(Z[0]), P(Z[1])), 'C': (P(Z[2]), P(Z[3]))}
def bd_problem(k):
    ch = bd_chart(k); return Problem(ch, ch, bd_vfun, name=f'BD_k{k}')

def spike_cr_chart(k, fresh=False):
    bc = {frozenset(('D', 'E')): {'cr': 'N', 'p2': 'D'}}
    return Chart([dict(name='H', ref=UNIT, labels=('B', 'D', 'E'), grade=1),
                  dict(name='C', ref=UNIT, labels=('B', 'D', 'C'), grade=1)], k, bc=bc, fresh=fresh)
def c_only_chart(k, fresh=False):
    return Chart([dict(name='C', ref=UNIT, labels=('B', 'D', 'C'), grade=1)], k, fresh=fresh)
def spike_vfun(Z):
    return {'B': (P(0), P(0)), 'D': (P(1), P(0)), 'E': (P(1), P(Z[0])), 'C': (P(Z[2]), P(Z[3]))}
def spike_problem(k, kp2=None):
    def tails(p0, r): return [spike_tail(p0, r)]
    return Problem(spike_cr_chart(k), c_only_chart(kp2 or k), spike_vfun, tails_fn=tails, name=f'SPIKE_k{k}')

# ---------------------------------------------------------------------------------------------------------------
# SPIKE2 (Lemma 9.4 in LEMMAS.md): long spikes a_x >= Xc := 1 + s_c * t ... we use a FIXED stub length factor:
#   E = (1, t), D' = (1 + w, t w), E' = (1 + w, t (1 + w)), w = W*t  (stub of length w and vertical width t, walls parallel to BA)
#   CR chart: pieces C (B,D,C), H (B,D,E), S1 (D, D', E'), S2 (D, E', E); Neumann on D'E'.  (core of the limit spike)
#   P2 chart: BD chart with A = X0 (1, t) (the shortest admissible spike; lambda_2 is non-increasing in |A| on the ray).
#   params (t, -, cx, cy)
SPIKE2_W = Fr(3)        # stub length / t
def spike2_cr_chart(k, fresh=False):
    bc = {frozenset(("D'", "E'")): {'cr': 'N', 'p2': 'D'}}
    return Chart([dict(name='C', ref=UNIT, labels=('B', 'D', 'C'), grade=1),
                  dict(name='H', ref=UNIT, labels=('B', 'D', 'E'), grade=1),
                  dict(name='S1', ref=UNIT, labels=('D', "D'", "E'"), grade=0),
                  dict(name='S2', ref=UNIT, labels=('D', "E'", 'E'), grade=0)], k, bc=bc, fresh=fresh)
def spike2_vfun(Z):
    t = P(Z[0]); w = t*float(SPIKE2_W)   # RIGOR: exact: SPIKE2_W = 3
    return {'B': (P(0), P(0)), 'D': (P(1), P(0)), 'E': (P(1), t), "D'": (w + 1, t*w), "E'": (w + 1, t*(w + 1)),
            'C': (P(Z[2]), P(Z[3]))}
def spike2_p2_vfun_factory(X0):
    def vf(Z):
        t = P(Z[0])
        return {'B': (P(0), P(0)), 'D': (P(1), P(0)), 'A': (P(X0), t*X0), 'C': (P(Z[2]), P(Z[3]))}
    return vf
class _Combo:
    """vertex function used by both charts: union of the label maps"""
    def __init__(self, f1, f2): self.f1, self.f2 = f1, f2
    def __call__(self, Z):
        d = dict(self.f1(Z)); d.update(self.f2(Z)); return d
def spike2_problem(k, X0, kp2=None):
    def tails(p0, r): return [spike_tail(p0, r)]
    vf = _Combo(spike2_vfun, spike2_p2_vfun_factory(X0))
    return Problem(spike2_cr_chart(k), bd_chart(kp2 or k), vf, tails_fn=tails, name=f'SPIKE2_k{k}_X0{X0}')

# ---------------------------------------------------------------------------------------------------------------
# FAR (Lemma 9.3 with a near-perpendicular cut): arm A with direction (u, 1), u = cot b_A, A = rho (u, 1), rho >= RA.
#   cut endpoints Q1 = Lc (u, 1) (on BA), Q2 = D + (Lc - k0 u)(u, 1) (on the parallel through D), k0 = 1/(1 + u0^2) dyadic-rounded
#   from the box centre, so Q1Q2 is perpendicular to the arm at u = u0.  CR chart: pieces C, S1 = (B, D, Q2), S2 = (B, Q2, Q1),
#   Neumann on Q1Q2.  Tail: lambda1 >= pi^2/|Q1Q2|^2 (chords of the cut-off arm parallel to the cut, Lemma 6.1).
#   P2 chart: BD chart with A = RA (u, 1)  (lambda_2 non-increasing along the ray).   params (u, -, cx, cy)
FAR_LC = Fr(3, 2)
def _k0_from(Z0):   # RIGOR-FUNCTION: heuristic: free dyadic choice of the cut slope from the box centre
    import math
    u0 = float(Z0.d.get((0, 0, 0, 0)).mid()) if isinstance(Z0, Poly) else float(Z0)
    k = 1.0/(1.0 + u0*u0)
    return Fr(round(k*4096), 4096)
def far_cr_chart(k, fresh=False):
    bc = {frozenset(('Q1', 'Q2')): {'cr': 'N', 'p2': 'D'}}
    return Chart([dict(name='C', ref=UNIT, labels=('B', 'D', 'C'), grade=1),
                  dict(name='S1', ref=UNIT, labels=('B', 'D', 'Q2'), grade=1),
                  dict(name='S2', ref=UNIT, labels=('B', 'Q2', 'Q1'), grade=None)], k, bc=bc, fresh=fresh)
def far_vfun_factory(RA):
    def vf(Z):
        u = P(Z[0]); k0 = _k0_from(Z[0]); L = float(FAR_LC)   # RIGOR: exact: FAR_LC = 3/2
        L2 = u*(-float(k0)) + L   # RIGOR: exact: k0 is a dyadic m/4096
        return {'B': (P(0), P(0)), 'D': (P(1), P(0)), 'Q1': (u*L, P(L)), 'Q2': (u*L2 + 1, L2),
                'A': (u*RA, P(RA)), 'C': (P(Z[2]), P(Z[3]))}
    return vf
def far_problem(k, RA=3.0, kp2=None):
    def tails(p0, r): return [cut_tail(p0[0], r[0], _k0_from(p0[0]))]
    return Problem(far_cr_chart(k), bd_chart(kp2 or k), far_vfun_factory(RA), tails_fn=tails, name=f'FAR_k{k}_RA{RA}')

# DEEP: both arms cut.  params (uA, -, uC, -): A = RA (uA, 1), C = RC (uC, -1)
def deep_cr_chart(k, fresh=False):
    bc = {frozenset(('Q1', 'Q2')): {'cr': 'N', 'p2': 'D'}, frozenset(('R1', 'R2')): {'cr': 'N', 'p2': 'D'}}
    return Chart([dict(name='S1', ref=UNIT, labels=('B', 'D', 'Q2'), grade=1),
                  dict(name='S2', ref=UNIT, labels=('B', 'Q2', 'Q1'), grade=None),
                  dict(name='T1', ref=UNIT, labels=('B', 'D', 'R2'), grade=1),
                  dict(name='T2', ref=UNIT, labels=('B', 'R2', 'R1'), grade=None)], k, bc=bc, fresh=fresh)
def deep_vfun_factory(RA, RC):
    def vf(Z):
        uA = P(Z[0]); uC = P(Z[2]); kA = float(_k0_from(Z[0])); kC = float(_k0_from(Z[2])); L = float(FAR_LC)   # RIGOR: exact: dyadic k0 and FAR_LC = 3/2
        LA = uA*(-kA) + L; LC = uC*(-kC) + L
        return {'B': (P(0), P(0)), 'D': (P(1), P(0)), 'Q1': (uA*L, P(L)), 'Q2': (uA*LA + 1, LA),
                'R1': (uC*L, P(-L)), 'R2': (uC*LC + 1, -LC), 'A': (uA*RA, P(RA)), 'C': (uC*RC, P(-RC))}
    return vf
def deep_problem(k, RA=3.0, RC=3.0, kp2=None):
    def tails(p0, r): return [cut_tail(p0[i], r[i], _k0_from(p0[i])) for i in (0, 2)]
    return Problem(deep_cr_chart(k), bd_chart(kp2 or k), deep_vfun_factory(RA, RC), tails_fn=tails, name=f'DEEP_k{k}')

# ---------------------------------------------------------------------------------------------------------------
# FAR2: as FAR (CR core with the near-perpendicular Neumann cut at L = FAR_LC), but the P2 upper bound uses the dart with
#   A = R (u, 1) and its arm truncated by a Dirichlet cut parallel to BD at height L2 (pieces C, (B, D, P2), (B, P2, P1)),
#   P1 = L2 (u, 1) on BA, P2 = D + (L2/R)(A - D) = (1 + L2 u - L2/R, L2) on DA.  For rho >= R the truncated arm of Q(rho) contains
#   that of Q(R) (A(R) on the segment B A(rho) and the cut is the same line), so lambda_2(Q(rho)) <= lambda_2(truncated Q(R)).
FAR2_L2 = 1.5
def far2_p2_chart(k, fresh=False):
    return Chart([dict(name='C', ref=UNIT, labels=('B', 'D', 'C'), grade=1),
                  dict(name='U1', ref=UNIT, labels=('B', 'D', 'P2'), grade=1),
                  dict(name='U2', ref=UNIT, labels=('B', 'P2', 'P1'), grade=None)], k, fresh=fresh)
def far2_vfun_factory(RA, L2=None):
    L2 = FAR2_L2 if L2 is None else L2
    base = far_vfun_factory(RA)
    def vf(Z):
        d = base(Z); u = P(Z[0])
        d['P1'] = (u*L2, P(L2)); d['P2'] = (u*L2 + (1 - L2/RA), P(L2))
        return d
    return vf
def far2_problem(k, RA=2.0, kp2=None):
    import numpy as np
    pb_old = far_problem(k, RA)
    return Problem(far_cr_chart(k), far2_p2_chart(kp2 or k), far2_vfun_factory(RA), tails_fn=pb_old.tails_fn, name=f'FAR2_k{k}_RA{RA}')

# ---------------------------------------------------------------------------------------------------------------
# FARL: the FAR chart with a cut level L chosen per box (cover5 FEM fallback of the far-arm families, LEMMAS 9.3 / 11.3).
#   A = a_y (u, 1) with a_y >= RA (any length); cut Q1 = L (u, 1) on BA, Q2 = D + (L - k0 u)(u, 1) on the parallel through D,
#   k0 = 1/(1 + u0^2) (dyadic) from the box centre (the cut is perpendicular to the arm at u = u0).  CR core: C, (B,D,Q2), (B,Q2,Q1)
#   with Neumann on Q1Q2; tail pi^2/|Q1Q2|^2 (Lemma 6.1); P2: the dart with A = RA (u, 1) (lambda_2 non-increasing in a_y).
#   VALIDITY (checked by the caller, fam5.farl_valid): on the whole box the C-triangle and D lie strictly on B's side of the cut
#   line, L - k0 u >= 0, and RA > L (the far part beyond the cut is exactly the arm tip, a triangle with sides on BA and DA).
_FARL_CACHE = {}
def farl_vfun_factory(RA, L, k0):
    def vf(Z):
        u = P(Z[0]); L2 = u*(-float(k0)) + float(L)   # RIGOR: exact: dyadic parameters, checked exactly representable in farl_problem
        return {'B': (P(0), P(0)), 'D': (P(1), P(0)), 'Q1': (u*float(L), P(float(L))), 'Q2': (u*L2 + 1, L2),   # RIGOR: exact: dyadic parameters, checked exactly representable
                'A': (u*float(RA), P(float(RA))), 'C': (P(Z[2]), P(Z[3]))}   # RIGOR: exact: dyadic parameters, checked exactly representable
    return vf
from nanguard import fin
def farl_problem(k, RA, L, k0):
    """RA, L, k0: Fractions (dyadic); cached per (k, RA, L, k0)"""
    key = (k, RA, L, k0)
    if key in _FARL_CACHE: return _FARL_CACHE[key]
    for x in (RA, L, k0):
        if not (isinstance(x, Fr) and Fr(float(x)) == x): raise ValueError('FARL parameters must be exactly representable rationals')   # RIGOR: exact: the check itself
    if not (k0 > 0 and L > 0 and RA > L): raise ValueError('FARL parameters violate Lemma 9.6 (k0 > 0, R > L > 0)')
    def tails(p0, r): return [cut_tail(p0[0], r[0], k0)]
    base = _FARL_CACHE.get(('base', k))
    if base is None:
        base = Problem(far_cr_chart(k), bd_chart(k), farl_vfun_factory(RA, L, k0), tails_fn=tails, name='FARL_base')
        _FARL_CACHE[('base', k)] = base
    import copy
    pb = copy.copy(base)                    # shares the (parameter-independent) reference matrices of the charts
    pb.farl = (RA, L, k0)                   # Lemma 9.6 parameters (exact rationals), for the containment tests
    pb.vfun = farl_vfun_factory(RA, L, k0); pb.tails_fn = tails; pb.name = f'FARL_k{k}_RA{float(RA)}_L{float(L)}_k0{float(k0)}'   # RIGOR: info: name string
    _FARL_CACHE[key] = pb
    return pb
