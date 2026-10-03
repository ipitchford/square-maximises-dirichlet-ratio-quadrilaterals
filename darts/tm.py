"""Order-2 Taylor models in NV normalised variables t in [-1,1]^NV, coefficients = arb balls, remainder = arb >= 0.
Semantics: f(t) in sum_{|alpha|<=2} c_alpha t^alpha + [-R, R] for every t in the box.
Products: all terms of degree 3,4 are bounded with |t^alpha| <= 1 and folded into R (the fold is explicit).
Reciprocal: f = f0 + f',  |f'| <= B < |f0|:  1/f = (1/f0)(1 - u + u^2) + tail, |tail| <= q^3/((1-q)|f0|), q = B/|f0|."""
from flint import arb, ctx
from nanguard import afin, amin, amax
ctx.prec = 106
import os, itertools
NV = 4
ORDER = int(os.environ.get('TM_ORDER', '3'))
MON = [m for d in range(ORDER + 1) for m in itertools.combinations_with_replacement(range(NV), d)]
IDX = {m: k for k, m in enumerate(MON)}
NC = len(MON)
DEG = [len(m) for m in MON]
ZERO = arb(0)

def _mulidx():
    out = []
    for a, ma in enumerate(MON):
        for b, mb in enumerate(MON):
            m = tuple(sorted(ma + mb)); out.append((a, b, IDX.get(m, -1)))
    return out
MULIDX = _mulidx()

class TM:
    __slots__ = ('c', 'R')
    def __init__(self, c=None, R=None):
        self.c = c if c is not None else [arb(0)]*NC
        self.R = R if R is not None else arb(0)
    @staticmethod
    def const(x):
        c = [arb(0)]*NC; c[0] = arb(x); return TM(c)
    @staticmethod
    def var(i, mid, rad):
        """p_i = mid + rad * t_i"""
        c = [arb(0)]*NC; c[0] = arb(mid); c[1+i] = arb(rad); return TM(c)
    def copy(self): return TM(list(self.c), self.R)
    def __add__(self, o):
        if not isinstance(o, TM): o = TM.const(o)
        return TM([a+b for a, b in zip(self.c, o.c)], self.R + o.R)
    __radd__ = __add__
    def __neg__(self): return TM([-a for a in self.c], self.R)
    def __sub__(self, o):
        if not isinstance(o, TM): o = TM.const(o)
        return self + (-o)
    def __rsub__(self, o): return (-self) + o
    def absbound(self):
        """upper bound for |f| over the box (crude: sum of |coefficients| + R)"""
        s = self.R
        for a in self.c: s = s + abs(a)
        return s
    def __mul__(self, o):
        if not isinstance(o, TM):
            o = arb(o); return TM([a*o for a in self.c], self.R*abs(o))
        c = [arb(0)]*NC; drop = arb(0)
        for a, b, k in MULIDX:
            if self.c[a] == 0 or o.c[b] == 0: continue
            p = self.c[a]*o.c[b]
            if k >= 0: c[k] = c[k] + p
            else: drop = drop + abs(p)
        R = drop + self.R*o.absbound() + o.R*self.absbound() + self.R*o.R
        return TM(c, R)
    __rmul__ = __mul__
    def nonconst_bound(self):
        s = self.R
        for a in self.c[1:]: s = s + abs(a)
        return s
    def recip(self):
        f0 = self.c[0]
        B = self.nonconst_bound()
        q = B/abs(f0)
        if not (q < 1): raise ValueError('reciprocal: TM not bounded away from 0')
        u = self.copy(); u.c = [arb(0)] + [x/f0 for x in self.c[1:]]; u.R = self.R/abs(f0)
        u2 = u*u
        res = (TM.const(1) - u + u2)*(1/f0)
        # tail bound q^3/((1-q)|f0|) is increasing in q: use the exact upper end of q (python-flint's ball ** 3 is NaN for a
        # ball containing 0, which happens when the non-constant part of the model is itself a ball around 0)
        qu = arb(q.upper())
        res.R = res.R + qu*qu*qu/((1 - qu)*abs(f0))
        return res
    def __truediv__(self, o):
        if not isinstance(o, TM): return self*(1/arb(o))
        return self*o.recip()
    def range(self):
        """arb enclosure of the range over the box: exact 1-D quadratic ranges for t_i^2, t_i terms; |c| for all other monomials"""
        for c_ in self.c: afin(c_)
        afin(self.R)
        lo = self.c[0] - self.R; hi = self.c[0] + self.R
        lo = arb(lo.lower()); hi = arb(hi.upper())
        used = {0}
        for i in range(NV):
            ia = IDX[(i, i)] if ORDER >= 2 else None; ib = IDX[(i,)]
            a = self.c[ia] if ia is not None else arb(0); b = self.c[ib]
            used.add(ib)
            if ia is not None: used.add(ia)
            # range of a t^2 + b t over t in [-1, 1]: the end values a -+ b, the value 0 (t = 0) and, if the vertex -b/(2a) may lie in
            # [-1, 1], the vertex value -b^2/(4a).  When the ball a contains 0 the vertex value is not formed (it would be a NaN
            # ball): then a t^2 lies in [min(0, a_lo), max(0, a_hi)] and |b t| <= |b|, which encloses every value.
            cands = [afin(a - b), afin(a + b)]
            if a.lower() <= 0 and a.upper() >= 0:
                if not (a == 0):
                    bm = abs(b).upper()
                    cands += [afin(arb(a.lower()) - bm), afin(arb(a.upper()) + bm)]
            elif (abs(b) - 2*abs(a)).lower() < 0:
                cands.append(afin(-(b*b)/(4*a)))
            mn = amin(arb(0), *cands).lower(); mx = amax(arb(0), *cands).upper()      # exact points; t = 0 gives the value 0
            lo = lo + mn; hi = hi + mx
        for k, c in enumerate(self.c):
            if k in used: continue
            a = abs(c); lo = lo - a; hi = hi + a
        return arb(lo.lower()), arb(hi.upper())
    def upper(self): return self.range()[1]
    def lower(self): return self.range()[0]
    def __repr__(self):
        return f"TM(c0={self.c[0]}, lin={[str(x) for x in self.c[1:1+NV]]}, R={self.R})"
