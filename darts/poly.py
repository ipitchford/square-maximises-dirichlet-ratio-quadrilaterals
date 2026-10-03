"""Exact multivariate polynomials in t = (t0..t3) with arb coefficients (no truncation); used for residual coefficient
functions.  P = dict{exponent tuple: arb}.  Conversion from TM: polynomial part (deg <= 2) + remainder R."""
from flint import arb
from tm import TM, MON, NV

class Poly:
    __slots__ = ('d',)
    def __init__(self, d=None): self.d = d if d is not None else {}
    @staticmethod
    def const(c): return Poly({(0,)*NV: arb(c)})
    @staticmethod
    def var(i, mid, rad):
        e = [0]*NV; e[i] = 1
        return Poly({(0,)*NV: arb(mid), tuple(e): arb(rad)})
    @staticmethod
    def mono(alpha):
        """alpha indexes tm.MON: 0 -> 1, 1..4 -> t_i, 5.. -> t_i t_j"""
        e = [0]*NV
        for i in MON[alpha]: e[i] += 1
        return Poly({tuple(e): arb(1)})
    def __add__(self, o):
        if not isinstance(o, Poly): o = Poly.const(o)
        d = dict(self.d)
        for k, v in o.d.items(): d[k] = d[k] + v if k in d else v
        return Poly(d)
    __radd__ = __add__
    def __neg__(self): return Poly({k: -v for k, v in self.d.items()})
    def __sub__(self, o): return self + (-o if isinstance(o, Poly) else Poly.const(-arb(o)))
    def __rsub__(self, o): return (-self) + o
    def __mul__(self, o):
        if not isinstance(o, Poly):
            o = arb(o); return Poly({k: v*o for k, v in self.d.items()})
        d = {}
        for k1, v1 in self.d.items():
            for k2, v2 in o.d.items():
                k = tuple(a+b for a, b in zip(k1, k2)); p = v1*v2
                d[k] = d[k] + p if k in d else p
        return Poly(d)
    __rmul__ = __mul__
    def absbound(self):
        s = arb(0)
        for v in self.d.values(): s = s + abs(v)
        return s
    @staticmethod
    def from_tm(T):
        d = {}
        for c, m in zip(T.c, MON):
            e = [0]*NV
            for i in m: e[i] += 1
            d[tuple(e)] = c
        return Poly(d), T.R
