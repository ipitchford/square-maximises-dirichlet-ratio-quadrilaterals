"""Multivariate Taylor models in p=(a,b,c,d) on the Euclidean ball |p|_2 <= rho0.

An MTM is  f(p) = sum_{|g|<=m} c_g p^g + R(p),  |R(p)| <= r |p|^{m+1},  with c_g arb balls and r an arb (its upper
endpoint is the bound).  Bounds use |p^g| <= sig(g) |p|^{|g|}, sig(g) = max_{|e|=1} |e^g| = prod (g_i/|g|)^{g_i/2},
which is submultiplicative: sig(g1+g2) <= sig(g1) sig(g2) (max of a product <= product of maxima).
"""
from flint import arb, ctx as flintctx
from localcert import sigma_arb
flintctx.prec = 256

_SIG = {}
def sig(g):
    v = _SIG.get(g)
    if v is None:
        v = sigma_arb(g); _SIG[g] = v
    return v

def deg(g): return g[0] + g[1] + g[2] + g[3]
ZERO4 = (0, 0, 0, 0)

class MTM:
    __slots__ = ('c', 'r', 'm', 'rho', '_prof', '_bydeg')
    def __init__(self, c, r, m, rho):
        self.c = {g: v for g, v in c.items() if not (v.is_exact() and v == 0)}
        self.r = arb(r); self.m = m; self.rho = arb(rho); self._prof = None; self._bydeg = None
    @staticmethod
    def const(x, m, rho): return MTM({ZERO4: arb(x)}, 0, m, rho)
    def bydeg(self):
        if self._bydeg is None:
            B = [[] for _ in range(self.m + 1)]
            for g, v in self.c.items(): B[deg(g)].append((g, v))
            self._bydeg = B
        return self._bydeg
    def profiles(self):
        """P[k] >= sup_{|e|=1} |sum_{|g|=k} c_g e^g|"""
        if self._prof is None:
            P = [arb(0) for _ in range(self.m + 1)]
            for g, v in self.c.items(): P[deg(g)] += abs(v) * sig(g)
            self._prof = P
        return self._prof
    def sup(self):
        P = self.profiles(); s = arb(0); rk = arb(1)
        for k in range(self.m + 1):
            s += P[k] * rk; rk *= self.rho
        return s + self.r * rk
    def __add__(self, o):
        if not isinstance(o, MTM): o = MTM.const(o, self.m, self.rho)
        c = dict(self.c)
        for g, v in o.c.items(): c[g] = c[g] + v if g in c else v
        return MTM(c, self.r + o.r, self.m, self.rho)
    __radd__ = __add__
    def __neg__(self): return MTM({g: -v for g, v in self.c.items()}, self.r, self.m, self.rho)
    def __sub__(self, o): return self + (-o)
    def __rsub__(self, o): return (-self) + o
    def scale(self, x):
        x = arb(x); return MTM({g: v * x for g, v in self.c.items()}, self.r * abs(x), self.m, self.rho)
    def __mul__(self, o):
        if not isinstance(o, MTM): return self.scale(o)
        m = self.m; rho = self.rho; A = self.bydeg(); B = o.bydeg(); out = {}
        for k1 in range(m + 1):
            for k2 in range(m - k1 + 1):
                for g1, v1 in A[k1]:
                    for g2, v2 in B[k2]:
                        g = (g1[0]+g2[0], g1[1]+g2[1], g1[2]+g2[2], g1[3]+g2[3]); t = v1 * v2
                        out[g] = out[g] + t if g in out else t
        PA = self.profiles(); PB = o.profiles(); drop = arb(0)
        for k1 in range(m + 1):
            for k2 in range(m + 1):
                if k1 + k2 > m: drop += PA[k1] * PB[k2] * rho ** (k1 + k2 - m - 1)
        r = drop + self.r * o.sup() + o.r * self.sup() + self.r * o.r * rho ** (m + 1)
        return MTM(out, r, m, rho)
    __rmul__ = __mul__
    def inv(self):
        """1/f for f with constant term c0 not containing 0: f = c0(1+g), 1/(1+g) = sum_{n<=m} (-g)^n + tail,
        |g| <= G|p|, |tail| <= (G|p|)^{m+1}/(1-G rho)."""
        m = self.m; rho = self.rho
        c0 = self.c.get(ZERO4)
        assert c0 is not None and not c0.contains(0), "inverse: constant term contains 0"
        i0 = 1 / c0
        g = MTM({k: v * i0 for k, v in self.c.items() if deg(k) > 0}, self.r * abs(i0), m, rho)
        P = g.profiles(); G = arb(0); rk = arb(1)
        for k in range(1, m + 1):
            G += P[k] * rk; rk *= rho
        G += g.r * rk                      # rk = rho^m
        q = G * rho
        assert q < 1, f"inverse: G*rho = {q} not < 1"
        total = MTM.const(1, m, rho); power = MTM.const(1, m, rho); mg = -g
        for n in range(1, m + 1):
            power = power * mg; total = total + power
        total = MTM(total.c, total.r + G ** (m + 1) / (1 - q), m, rho)
        return total.scale(i0)
    def zero_degree(self, k, tol=None):
        """Remove degree-k coefficients that are known to vanish exactly (by a proved identity); each must contain 0."""
        for g in [g for g in self.c if deg(g) == k]:
            assert self.c[g].contains(0), f"zero_degree({k}): coefficient {g} = {self.c[g]} does not contain 0"
            del self.c[g]
        self._prof = None; self._bydeg = None
        return self
    def degree_coeffs(self, k): return {g: v for g, v in self.c.items() if deg(g) == k}
