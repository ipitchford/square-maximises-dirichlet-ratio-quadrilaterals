"""Exact Taylor-coefficient engine, centred reference square [-1,1]^2 (coordinates u,v), Euclidean ball |p|_2<=rho0.
Vertices (Endo-Osting displacements, shifted so that p=0 is the unit square [0,1]^2 in x-space):
  v1=(0,0)+d1, v2=(1,0)+d2, v3=(1,1)+d3, v4=(0,1)+d4.
Bilinear map Phi(u,v) = sum N_i v_i with N1=(1-u)(1-v)/4, N2=(1+u)(1-v)/4, N3=(1+u)(1+v)/4, N4=(1-u)(1+v)/4.
Modes g = 2 sin(m pi (u+1)/2) sin(n pi (v+1)/2): M(0)=I, K(0)=diag(pi^2(m^2+n^2)).
Taylor model: exact polynomial in (u,v,a,b,c,d) with p-degree <= m, plus remainder |R| <= r |p|^{m+1} on [-1,1]^2 x ball.
Bounds used for polynomials: |u^a v^b| <= 1 ; |p^gamma| <= sigma(gamma) |p|^{|gamma|}, sigma(gamma)=prod (gamma_i/|gamma|)^{gamma_i/2}.
"""
from fractions import Fraction as Fr
from math import comb, isqrt
from flint import fmpq_mpoly_ctx, fmpq_mpoly, fmpq, arb, arb_mat, Ordering, ctx as flintctx
import time

flintctx.prec = 256
CTX = fmpq_mpoly_ctx.get(('u','v','a','b','c','d'), Ordering.lex)
U,V,A,B_,C_,D_ = CTX.gens()
qtr = fmpq(1,4); half = fmpq(1,2)

def pdeg(exp): return int(exp[2])+int(exp[3])+int(exp[4])+int(exp[5])
def gam(exp): return (int(exp[2]),int(exp[3]),int(exp[4]),int(exp[5]))

_sigma_cache = {}
def sigma(g):
    """rational upper bound of max_{|e|_2=1} |e^g| = prod (g_i/k)^{g_i/2}."""
    if g in _sigma_cache: return _sigma_cache[g]
    k = sum(g)
    if k == 0: val = Fr(1)
    else:
        # exact: sqrt( prod (g_i/k)^{g_i} ) ; compute rational upper bound of sqrt
        prod = Fr(1)
        for gi in g:
            if gi > 0: prod *= Fr(gi, k)**gi
        scale = 10**24
        q = isqrt(prod.numerator*scale*scale // prod.denominator) + 1
        val = Fr(q, scale)
    _sigma_cache[g] = val; return val

def trunc(poly, m):
    keep = {}; drop = {}
    for exp, co in poly.terms():
        (keep if pdeg(exp) <= m else drop)[exp] = co
    return CTX.from_dict(keep), CTX.from_dict(drop)

def bound_poly(poly, rho0, shift=0):
    """sum |c| sigma(gamma) rho0^{|gamma|-shift}  >= sup |poly| / |p|^shift  (requires |gamma|>=shift on all terms)."""
    s = Fr(0); r0 = Fr(rho0)
    for exp, co in poly.terms():
        g = gam(exp); k = sum(g) - shift
        assert k >= 0
        s += abs(Fr(int(co.p), int(co.q))) * sigma(g) * r0**k
    return s

class TM:
    __slots__ = ('poly','r','m','rho0','_sup')
    def __init__(self, poly, r, m, rho0):
        self.poly = poly; self.r = Fr(r); self.m = m; self.rho0 = Fr(rho0); self._sup = None
    def sup(self):
        if self._sup is None: self._sup = bound_poly(self.poly, self.rho0) + self.r*self.rho0**(self.m+1)
        return self._sup
    @staticmethod
    def const(c, m, rho0): return TM(CTX.from_dict({(0,)*6: fmpq(Fr(c).numerator, Fr(c).denominator)}), 0, m, rho0)
    @staticmethod
    def exact(poly, m, rho0):
        keep, drop = trunc(poly, m)
        return TM(keep, bound_poly(drop, rho0, shift=m+1), m, rho0)
    def __add__(self, o):
        if isinstance(o, TM): return TM(self.poly+o.poly, self.r+o.r, self.m, self.rho0)
        return self + TM.const(o, self.m, self.rho0)
    __radd__ = __add__
    def __neg__(self): return TM(-self.poly, self.r, self.m, self.rho0)
    def __sub__(self, o): return self + (-o)
    def scale(self, c):
        c = Fr(c); return TM(self.poly*fmpq(c.numerator, c.denominator), abs(c)*self.r, self.m, self.rho0)
    def __mul__(self, o):
        if not isinstance(o, TM): return self.scale(o)
        assert o.m == self.m
        prod = self.poly * o.poly
        keep, drop = trunc(prod, self.m)
        r = bound_poly(drop, self.rho0, shift=self.m+1)
        r += self.r * o.sup() + o.r * self.sup() + self.r*o.r*self.rho0**(self.m+1)
        return TM(keep, r, self.m, self.rho0)
    __rmul__ = __mul__

def build_map(m, rho0):
    a,b,c,d = A,B_,C_,D_
    d1 = (half*(a+d-b), half*(-a+d+c)); d2 = (half*(-a+d+b), half*(-a-d-c))
    d3 = (half*(-a-d-b), half*(a-d+c)); d4 = (half*(a-d+b), half*(a+d-c))
    vx = [d1[0], 1+d2[0], 1+d3[0], d4[0]]; vy = [d1[1], d2[1], 1+d3[1], 1+d4[1]]
    N = [(1-U)*(1-V)*qtr, (1+U)*(1-V)*qtr, (1+U)*(1+V)*qtr, (1-U)*(1+V)*qtr]
    X = sum(N[i]*vx[i] for i in range(4)); Y = sum(N[i]*vy[i] for i in range(4))
    Xu, Xv, Yu, Yv = X.derivative(0), X.derivative(1), Y.derivative(0), Y.derivative(1)
    J = Xu*Yv - Xv*Yu
    adj = [[Yv, -Xv], [-Yu, Xu]]
    B = [[adj[i][0]*adj[j][0] + adj[i][1]*adj[j][1] for j in range(2)] for i in range(2)]
    dB = [[[B[i][j].derivative(al) for j in range(2)] for i in range(2)] for al in range(2)]
    dJ = [J.derivative(0), J.derivative(1)]
    E = lambda pol: TM.exact(pol, m, rho0)
    return dict(J=E(J), Jpoly=J, B=[[E(B[i][j]) for j in range(2)] for i in range(2)],
                dB=[[[E(dB[al][i][j]) for j in range(2)] for i in range(2)] for al in range(2)], dJ=[E(dJ[0]), E(dJ[1])], vx=vx, vy=vy)

def j1_bound(Jpoly, rho0):
    """J = 1/4 (1 + j1); return exact rational ghat with sup|j1| <= ghat |p| on domain, via linear part (CS at corners) + quadratic part (Frobenius)."""
    j1 = Jpoly*4 - 1
    lin = {}; quad = {}
    for exp, co in j1.terms():
        g = gam(exp); k = sum(g); e = (int(exp[0]), int(exp[1]))
        (lin if k==1 else quad)[(e, g)] = Fr(int(co.p), int(co.q))
        assert k in (1,2)
    # linear: f_i(u,v) affine ; sup sqrt(sum f_i^2) at corners
    best = Fr(0)
    for (uu,vv) in [(-1,-1),(1,-1),(-1,1),(1,1)]:
        ss = Fr(0)
        for i in range(4):
            val = sum(cf*(uu**e[0])*(vv**e[1]) for (e,g),cf in lin.items() if g[i]==1)
            ss += val*val
        best = max(best, ss)
    scale=10**24; s1 = Fr(isqrt(best.numerator*scale*scale//best.denominator)+1, scale)
    # quadratic: Frobenius of symmetric matrix with entries F_ii = f_ii, F_ij = f_ij/2 ; sup of each |f| at corners (affine in u,v)
    tot = Fr(0)
    for i in range(4):
        for j in range(i,4):
            sup = Fr(0)
            for (uu,vv) in [(-1,-1),(1,-1),(-1,1),(1,1)]:
                val = Fr(0)
                for (e,g),cf in quad.items():
                    if (i==j and g[i]==2) or (i!=j and g[i]==1 and g[j]==1):
                        val += cf*(uu**e[0])*(vv**e[1])
                sup = max(sup, abs(val))
            tot += sup*sup if i==j else sup*sup/2
    s2 = Fr(isqrt(tot.numerator*scale*scale//tot.denominator)+1, scale)
    return s1 + s2*Fr(rho0), s1, s2

def jinv_TM(mp, k, m, rho0):
    """TM for (4J)^{-k} = (1+j1)^{-k} = sum_n C(k+n-1,n)(-j1)^n."""
    r0 = Fr(rho0)
    ghat, s1, s2 = j1_bound(mp['Jpoly'], rho0)
    assert ghat*r0 < 1, float(ghat*r0)
    j1 = TM.exact(mp['Jpoly']*4 - 1, m, rho0); mj1 = -j1
    total = TM.const(1, m, rho0); power = TM.const(1, m, rho0)
    for n in range(1, m+1):
        power = power * mj1
        total = total + power.scale(comb(k+n-1, n))
    tail = Fr(0); n = m+1
    while True:
        term = comb(k+n-1, n) * ghat**n * r0**(n-m-1); tail += term
        if term < Fr(1, 10**40) or n > 4000: break
        n += 1
    q = ghat*r0*Fr(k+n, n+1); assert q < 1
    tail += term*q/(1-q)
    total.r += tail
    return total, ghat

def build_tensors(m, rho0, mLL=None):
    """TK[(al,be)] = B^{al be}/J ; J ; TLL[(i,j)] = Coef_i Coef_j J with Lg = sum_i Coef_i Der_i g."""
    if mLL is None: mLL = m
    t0 = time.time()
    mp = build_map(m, rho0)
    J1, ghat = jinv_TM(mp, 1, m, rho0)          # (4J)^{-1}
    TK = {}
    for al in range(2):
        for be in range(al,2):
            TK[(al,be)] = (J1 * mp['B'][al][be]).scale(4)      # B/J = 4 B (4J)^{-1}
    mpL = build_map(mLL, rho0) if mLL != m else mp
    J2, _ = jinv_TM(mpL, 2, mLL, rho0); J3, _ = jinv_TM(mpL, 3, mLL, rho0)   # (4J)^{-2}, (4J)^{-3}
    BL, dBL, dJL = mpL['B'], mpL['dB'], mpL['dJ']
    Dc = []
    for be in range(2):
        t = None
        for al in range(2):
            term = (J2 * dBL[al][al][be]).scale(16) - (J3 * (BL[al][be] * dJL[al])).scale(64)   # dB/J^2 - B dJ/J^3
            t = term if t is None else t + term
        Dc.append(t)
    Ecoef = {}
    for al in range(2):
        for be in range(al,2):
            Ecoef[(al,be)] = (J2 * BL[al][be]).scale(16)
    Coef = [Dc[0], Dc[1], Ecoef[(0,0)], Ecoef[(0,1)].scale(2), Ecoef[(1,1)]]
    TLL = {}
    for i in range(5):
        for j in range(i,5):
            TLL[(i,j)] = (Coef[i]*Coef[j]) * mpL['J']
    print(f"[tensors] m={m} mLL={mLL} rho0={float(rho0)} built in {time.time()-t0:.1f}s; ghat={float(ghat):.4f}; rem TK={[f'{float(t.r):.3g}' for t in TK.values()]} rem TLL max={max(float(t.r) for t in TLL.values()):.3g}; at rho0: TK {max(float(t.r) for t in TK.values())*float(rho0)**(m+1):.2e}, TLL {max(float(t.r) for t in TLL.values())*float(rho0)**(mLL+1):.2e}", flush=True)
    return dict(TK=TK, J=mpL['J'], Jm=mp['J'], TLL=TLL, ghat=ghat, map=mp, mapL=mpL, Coef=Coef)

# ---------------- trig tables on [-1,1]: Ic[a][k]=int x^a cos(k pi (x+1)/2), Is[a][k]=int x^a sin(k pi (x+1)/2) -----
def trig_tables(amax, kmax):
    pi = arb.pi()
    Ic = [[arb(0)]*(kmax+1) for _ in range(amax+1)]; Is = [[arb(0)]*(kmax+1) for _ in range(amax+1)]
    for a in range(amax+1):
        Ic[a][0] = (arb(1) - arb(-1)**(a+1))/(a+1); Is[a][0] = arb(0)
    for k in range(1, kmax+1):
        w = k*pi/2
        # int_{-1}^1 x^a cos(w x + w) dx : boundary sin(w x + w): at 1 -> sin(k pi)=0 ; at -1 -> sin 0 = 0
        # int x^a sin(w x + w) dx = [-x^a cos(wx+w)/w] + (a/w) int x^{a-1} cos(wx+w) ; boundary: -( cos(k pi) - (-1)^a cos 0 )/w
        for a in range(amax+1):
            if a == 0:
                Ic[0][k] = arb(0)
                Is[0][k] = -(arb(-1)**k - 1)/w
            else:
                Ic[a][k] = -(a/w)*Is[a-1][k]
                Is[a][k] = -(arb(-1)**k - arb(-1)**a)/w + (a/w)*Ic[a-1][k]
    return Ic, Is

def factor_int(a, f1, f2, Ic, Is):
    t1,m1 = f1; t2,m2 = f2; dm = abs(m1-m2); sm = m1+m2
    if t1=='s' and t2=='s': return (Ic[a][dm] - Ic[a][sm])/2
    if t1=='c' and t2=='c': return (Ic[a][dm] + Ic[a][sm])/2
    if t1=='s': ms, mc = m1, m2
    else: ms, mc = m2, m1
    sgn = 1 if ms >= mc else -1
    return (Is[a][sm] + sgn*Is[a][dm])/2

def deriv_desc(mode, der):
    """g=2 sin(m pi(u+1)/2) sin(n pi(v+1)/2); derivatives w.r.t. u,v on [-1,1]."""
    m,n = mode; pi = arb.pi(); wm = m*pi/2; wn = n*pi/2
    if der is None: return (arb(2), ('s',m), ('s',n))
    if der == 0: return (2*wm, ('c',m), ('s',n))
    if der == 1: return (2*wn, ('s',m), ('c',n))
    if der == 2: return (-2*wm*wm, ('s',m), ('s',n))
    if der == 3: return (2*wm*wn, ('c',m), ('c',n))
    if der == 4: return (-2*wn*wn, ('s',m), ('s',n))
    raise ValueError

class Contractor:
    """Contract a tensor TM (order tm.m, remainder r|p|^{m+1}) against trig products, returning a Taylor model of order
    mtrunc <= tm.m.  Terms of p-degree in (mtrunc, tm.m] are contracted EXACTLY and then folded into the remainder with
    |p^g| <= sigma(g) rho0^{|g|-mtrunc-1} |p|^{mtrunc+1}; the original remainder becomes r rho0^{m-mtrunc} |p|^{mtrunc+1}.
    (Fix of 2026-10-01: the previous version discarded the (mtrunc, m] terms without bounding them.)"""
    def __init__(self, tm, mtrunc):
        assert mtrunc <= tm.m
        self.m = mtrunc; self.mfull = tm.m; self.rho0 = tm.rho0; self.r = tm.r
        gammas = {}; abs_ = {}; rows = []
        for exp, co in tm.poly.terms():
            g = gam(exp); ab = (int(exp[0]), int(exp[1]))
            gammas.setdefault(g, len(gammas)); abs_.setdefault(ab, len(abs_)); rows.append((g, ab, co))
        self.gammas = list(gammas.keys()); self.abs = list(abs_.keys())
        Mx = arb_mat(len(self.gammas), len(self.abs))
        for g, ab, co in rows: Mx[gammas[g], abs_[ab]] += arb(co.p)/arb(co.q)
        self.Mx = Mx
        self.amax = max(ab[0] for ab in self.abs); self.bmax = max(ab[1] for ab in self.abs)
        r0 = arb(self.rho0.numerator)/arb(self.rho0.denominator)
        self.fold = {g: (arb(sigma(g).numerator)/arb(sigma(g).denominator)) * r0**(sum(g) - mtrunc - 1)
                     for g in self.gammas if sum(g) > mtrunc}
        self.rscale = r0**(self.mfull - mtrunc)
    def contract(self, Ic, Is, modek, derk, model, derl):
        ck, fuk, fvk = deriv_desc(modek, derk); cl, ful, fvl = deriv_desc(model, derl)
        const = ck*cl
        Iu = [factor_int(a, fuk, ful, Ic, Is) for a in range(self.amax+1)]
        Iv = [factor_int(b, fvk, fvl, Ic, Is) for b in range(self.bmax+1)]
        vec = arb_mat(len(self.abs), 1)
        for idx, (a,b) in enumerate(self.abs): vec[idx,0] = Iu[a]*Iv[b]
        res = self.Mx * vec
        out = {}; folded = arb(0)
        for i, g in enumerate(self.gammas):
            val = const*res[i,0]
            if sum(g) <= self.m: out[g] = val
            else: folded += abs(val) * self.fold[g]
        # |int_{[-1,1]^2} R Theta| <= r |p|^{m+1} * 4 |const| (|trig| <= 1, area 4)
        rem = abs(const) * 4 * arb(self.r.numerator)/arb(self.r.denominator) * self.rscale + folded
        return out, rem

if __name__ == "__main__":
    T = build_tensors(6, Fr(12,100), mLL=8)
    Ic, Is = trig_tables(20, 24)
    # checks: int_{-1}^{1} sin^2(pi(x+1)/2) dx = 1 ; via factor_int a=0
    print("int sin^2 =", factor_int(0, ('s',1), ('s',1), Ic, Is), " int x sin(pi(x+1)/2) =", Is[1][1], " exact 0")
    print("int x^2 cos(2 pi (x+1)/2)=", Ic[2][2], " exact 2/(pi^2)*? :", 2/arb.pi()**2)
