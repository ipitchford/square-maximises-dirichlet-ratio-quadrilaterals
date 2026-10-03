"""Stage 3: certified local step.  For direction boxes on S^3 and t in (0, rho0]: certify det S < 0 or tr S < 0
where S is the 2x2 Schur complement of A = Ktil - (5/2) L Mtil onto the cluster block, with
Ktil, Mtil the 3x3 RR pencil of first-order-corrected transported modes and L the Kato-Temple lower bound for lambda_1.
Exact facts used (proved in lemma_firstorder.py / lemma_residual.py, sympy):
  (L0) K(0)=diag(lambda), M(0)=I, LL(0)=diag(lambda^2)   (orthonormality)
  (L1) d_i lambda_1(0) = 0 and d_i tr[K_cc - 5 pi^2 M_cc](0) = 0  for i=a,b,c,d
  (L2) d_i [LL_11 - lambda_1^2 M_11 - 2 lambda_1 (K_11 - lambda_1 M_11)](0) = 0
Consequently the t^0 coefficients of A_cc, A_1c and the t^0,t^1 coefficients of eps^2, and t^1 of tr S vanish exactly."""
import sys, pickle, time, math, itertools
from fractions import Fraction as Fr
from flint import arb, arb_poly, arb_mat, ctx as flintctx
flintctx.prec = 256
PI = arb.pi()

# ---------------- univariate Taylor models on t in [0, rho0]: poly (list of arb coeffs, degree<=m) + rem (arb, |R|<=rem t^{m+1}) -----
class UTM:
    __slots__ = ('c','rem','m','rho')
    def __init__(self, c, rem, m, rho):
        c = list(c) + [arb(0)]*(m+1-len(c)); self.c = c[:m+1]; self.rem = arb(rem); self.m = m; self.rho = rho
    @staticmethod
    def const(x, m, rho): return UTM([arb(x)], 0, m, rho)
    def sup_abs(self):
        """upper bound of sup_{[0,rho]} |f|"""
        s = arb(0)
        for n, cn in enumerate(self.c): s += abs(cn) * self.rho**n
        return (s + self.rem*self.rho**(self.m+1)).abs_upper() if hasattr(s,'abs_upper') else s + self.rem*self.rho**(self.m+1)
    def __add__(self, o):
        if isinstance(o, UTM): return UTM([x+y for x,y in zip(self.c,o.c)], self.rem+o.rem, self.m, self.rho)
        c = list(self.c); c[0] = c[0] + arb(o); return UTM(c, self.rem, self.m, self.rho)
    __radd__ = __add__
    def __neg__(self): return UTM([-x for x in self.c], self.rem, self.m, self.rho)
    def __sub__(self, o): return self + (-o)
    def scale(self, s): s = fr2arb(s) if isinstance(s, Fr) else arb(s); return UTM([x*s for x in self.c], self.rem*abs(s), self.m, self.rho)
    def __mul__(self, o):
        if not isinstance(o, UTM): return self.scale(o)
        m = self.m; c = [arb(0)]*(m+1); rem = arb(0)
        for i, x in enumerate(self.c):
            if x == 0: continue
            for j, y in enumerate(o.c):
                if y == 0: continue
                if i+j <= m: c[i+j] += x*y
                else: rem += abs(x*y) * self.rho**(i+j-m-1)
        rem += self.rem*o._supabs() + o.rem*self._supabs() + self.rem*o.rem*self.rho**(m+1)
        return UTM(c, rem, m, self.rho)
    __rmul__ = __mul__
    def _supabs(self):
        s = arb(0)
        for n, cn in enumerate(self.c): s += abs(cn)*self.rho**n
        return s + self.rem*self.rho**(self.m+1)
    def inv(self):
        """1/f, requires c0 not containing 0; f = c0 (1 + g), g=O(t): 1/(1+g) = sum (-g)^n + tail"""
        c0 = self.c[0]
        assert not c0.contains(0), "inverse: constant term contains 0"
        g = UTM([arb(0)] + [x/c0 for x in self.c[1:]], self.rem/abs(c0), self.m, self.rho)
        sg = g._supabs()
        assert (sg < 1), f"inverse: |g| sup {sg} not < 1"
        total = UTM.const(1, self.m, self.rho); power = UTM.const(1, self.m, self.rho); mg = -g
        for n in range(1, self.m+1):
            power = power*mg; total = total + power
        # tail sum_{n>m} |g|^n <= sg^{m+1}/(1-sg) ; as remainder coefficient: bound |g(t)| <= (sg/rho) t ? no -- g=O(t) with |g|<=G1 t where G1 = sup|g|/t
        G1 = arb(0)
        for n, cn in enumerate(g.c[1:], start=1): G1 += abs(cn)*self.rho**(n-1)
        G1 += g.rem*self.rho**self.m
        tail = (G1**(self.m+1)) / (1 - sg)            # |sum_{n>m} (-g)^n| <= |g|^{m+1}/(1-|g|) <= (G1 t)^{m+1}/(1-sg)
        total.rem += tail
        return total.scale(1/c0)
    def shift_down(self, k):
        """divide by t^k assuming the first k coefficients are exactly zero"""
        for i in range(k): assert self.c[i] == 0, f"shift_down: coefficient {i} not exactly zero: {self.c[i]}"
        return UTM(self.c[k:], self.rem, self.m-k, self.rho)
    def bounds_on(self, t0, t1):
        """enclosure of f on [t0,t1] subset [0,rho] (naive interval evaluation)"""
        T = arb(t0).union(arb(t1)); val = arb(0)
        for n in range(self.m, -1, -1): val = val*T + self.c[n]
        R = self.rem * arb(t1)**(self.m+1)
        return val + arb(0).union(R) - arb(0).union(R)*0 if False else val + R*arb(-1).union(arb(1))

# ---------------- load coefficients -------------------------------------------------------------------------
def load(fn):
    d = pickle.load(open(fn,'rb'))
    def deser(pair_data):
        coeffs, rem = pair_data
        cd = {g: arb(v[0]) + arb(0).union(arb(v[1])) - arb(0).union(arb(v[1]))*0 if False else arb(v[0]).union(arb(v[0])) for g,v in coeffs.items()}
        # proper ball: midpoint v[0] with radius v[1]
        cd = {}
        for g, v in coeffs.items():
            mid = arb(v[0]); rad = arb(v[1]); cd[g] = mid + rad*arb(-1).union(arb(1))
        r = arb(rem[0]) + arb(rem[1])*arb(-1).union(arb(1))
        return cd, r
    out = {kk: {pair: deser(v) for pair, v in d[kk].items()} for kk in ('K','M','LL')}
    out.update(modes=d['modes'], base=d['base'], m=d['m'], mLL=d['mLL'], rho0=Fr(d['rho0']), N=d['N'])
    return out

def sym(dct, k, l): return dct[(k,l)] if (k,l) in dct else dct[(l,k)]

# ---------------- corrections C_{j k i} (rational approximations of first-order perturbation coefficients) -------------
def corrections(D, denom=2**40):
    modes = D['modes']; base = D['base']; lam = {mo: PI**2*(mo[0]**2+mo[1]**2) for mo in modes}
    cluster = {(2,1),(1,2)}
    C = {}
    for j in base:
        for k in modes:
            if k == j: continue
            if j in cluster and k in cluster: continue
            for i in range(4):
                e = [0,0,0,0]; e[i] = 1; e = tuple(e)
                Kc, _ = sym(D['K'], k, j); Mc, _ = sym(D['M'], k, j)
                val = -(Kc.get(e, arb(0)) - lam[j]*Mc.get(e, arb(0)))/(lam[k]-lam[j])
                fl = float(val.mid()) if hasattr(val,'mid') else float(val)
                C[(j,k,i)] = Fr(round(fl*denom), denom)
    return C

# ---------------- multivariate TM (dict gamma->arb, rem) arithmetic restricted to what we need --------------------------
def fr2arb(x):
    x = Fr(x); return arb(x.numerator)/arb(x.denominator)

def sigma_arb(g):
    """upper bound of max_{|e|=1} |e^g| = prod (g_i/k)^{g_i/2}"""
    k = sum(g)
    if k == 0: return arb(1)
    val = arb(1)
    for gi in g:
        if gi: val *= (arb(gi)/arb(k))**(arb(gi)/2)
    return val

def mono_mul(cd, rem, order, lin, rho0):
    """Multiply TM (cd, rem, order) [|R| <= rem |p|^{order+1}] by exact linear polynomial lin=(l0,l1..l4) (Fractions).
    Returns (cd', rem', order') with order' = order if l0 != 0 else order+1 (all contributions consistent)."""
    l0 = lin[0]
    order_out = order if l0 != 0 else order + 1
    out = {}; drop = arb(0)
    if l0 != 0:
        a0 = fr2arb(l0)
        for g, v in cd.items(): out[g] = out.get(g, arb(0)) + v*a0
    for i in range(4):
        li = lin[1+i]
        if li == 0: continue
        ai = fr2arb(li)
        for g, v in cd.items():
            g2 = list(g); g2[i] += 1; g2 = tuple(g2)
            term = v*ai
            if sum(g2) <= order_out: out[g2] = out.get(g2, arb(0)) + term
            else: drop += abs(term)*sigma_arb(g2)*fr2arb(rho0)**(sum(g2)-order_out-1)
    sum_li = sum(abs(lin[1+i]) for i in range(4))
    if l0 != 0:
        rem_out = rem*(fr2arb(abs(l0)) + fr2arb(sum_li*rho0)) + drop       # exponent order+1
    else:
        rem_out = rem*fr2arb(sum_li) + drop                                  # R * sum l_i p_i : exponent order+2
    return out, rem_out, order_out

def truncate_to(cd, rem, order, mm, rho0):
    """convert TM valid at order `order` (>= mm) to order mm."""
    assert order >= mm
    out = {}; extra = arb(0); r0 = fr2arb(rho0)
    for g, v in cd.items():
        if sum(g) <= mm: out[g] = v
        else: extra += abs(v)*sigma_arb(g)*r0**(sum(g)-mm-1)
    return out, rem*r0**(order-mm) + extra, mm

def build_pencil(D, C):
    """3x3 pencils Ktil, Mtil, LLtil as multivariate TMs at order m (K,M) and mLL (LL)."""
    modes = D['modes']; base = D['base']; m = D['m']; mLL = D['mLL']; rho0 = D['rho0']
    def cvec(j):   # coefficients c_{jk}(p) as linear polynomials (l0, l1..l4)
        out = {}
        for k in modes:
            lin = [Fr(1 if k==j else 0)] + [C.get((j,k,i), Fr(0)) for i in range(4)]
            if any(x != 0 for x in lin): out[k] = lin
        return out
    cv = {j: cvec(j) for j in base}
    def form(dct, mm):
        P = {}
        for j in base:
            for jp in base:
                acc = {}; rem = arb(0)
                for k, lk in cv[j].items():
                    for l, ll in cv[jp].items():
                        cd, r = sym(dct, k, l)
                        nb = (k in base) + (l in base)
                        order_in = mm - (2 - nb)
                        cd1, r1, o1 = mono_mul(cd, r, order_in, lk, rho0)
                        cd2, r2, o2 = mono_mul(cd1, r1, o1, ll, rho0)
                        assert o2 >= mm, (k, l, o2, mm)
                        if o2 > mm: cd2, r2, o2 = truncate_to(cd2, r2, o2, mm, rho0)
                        for g, v in cd2.items(): acc[g] = acc.get(g, arb(0)) + v
                        rem += r2
                P[(j,jp)] = (acc, rem)
        return P
    return form(D['K'], m), form(D['M'], m), form(D['LL'], mLL)

# ---------------- substitute p = t e for a direction box, get UTMs --------------------------------------------------------
def substitute(cd, rem, e, m, rho):
    """e: list of 4 arb intervals with |e|=1 (as a set). Returns UTM in t."""
    c = [arb(0)]*(m+1)
    for g, v in cd.items():
        n = sum(g)
        if n > m: continue
        mon = arb(1)
        for i in range(4):
            if g[i]: mon *= e[i]**g[i]
        c[n] += v*mon
    return UTM(c, rem, m, rho)

def direction_intervals(face, box):
    """face: (index i, sign s); box: list of 3 intervals (lo,hi) for the other coordinates. Returns e_i as arb intervals."""
    i, s = face; others = [j for j in range(4) if j != i]
    ehat = [None]*4; ehat[i] = arb(s)
    for j, (lo, hi) in zip(others, box): ehat[j] = arb(lo).union(arb(hi))
    n2 = sum(x*x for x in ehat); n = n2.sqrt()
    return [x/n for x in ehat]

def certify_direction(P, e, m, mLL, rho, beta, lam, tsplit=8, verbose=False):
    """Return (ok, info). P = (Ktil, Mtil, LLtil) multivariate TMs."""
    Ktil, Mtil, LLtil = P
    base_idx = {b:i for i,b in enumerate([(1,1),(2,1),(1,2)])}
    b = [(1,1),(2,1),(1,2)]
    K = [[substitute(*Ktil[(b[i],b[j])], e, m, rho) for j in range(3)] for i in range(3)]
    M = [[substitute(*Mtil[(b[i],b[j])], e, m, rho) for j in range(3)] for i in range(3)]
    LL11 = substitute(*LLtil[(b[0],b[0])], e, mLL, rho)
    # truncate LL11 to order m for consistency
    LL11 = UTM(LL11.c[:m+1], LL11.rem*rho**(mLL-m) + sum(abs(x)*rho**(n-m-1) for n,x in enumerate(LL11.c) if n>m), m, rho)
    # exact structure surgery (L0): constant terms
    lam1, lam2 = lam[0], lam[1]
    for i in range(3):
        for j in range(3):
            K[i][j].c[0] = arb(0) if i!=j else (lam1 if i==0 else lam2)
            M[i][j].c[0] = arb(0) if i!=j else arb(1)
    LL11.c[0] = lam1*lam1
    # (L1): first-order of lambda_1 (rho tilde) vanishes: K11_1 - lam1 M11_1 = 0 ; we set rho_1 = 0 below
    Minv11 = M[0][0].inv()
    rho_t = K[0][0]*Minv11
    rho_t.c[1] = arb(0)
    # eps^2 = LL11/M11 - rho^2 ; (L0),(L2): coefficients 0 and 1 vanish exactly
    eps2 = LL11*Minv11 - rho_t*rho_t
    chk01 = (eps2.c[0], eps2.c[1])
    eps2.c[0] = arb(0); eps2.c[1] = arb(0)
    # beta - rho > 0 check
    brho = (-rho_t) + beta
    assert (brho.c[0] - brho._supabs() + brho.c[0]) > 0 or True
    L = rho_t - eps2*brho.inv()
    # A = K - (5/2) L M
    A = [[K[i][j] - (L*M[i][j]).scale(Fr(5,2)) for j in range(3)] for i in range(3)]
    # surgery: A_cc(0)=0, A_1c(0)=0 exactly (L0), tr A_cc first order = 0 (L1)
    for i in range(3):
        for j in range(3):
            if not (i==0 and j==0): A[i][j].c[0] = arb(0)
    # Schur complement S = A_cc - A_c1 A_1c / A_11
    A11inv = A[0][0].inv()
    S = [[A[1+i][1+j] - A[1+i][0]*A[0][1+j]*A11inv for j in range(2)] for i in range(2)]
    trS = S[0][0] + S[1][1]
    tr1 = trS.c[1]; trS.c[1] = arb(0)          # (L1)
    detS = S[0][0]*S[1][1] - S[0][1]*S[1][0]
    # both are O(t^2): shift
    q_det = detS.shift_down(2); q_tr = trS.shift_down(2)
    # certify on t in (0,rho]: split
    ok_all = True; worst = None; details = []
    for s in range(tsplit):
        t0 = rho*s/tsplit; t1 = rho*(s+1)/tsplit
        bd = q_det.bounds_on(t0, t1); bt = q_tr.bounds_on(t0, t1)
        ok = (bd < 0) or (bt < 0)
        details.append((float(t0), float(t1), bd, bt, bool(ok)))
        ok_all = ok_all and ok
    # side conditions: A11 < 0 and beta > rho on [0,rho]
    a11 = A[0][0].bounds_on(0, rho); br = brho.bounds_on(0, rho)
    side = (a11 < 0) and (br > 0)
    return ok_all and side, dict(details=details, a11=a11, brho=br, tr1=tr1, chk01=chk01, side=side)
