"""Independent re-certification of global boxes with the affine-transport certificate (exact rational reference integrals).
Different code, different test functions (own float RR), exact rational arithmetic (python-flint fmpq_mpoly) for every integral
and for the box polynomials; the only floats are used to CHOOSE test functions and the flux potential (any choice is admissible).
Certificate: xi(Q_p) < 5/2 for all p in the box (where Q_p is convex) if, with D = det A_1 det A_2 > 0,
  At := D^2 (K - (5/2) L M)  [cleared denominators]  has At11 < 0 and S' = -At11 At_cc + At_c1 At_1c has a negative eigenvalue,
  plus beta*D*M11 - Kt11 > 0, M > 0 (Sylvester minors).  Lehmann-Goerisch: L = (beta A1 - A0)/(beta A2 - A1)."""
import sys, json, math, time
from fractions import Fraction as Fr
import numpy as np
from flint import fmpq_mpoly_ctx, fmpq, fmpq_mpoly, Ordering, arb
sys.path.insert(0, '../explore')
from lg_proto import Basis as FBasis, psi_basis, quad as fquad


C2 = fmpq_mpoly_ctx.get(('x', 'y'), Ordering.lex); X, Y = C2.gens()
C3 = fmpq_mpoly_ctx.get(('x', 'y', 't'), Ordering.lex); X3, Y3, T3 = C3.gens()
C4 = fmpq_mpoly_ctx.get(('d1', 'd2', 'd3', 'd4'), Ordering.lex); DV4 = C4.gens()
def q(v): v = Fr(v); return fmpq(v.numerator, v.denominator)

def verts_exact(p):
    a, b, c, d = [Fr(x) for x in p]; h = Fr(1, 2)
    return [(h*(a+d-b), h*(-a+d+c)), (1+h*(-a+d+b), h*(-a-d-c)), (1+h*(-a-d-b), 1+h*(a-d+c)), (h*(a-d+b), 1+h*(a+d-c))]
# vertex derivative table d v_j / d p_i (exact)
DVT = [[(Fr(1,2), Fr(-1,2)), (Fr(-1,2), 0), (0, Fr(1,2)), (Fr(1,2), Fr(1,2))],
       [(Fr(-1,2), Fr(-1,2)), (Fr(1,2), 0), (0, Fr(-1,2)), (Fr(1,2), Fr(-1,2))],
       [(Fr(-1,2), Fr(1,2)), (Fr(-1,2), 0), (0, Fr(1,2)), (Fr(-1,2), Fr(-1,2))],
       [(Fr(1,2), Fr(1,2)), (Fr(1,2), 0), (0, Fr(-1,2)), (Fr(-1,2), Fr(1,2))]]   # [vertex][param] = (dx, dy)

def legendre_coeffs(k):
    """exact rational coefficients of P_k (ascending powers)"""
    P = [[Fr(1)], [Fr(0), Fr(1)]]
    for n in range(1, k):
        a = [Fr(0)] + [Fr(2*n+1, n+1)*c for c in P[n]]; b = [Fr(n, n+1)*c for c in P[n-1]] + [Fr(0)]*2
        P.append([a[i] - (b[i] if i < len(b) else 0) for i in range(len(a))])
    return P[k]
def leg_poly(k, var, c, s):
    """P_k((var - c)/s) as an fmpq_mpoly"""
    z = (var - q(c)) * q(Fr(1)/Fr(s)); out = var*0; zp = var*0 + 1
    for coef in legendre_coeffs(k):
        out = out + zp*q(coef); zp = zp*z
    return out

def tri_integral(f, P0, P1, P2):
    """exact integral of polynomial f(x,y) over triangle P0P1P2"""
    xi, eta = X, Y   # reuse context vars as (xi, eta)
    gx = q(P0[0]) + xi*q(P1[0]-P0[0]) + eta*q(P2[0]-P0[0]); gy = q(P0[1]) + xi*q(P1[1]-P0[1]) + eta*q(P2[1]-P0[1])
    g = f.compose(gx, gy)
    tot = Fr(0)
    for exp, co in g.terms():
        i, j = int(exp[0]), int(exp[1])
        tot += Fr(int(co.p), int(co.q)) * Fr(math.factorial(i)*math.factorial(j), math.factorial(i+j+2))
    det = abs((P1[0]-P0[0])*(P2[1]-P0[1]) - (P1[1]-P0[1])*(P2[0]-P0[0]))
    return tot*det

class XCert:
    def __init__(self, p0, n=8, kpsi=12, shift=None):
        self.p0 = [Fr(v) for v in p0]; V = verts_exact(self.p0)
        if shift is None:
            # choose the diagonal (v1v3: shift 0, v2v4: shift 1) maximising the smaller triangle area at p0
            def ar(P0, P1, P2): return abs((P1[0]-P0[0])*(P2[1]-P0[1]) - (P1[1]-P0[1])*(P2[0]-P0[0]))
            m0 = min(ar(V[0], V[1], V[2]), ar(V[0], V[2], V[3])); m1 = min(ar(V[1], V[2], V[3]), ar(V[1], V[3], V[0]))
            shift = 0 if m0 >= m1 else 1
        self.shift = shift
        V = V[shift:] + V[:shift]; self.V = V; self.n = n
        fV = np.array([[float(a), float(b)] for a, b in V])
        lo = fV.min(0); hi = fV.max(0)
        self.c = [Fr(float((lo[i]+hi[i])/2)) for i in range(2)]; self.s = [Fr(float((hi[i]-lo[i])/2)) for i in range(2)]
        # float RR with the same scaled-Legendre bubble basis (lg_proto.Basis uses identical c, s)
        FB = FBasis(fV, n); FB.c = np.array([float(v) for v in self.c]); FB.s = np.array([float(v) for v in self.s])
        xq, wq = fquad(fV, n + 8)
        Phi, Gx, Gy, Lap = FB.eval(xq)
        K = (Gx*wq)@Gx.T + (Gy*wq)@Gy.T; M = (Phi*wq)@Phi.T
        lam, Qm = np.linalg.eigh(M); keep = lam > 1e-12*lam.max(); Tm = Qm[:, keep]/np.sqrt(lam[keep])
        th, U = np.linalg.eigh(Tm.T@K@Tm); Cf = (Tm@U)[:, :3]
        self.th_f = th[:3]
        # exact polynomials: bubble and P_k
        b0 = X*0 + 1
        for i in range(4):
            a = V[i]; e = (V[(i+1) % 4][0]-a[0], V[(i+1) % 4][1]-a[1])
            b0 = b0 * (q(e[0])*(Y - q(a[1])) - q(e[1])*(X - q(a[0])))
        Lx = [leg_poly(k, X, self.c[0], self.s[0]) for k in range(n+1)]; Ly = [leg_poly(k, Y, self.c[1], self.s[1]) for k in range(n+1)]
        basis = [Lx[i]*Ly[j] for (i, j) in FB.idx]
        self.u = []
        for k in range(3):
            P = X*0
            for m, bm in enumerate(basis): P = P + bm*q(Fr(float(Cf[m, k])))
            self.u.append(b0*P)
        self.gx = [u.derivative(0) for u in self.u]; self.gy = [u.derivative(1) for u in self.u]
        # tau^ = F e, F(x) = int_0^{t*(x)} u1(x + t e) dt, base line through v2, v4; e = v3 - v1
        e = (V[2][0]-V[0][0], V[2][1]-V[0][1]); self.e = e
        nb = (-(V[3][1]-V[1][1]), V[3][0]-V[1][0])          # normal of the base line v2v4
        # t*(x): nb.(x + t e - v2) = 0  ->  t* = -nb.(x - v2)/nb.e
        nbe = nb[0]*e[0] + nb[1]*e[1]
        tstar = (X - q(V[1][0]))*q(-nb[0]/nbe) + (Y - q(V[1][1]))*q(-nb[1]/nbe)
        u13 = self.u[0].compose(X3 + T3*q(e[0]), Y3 + T3*q(e[1]), ctx=C3) if False else None
        # compose into 3 variables: substitute x -> x + t e_x, y -> y + t e_y
        d = {}
        for exp, co in self.u[0].terms(): d[(int(exp[0]), int(exp[1]), 0)] = co
        u1_3 = C3.from_dict(d).compose(X3 + T3*q(e[0]), Y3 + T3*q(e[1]), T3)
        # antiderivative in t, evaluated at t*(x,y)
        anti = {}
        for exp, co in u1_3.terms():
            i, j, k = (int(v) for v in exp); anti[(i, j, k+1)] = co * fmpq(1, k+1)
        A3 = C3.from_dict(anti)
        ts3 = C3.from_dict({(int(a), int(b), 0): c for (a, b), c in [((ex[0], ex[1]), co) for ex, co in tstar.terms()]})
        F3 = A3.compose(X3, Y3, ts3)
        Fd = {}
        for exp, co in F3.terms(): Fd[(int(exp[0]), int(exp[1]))] = co
        F = C2.from_dict(Fd)
        self.tau = (F*q(e[0]), F*q(e[1]))
        # check div tau = -u1 exactly
        assert self.tau[0].derivative(0) + self.tau[1].derivative(1) + self.u[0] == 0, "div tau != -u1"
        # psi by float least squares (minimise ||tau + curl psi||^2 on Q0), then exact rational coefficients
        tau_f = np.array([[float(Fr(str(self.tau[a](q(Fr(float(x))), q(Fr(float(yv))))))) for a in range(2)] for x, yv in xq]) if False else None
        taux = [float(Fr(int(co.p), int(co.q))) for co in []]
        # float evaluation of tau at nodes via numpy polynomial evaluation
        def fe(poly, pts):
            out = np.zeros(len(pts))
            for exp, co in poly.terms():
                out += float(Fr(int(co.p), int(co.q))) * pts[:, 0]**int(exp[0]) * pts[:, 1]**int(exp[1])
            return out
        tfx = fe(self.tau[0], xq); tfy = fe(self.tau[1], xq)
        cx, cy = psi_basis(xq, np.array([float(v) for v in self.c]), np.array([float(v) for v in self.s]), kpsi)
        G = (cx*wq)@cx.T + (cy*wq)@cy.T; rhs = -((cx*wq)@tfx + (cy*wq)@tfy)
        coef = np.linalg.lstsq(G, rhs, rcond=None)[0]
        # exact psi with the same basis ordering as lg_proto.psi_basis
        psi = X*0; m = 0
        Lxk = [leg_poly(k, X, self.c[0], self.s[0]) for k in range(kpsi+1)]; Lyk = [leg_poly(k, Y, self.c[1], self.s[1]) for k in range(kpsi+1)]
        for i in range(kpsi+1):
            for j in range(kpsi+1-i):
                if i == 0 and j == 0: continue
                psi = psi + Lxk[i]*Lyk[j]*q(Fr(float(coef[m]))); m += 1
        self.curl = (psi.derivative(1), -psi.derivative(0))
        # exact reference integrals on K1=(v1,v2,v3), K2=(v1,v3,v4)
        self.tris = [(0, 1, 2), (0, 2, 3)]; self.ref = []
        g = [(self.gx[k], self.gy[k]) for k in range(3)]
        for (ia, ib, ic) in self.tris:
            P0, P1, P2 = V[ia], V[ib], V[ic]; R = {}
            R['I'] = {(a, b, i, j): tri_integral(g[i][a]*g[j][b], P0, P1, P2) for a in range(2) for b in range(2) for i in range(3) for j in range(3) if j >= i}
            R['J'] = {(i, j): tri_integral(self.u[i]*self.u[j], P0, P1, P2) for i in range(3) for j in range(3) if j >= i}
            R['TT'] = {(a, b): tri_integral(self.tau[a]*self.tau[b], P0, P1, P2) for a in range(2) for b in range(2)}
            R['TC'] = {(a, b): tri_integral(self.tau[a]*self.curl[b], P0, P1, P2) for a in range(2) for b in range(2)}
            R['CC'] = {(a, b): tri_integral(self.curl[a]*self.curl[b], P0, P1, P2) for a in range(2) for b in range(2)}
            MK0 = [[P1[0]-P0[0], P2[0]-P0[0]], [P1[1]-P0[1], P2[1]-P0[1]]]
            det0 = MK0[0][0]*MK0[1][1] - MK0[0][1]*MK0[1][0]
            R['MK0inv'] = [[MK0[1][1]/det0, -MK0[0][1]/det0], [-MK0[1][0]/det0, MK0[0][0]/det0]]
            R['idx'] = (ia, ib, ic)
            self.ref.append(R)

    def box_polys(self, beta):
        """exact polynomials in delta (4 vars): D, Kt (=D K), M, D*A2, and the cleared criterion matrix At (3x3)."""
        one = C4.from_dict({(0, 0, 0, 0): fmpq(1)})
        dets = []; Kparts = []; Mparts = []; A2parts = []
        for R in self.ref:
            ia, ib, ic = R['idx']
            # M_K(p) = M_K(p0) + sum_i delta_i dM_i ; A_K = M_K(p) MK0inv
            def edge(j, k, comp):
                base = self.V[k][comp] - self.V[j][comp]
                poly = one*q(base)
                DVS = DVT[self.shift:] + DVT[:self.shift]
                for i in range(4): poly = poly + DV4[i]*q(DVS[k][i][comp] - DVS[j][i][comp])
                return poly
            MK = [[edge(ia, ib, 0), edge(ia, ic, 0)], [edge(ia, ib, 1), edge(ia, ic, 1)]]
            Mi = R['MK0inv']
            A = [[MK[r][0]*q(Mi[0][c]) + MK[r][1]*q(Mi[1][c]) for c in range(2)] for r in range(2)]
            det = A[0][0]*A[1][1] - A[0][1]*A[1][0]
            adj = [[A[1][1], -A[0][1]], [-A[1][0], A[0][0]]]
            GA = [[adj[r][0]*adj[s][0] + adj[r][1]*adj[s][1] for s in range(2)] for r in range(2)]   # adj adj^T
            G = [[A[0][r]*A[0][s] + A[1][r]*A[1][s] for s in range(2)] for r in range(2)]           # A^T A
            Kp = {}
            for i in range(3):
                for j in range(i, 3):
                    Kp[(i, j)] = sum((GA[a][b]*q(R['I'][(a, b, i, j)]) for a in range(2) for b in range(2)), one*0)
            Mp = {(i, j): det*q(R['J'][(i, j)]) for i in range(3) for j in range(i, 3)}
            tt = sum((G[a][b]*q(R['TT'][(a, b)]) for a in range(2) for b in range(2)), one*0)
            tc = sum((G[a][b]*q(R['TC'][(a, b)]) for a in range(2) for b in range(2)), one*0)
            cc = sum((G[a][b]*q(R['CC'][(a, b)]) for a in range(2) for b in range(2)), one*0)
            dets.append(det); Kparts.append(Kp); Mparts.append(Mp); A2parts.append((tt, tc, cc))
        D = dets[0]*dets[1]
        Kt = {k: Kparts[0][k]*dets[1] + Kparts[1][k]*dets[0] for k in Kparts[0]}            # D*K
        Mm = {k: Mparts[0][k] + Mparts[1][k] for k in Mparts[0]}
        (tt1, tc1, cc1), (tt2, tc2, cc2) = A2parts
        DA2 = (dets[0]*tt1 + tc1*2)*D + cc1*dets[1] + (dets[1]*tt2 + tc2*2)*D + cc2*dets[0]   # D*A2
        kap = fmpq(5, 2)
        if beta == "KS":
            # pointwise Krahn-Szego shift beta(p) = Kq / S(p) <= 2 pi j01^2 / |Q_p| <= lambda_2(Q_p), S(p) = |Q_p| = 1 - a^2 - d^2;
            # all LG quantities multiplied by S(p) > 0 (checked by the caller), so signs are preserved.
            Kq = fmpq(363368, 10000)                    # < 2*pi*j01^2 = 36.33681...
            S = 1 - (DV4[0] + q(self.p0[0]))**2 - (DV4[3] + q(self.p0[3]))**2
            g1 = DA2*Kq - D*Mm[(0, 0)]*S                # S D (beta A2 - A1)
            g2 = D*Mm[(0, 0)]*Kq - Kt[(0, 0)]*S         # S D (beta A1 - A0)
            At = {k: Kt[k]*g1 - Mm[k]*g2*D*kap for k in Kt}
            return dict(D=D, dets=dets, Kt=Kt, M=Mm, DA2=DA2, g1=g1, g2=g2, At=At, S=S)
        b = q(beta)
        g1 = DA2*b - D*Mm[(0, 0)]                      # D (beta A2 - A1)
        g2 = D*Mm[(0, 0)]*b - Kt[(0, 0)]               # D (beta A1 - A0)
        At = {k: Kt[k]*g1 - Mm[k]*g2*D*kap for k in Kt}   # D^2 (K - kappa L M)
        return dict(D=D, dets=dets, Kt=Kt, M=Mm, DA2=DA2, g1=g1, g2=g2, At=At)

def upper_bound(poly, r):
    """rigorous upper bound of poly(delta) over |delta_i| <= r_i (exact rational coefficients): c0 + sum |c| r^gamma"""
    tot = Fr(0)
    for exp, co in poly.terms():
        c = Fr(int(co.p), int(co.q)); g = [int(v) for v in exp]
        if sum(g) == 0: tot += c
        else: tot += abs(c) * math.prod(Fr(r[i])**g[i] for i in range(4))
    return tot
def lower_bound(poly, r): return -upper_bound(-poly, r)

def beta_box(lo, hi):
    """rigorous lower bound for lambda_2 over the box (convex Q_p): max(Krahn-Szego, axis rectangle, AC with rectangle lambda_1)."""
    lo = [Fr(v) for v in lo]; hi = [Fr(v) for v in hi]
    a2min = min(x*x for x in (lo[0], hi[0])) if not (lo[0] <= 0 <= hi[0]) else Fr(0)
    d2min = min(x*x for x in (lo[3], hi[3])) if not (lo[3] <= 0 <= hi[3]) else Fr(0)
    area_max = 1 - a2min - d2min
    j01 = arb("2.404825557695772")                  # lower bound of j_{0,1} = 2.40482555769577276...
    ks = 2*arb.pi()*j01*j01/arb(area_max.numerator)*arb(area_max.denominator)
    mabs = lambda l, h: max(abs(l), abs(h))
    Aup = 1 - lo[0] + max(mabs(lo[1], hi[1]), mabs(lo[3], hi[3])); Bup = 1 + hi[0] + max(mabs(lo[2], hi[2]), mabs(lo[3], hi[3]))
    A2 = arb(Aup.numerator)/arb(Aup.denominator); A2 = A2*A2; B2 = arb(Bup.numerator)/arb(Bup.denominator); B2 = B2*B2
    pi2 = arb.pi()**2
    rect = min((pi2*(4/A2 + 1/B2)).lower(), (pi2*(1/A2 + 4/B2)).lower())
    l1rect = pi2*(1/A2 + 1/B2)
    Dmax2 = A2 + B2                                  # diameter^2 <= A^2 + B^2 (axis rectangle)
    ac = l1rect + 3*pi2/Dmax2
    best = max(float(ks.lower()), float(rect), float(ac.lower()))
    return Fr(best).limit_denominator(10**6) - Fr(1, 10**6)   # rational slightly below

EXACT_FALLBACK = False
def certify(lo, hi, n=8, kpsi=12, verbose=False):
    lo = [Fr(v) for v in lo]; hi = [Fr(v) for v in hi]
    # expand around a SIMPLE rational point near the midpoint (small denominators keep the exact arithmetic cheap);
    # the radius is enlarged so that p0 + [-r, r] still contains the whole box (validity is unaffected).
    p0 = []
    for l, h in zip(lo, hi):
        m = (l+h)/2; w = h-l
        q = m.limit_denominator(max(16, int(64/w) if w > 0 else 10**6)) if w > 0 else m
        p0.append(q)
    r = [max(h-c, c-l) for l, h, c in zip(lo, hi, p0)]
    r = [Fr(math.ceil(x * 2**24), 2**24) for x in r]          # round radii UP to simple dyadics (still cover the box)
    t0 = time.time(); xc = XCert(p0, n, kpsi); t1 = time.time()
    beta = beta_box(lo, hi)
    P = xc.box_polys(beta); t2 = time.time()
    res = {}
    res['D>0'] = lower_bound(P['dets'][0], r) > 0 and lower_bound(P['dets'][1], r) > 0
    res['g2>0'] = lower_bound(P['g2'], r) > 0                  # beta A1 > A0  (times D>0)
    M = P['M']
    m11 = M[(0, 0)]; m2 = M[(0, 0)]*M[(1, 1)] - M[(0, 1)]*M[(0, 1)]
    m3 = (M[(0, 0)]*(M[(1, 1)]*M[(2, 2)] - M[(1, 2)]*M[(1, 2)]) - M[(0, 1)]*(M[(0, 1)]*M[(2, 2)] - M[(1, 2)]*M[(0, 2)])
          + M[(0, 2)]*(M[(0, 1)]*M[(1, 2)] - M[(1, 1)]*M[(0, 2)]))
    res['M>0'] = all(lower_bound(z, r) > 0 for z in (m11, m2, m3))
    At = P['At']; a11 = At[(0, 0)]
    res['A11<0'] = upper_bound(a11, r) < 0
    # interval enclosures of the degree-12 entries (exact rationals), then exact interval arithmetic for the Schur complement
    I = {k: (lower_bound(v, r), upper_bound(v, r)) for k, v in At.items()}
    def imul(x, y):
        c = [x[0]*y[0], x[0]*y[1], x[1]*y[0], x[1]*y[1]]; return (min(c), max(c))
    def iadd(x, y): return (x[0]+y[0], x[1]+y[1])
    def ineg(x): return (-x[1], -x[0])
    def isq(x): return (Fr(0) if x[0] <= 0 <= x[1] else min(x[0]**2, x[1]**2), max(x[0]**2, x[1]**2))
    na = ineg(I[(0, 0)])
    S11 = iadd(imul(na, I[(1, 1)]), isq(I[(0, 1)]))
    S22 = iadd(imul(na, I[(2, 2)]), isq(I[(0, 2)]))
    S12 = iadd(imul(na, I[(1, 2)]), imul(I[(0, 1)], I[(0, 2)]))
    neg = None
    for name, iv in (('S11', S11), ('S22', S22), ('trS', iadd(S11, S22))):
        if iv[1] < 0: neg = name + '(iv)'; break
    if neg is None:
        detS = iadd(imul(S11, S22), ineg(isq(S12)))
        if detS[1] < 0: neg = 'detS(iv)'
    if neg is None and EXACT_FALLBACK:
        # exact degree-24 fallback (no detS: too large)
        s11 = -a11*At[(1, 1)] + At[(0, 1)]*At[(0, 1)]
        s22 = -a11*At[(2, 2)] + At[(0, 2)]*At[(0, 2)]
        for name, poly in (('S11', s11), ('S22', s22), ('trS', s11 + s22)):
            if upper_bound(poly, r) < 0: neg = name; break
    res['S has neg'] = neg
    ok = res['D>0'] and res['g2>0'] and res['M>0'] and res['A11<0'] and neg is not None
    if verbose: print(f"  setup {t1-t0:.1f}s, polys {t2-t1:.1f}s, beta={float(beta):.3f}, th_f={np.round(xc.th_f,4)}, checks={res}")
    return ok, res, float(beta)

if __name__ == "__main__":
    # smoke test on a generic box
    ok, res, beta = certify([Fr(29,100), Fr(19,100), Fr(39,100), Fr(9,100)], [Fr(31,100), Fr(21,100), Fr(41,100), Fr(11,100)], verbose=True)
    print("generic box certified:", ok)


# ------------------------------------------------------------------------------------------------------------------
# Adaptive range subdivision: ONE certificate (test functions, flux, beta) for the whole box; only the rigorous range
# bounding of the exact criterion polynomials is subdivided (Taylor shift to each sub-box centre, |coef| r^gamma bound,
# exact interval Schur complement).  Valid because all quantities are exact polynomials in delta = p - p0 and beta is a
# lower bound for lambda_2 on the whole box.
def _dyadic_up(x): return Fr(math.ceil(x * 2**30), 2**30)

def certify_adaptive(lo, hi, n=8, kpsi=12, maxdepth=10, maxleaves=4000, verbose=False):
    lo = [Fr(v) for v in lo]; hi = [Fr(v) for v in hi]
    p0 = []
    for l, h in zip(lo, hi):
        m = (l+h)/2; w = h-l
        p0.append(m.limit_denominator(max(16, int(64/w))) if w > 0 else m)
    t0 = time.time(); xc = XCert(p0, n, kpsi); xc.p0 = p0; beta = beta_box(lo, hi)
    variants = []
    for bmode in (beta, "KS"):
        P = xc.box_polys(bmode); M = P['M']; At = P['At']
        m2 = M[(0, 0)]*M[(1, 1)] - M[(0, 1)]*M[(0, 1)]
        m3 = (M[(0, 0)]*(M[(1, 1)]*M[(2, 2)] - M[(1, 2)]*M[(1, 2)]) - M[(0, 1)]*(M[(0, 1)]*M[(2, 2)] - M[(1, 2)]*M[(0, 2)])
              + M[(0, 2)]*(M[(0, 1)]*M[(1, 2)] - M[(1, 1)]*M[(0, 2)]))
        pos = {'D1': P['dets'][0], 'D2': P['dets'][1], 'g2': P['g2'], 'M11': M[(0, 0)], 'M2': m2, 'M3': m3}
        if 'S' in P: pos['S'] = P['S']
        variants.append((pos, {k: At[k] for k in At}))
    stats = {'leaves': 0, 'failed': 0, 'corner_fail': False}
    # cheap exact pre-check at the centre and the 16 corners: if no variant passes at some point, range subdivision is futile
    import itertools as _it
    def point_ok(pos, ent, d):
        ev = lambda poly: poly(*[q(x) for x in d])
        if not all(ev(v) > 0 for v in pos.values()): return False
        a = {k: ev(v) for k, v in ent.items()}
        if not a[(0, 0)] < 0: return False
        S11 = -a[(0, 0)]*a[(1, 1)] + a[(0, 1)]**2; S22 = -a[(0, 0)]*a[(2, 2)] + a[(0, 2)]**2; S12 = -a[(0, 0)]*a[(1, 2)] + a[(0, 1)]*a[(0, 2)]
        return S11 < 0 or S22 < 0 or S11 + S22 < 0 or S11*S22 - S12*S12 < 0
    pts = [[Fr(0)]*4] + [[(l if b == 0 else h) - c for l, h, c, b in zip(lo, hi, p0, bits)] for bits in _it.product((0, 1), repeat=4)]
    for d in pts:
        if not any(point_ok(pos, ent, d) for pos, ent in variants):
            stats['corner_fail'] = True
            if verbose: print("  pointwise failure at a corner/centre: box must be split", flush=True)
            return False, stats, float(beta)
    def bounds_on(polys, c, r, want):
        out = {}
        for k, v in polys.items():
            sv = v.compose(*[DV4[i] + q(c[i]) for i in range(4)])
            out[k] = (lower_bound(sv, r) if 'lo' in want else None, upper_bound(sv, r) if 'up' in want else None)
        return out
    def leaf_ok(dl, dh):
        return any(leaf_ok_v(dl, dh, pos, ent) for pos, ent in variants)
    def leaf_ok_v(dl, dh, pos, ent):
        c = [(a+b)/2 for a, b in zip(dl, dh)]; r = [_dyadic_up((b-a)/2) for a, b in zip(dl, dh)]
        pb = bounds_on(pos, c, r, ('lo',))
        if not all(v[0] > 0 for v in pb.values()): return False
        I = bounds_on(ent, c, r, ('lo', 'up'))
        if not (I[(0, 0)][1] < 0): return False
        def imul(x, y):
            cc = [x[0]*y[0], x[0]*y[1], x[1]*y[0], x[1]*y[1]]; return (min(cc), max(cc))
        def iadd(x, y): return (x[0]+y[0], x[1]+y[1])
        def ineg(x): return (-x[1], -x[0])
        def isq(x): return (Fr(0) if x[0] <= 0 <= x[1] else min(x[0]**2, x[1]**2), max(x[0]**2, x[1]**2))
        na = ineg(I[(0, 0)])
        S11 = iadd(imul(na, I[(1, 1)]), isq(I[(0, 1)])); S22 = iadd(imul(na, I[(2, 2)]), isq(I[(0, 2)]))
        S12 = iadd(imul(na, I[(1, 2)]), imul(I[(0, 1)], I[(0, 2)]))
        if S11[1] < 0 or S22[1] < 0 or iadd(S11, S22)[1] < 0: return True
        return iadd(imul(S11, S22), ineg(isq(S12)))[1] < 0
    def rec(dl, dh, depth):
        if stats['leaves'] >= maxleaves: stats['failed'] += 1; return False
        if leaf_ok(dl, dh): stats['leaves'] += 1; return True
        if depth >= maxdepth: stats['leaves'] += 1; stats['failed'] += 1; return False
        w = [b-a for a, b in zip(dl, dh)]; d = w.index(max(w)); mid = (dl[d]+dh[d])/2
        h1 = list(dh); h1[d] = mid; l2 = list(dl); l2[d] = mid
        return rec(dl, h1, depth+1) and rec(l2, dh, depth+1)
    ok = rec([l-c for l, c in zip(lo, p0)], [h-c for h, c in zip(hi, p0)], 0)
    if verbose: print(f"  adaptive: ok={ok} leaves={stats['leaves']} beta={float(beta):.3f} {time.time()-t0:.1f}s", flush=True)
    return ok, stats, float(beta)
