"""RIGOROUS box certificate for xi = lambda2/lambda1 < 5/2 on quadrilaterals in the BD chart.

Chart: B=(0,0), D=(1,0), A=(ax,ay) (ay>0), C=(cx,cy) (cy<0); Q(p) = int(ABD u CBD u BD).  Box P = p0 + r*[-1,1]^4,
p = (ax, ay, cx, cy), t = (p - p0)/r.  For every p in P (floats enter only with explicit rigorous error bounds):
 (1) CR pencil (K(p), M(p)) on the affine images of a fixed exact reference mesh:  K = sum_s g_s(p) Khat_s,
     M = J_A Mhat_A + J_C Mhat_C  (refexact.py).  Test vector x(t) = x0 + sum_i t_i x_i (FIXED float vectors: discrete
     eigenvector at p0 and r_i times its p_i-derivative; any fixed vectors are admissible).
 (2) Temple:  lambda_{1,h}(p) >= rho - eps^2/(sigma2 - rho),  rho = RQ(x(t)),  eps^2 = ||K x - rho M x||^2_{M^-1}/||x||_M^2,
     sigma2 <= lambda_{2,h}(p) on P (sigma.py: verified LDL^T inertia + strain transfer).
 (3) Liu (2015), interpolation constant 0.1893 (Liu-Kikuchi):  lambda_1(Q(p)) >= lam/(1 + (0.1893 h)^2 lam),  h = max
     element diameter over P.
 (4) Rayleigh-Ritz upper bound for lambda_2(Q(p)): conforming P2 vectors y_a(t) = y_a + sum_i t_i y_{a,i}, a = 1,2:
     N(p) = K2 - (5/2) L1 M2 negative definite (N11 < 0, det N > 0)  =>  xi(Q(p)) < 5/2.
Scalar p-dependence: order-2 Taylor models (tm.py); residual: Taylor coefficient vectors with rigorous error vectors."""
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh, splu
from flint import arb
from tm import TM, NV, MON, NC
from assemble import Chart
from poly import Poly

U = 2.0**-53
def gamma(m): return m*U/(1 - m*U)
LIU = arb('0.1893')
SAFE = 1.01

def ball(v, r):
    return arb(float(v), float(r)*(1 + 1e-12) + 1e-300)

def rquad(x, K, y):
    """enclosure of x^T Kex y; K CSR float = correctly rounded exact matrix (|Kex - K| <= 1.0000001 u |K|)"""
    z = K @ y; s = float(x @ z)
    absb = float(np.abs(x) @ (abs(K) @ np.abs(y)))
    m = int(np.max(np.diff(K.indptr))) if K.nnz else 1
    return ball(s, (gamma(len(x)) + gamma(m) + 2.0000001*U)*absb*SAFE)

def rquad_diag(x, d, y):
    s = float(np.sum(x*d*y)); absb = float(np.sum(np.abs(x)*np.abs(d)*np.abs(y)))
    return ball(s, (gamma(len(x)) + 3*U)*absb*SAFE)

def geometry_tm(p0, r):
    if not (p0[1] - r[1] > 0 and p0[3] + r[3] < 0): raise ValueError('box crosses BD')
    ax = TM.var(0, p0[0], r[0]); ay = TM.var(1, p0[1], r[1]); cx = TM.var(2, p0[2], r[2]); cy = TM.var(3, p0[3], r[3])
    gs = []; Js = []
    for x, eta in ((ax, ay), (cx, -cy)):
        ie = eta.recip(); gs += [eta + x*x*ie, -(x*ie), ie]; Js.append(eta)
    return gs, Js          # gs: 6 TMs (A11,A12,A22,C11,C12,C22), Js: [J_A, J_C]

def geometry_derivs(p0):
    """first and second derivatives of g_s (6) and J_X (2) w.r.t. p at p0 (floats; used only to build test vectors)"""
    dg = np.zeros((4, 6)); dJ = np.zeros((4, 2)); d2g = np.zeros((4, 4, 6)); d2J = np.zeros((4, 4, 2))
    for X, (ix, iy, sgn) in enumerate(((0, 1, 1.0), (2, 3, -1.0))):
        x = p0[ix]; e = sgn*p0[iy]
        # g = (e + x^2/e, -x/e, 1/e);  de/dp_iy = sgn
        dg[ix, 3*X:3*X+3] = [2*x/e, -1/e, 0.0]
        dg[iy, 3*X:3*X+3] = sgn*np.array([1 - x*x/(e*e), x/(e*e), -1/(e*e)])
        dJ[iy, X] = sgn
        d2g[ix, ix, 3*X:3*X+3] = [2/e, 0.0, 0.0]
        d2g[ix, iy, 3*X:3*X+3] = sgn*np.array([-2*x/(e*e), 1/(e*e), 0.0]); d2g[iy, ix, 3*X:3*X+3] = d2g[ix, iy, 3*X:3*X+3]
        d2g[iy, iy, 3*X:3*X+3] = np.array([2*x*x/e**3, -2*x/e**3, 2/e**3])     # sgn^2 = 1
    return dg, dJ, d2g, d2J

def g_at(p0):
    out = []; Js = []
    for (ix, iy, sgn) in ((0, 1, 1.0), (2, 3, -1.0)):
        x = p0[ix]; e = sgn*p0[iy]; out += [e + x*x/e, -x/e, 1/e]; Js.append(e)
    return np.array(out), np.array(Js)

def eig_derivs(Kt, Mt, gs, Js, dg, dJ, d2g, d2J, lam, x, r, order=2):
    """Scaled Taylor vectors of the eigenvector x(p) of K(p) x = lam(p) M(p) x (normalised x^T M x = 1):
    returns [r_i x_i (i=0..3)] + [r_i r_j x_ij (i<j), r_i^2 x_ii / 2 (i=j)] in the order of tm.MON degree-2 monomials.
    (Any vectors are admissible in the certificate; these are chosen to make the residual O(|r|^3).)"""
    n = len(x)
    def Mop(J): return sum(j*(sp.diags(Ms) if isinstance(Ms, np.ndarray) else Ms) for j, Ms in zip(J, Mt) if j != 0) if np.any(np.asarray(J) != 0) else sp.csr_matrix((n, n))
    def Kop(G): return sum(g*Ks for g, Ks in zip(G, Kt) if g != 0) if np.any(np.asarray(G) != 0) else sp.csr_matrix((n, n))
    K = Kop(gs).tocsc(); M = Mop(Js).tocsc()
    Mx = M @ x
    Bm = sp.bmat([[K - lam*M, sp.csc_matrix(Mx[:, None])], [sp.csc_matrix(Mx[None, :]), None]]).tocsc()
    lu = splu(Bm)
    Ki = [Kop(dg[i]) for i in range(4)]; Mi = [Mop(dJ[i]) for i in range(4)]
    xi = []; li = []
    for i in range(4):
        lam_i = float(x @ (Ki[i] @ x) - lam*(x @ (Mi[i] @ x)))
        rhs = -(Ki[i] @ x - lam_i*Mx - lam*(Mi[i] @ x))
        sol = lu.solve(np.concatenate([rhs, [-0.5*float(x @ (Mi[i] @ x))]]))
        xi.append(sol[:n]); li.append(lam_i)
    out = [r[i]*xi[i] for i in range(4)]
    if order < 2: return out
    for i in range(4):
        for j in range(i, 4):
            Kij = Kop(d2g[i, j]); Mij = Mop(d2J[i, j])
            # (K - lam M) x_ij = -[(K_i - l_i M - lam M_i) x_j + (K_j - l_j M - lam M_j) x_i
            #                    + (K_ij - l_ij M - l_i M_j - l_j M_i - lam M_ij) x]
            Ti = Ki[i] @ xi[j] - li[i]*(M @ xi[j]) - lam*(Mi[i] @ xi[j])
            Tj = Ki[j] @ xi[i] - li[j]*(M @ xi[i]) - lam*(Mi[j] @ xi[i])
            T0 = Kij @ x - li[i]*(Mi[j] @ x) - li[j]*(Mi[i] @ x) - lam*(Mij @ x)
            lij = float(x @ (Ti + Tj + T0))
            rhs = -(Ti + Tj + T0 - lij*Mx)
            # normalisation: x_ij^T M x = -(x_i^T M x_j + x_i^T M_j x + x_j^T M_i x + x^T M_ij x / 2)
            nrm = -(xi[i] @ (M @ xi[j]) + xi[i] @ (Mi[j] @ x) + xi[j] @ (Mi[i] @ x) + 0.5*(x @ (Mij @ x)))
            sol = lu.solve(np.concatenate([rhs, [nrm]]))
            fac = r[i]*r[j] if i != j else 0.5*r[i]*r[i]
            out.append(fac*sol[:n])
    return out

def hmax_box(ch, p0, r):
    ev = ch.R['evecs']; h = 0.0
    for (ix, iy) in ((0, 1), (2, 3)):
        xl, xh = p0[ix] - r[ix], p0[ix] + r[ix]; yh = abs(p0[iy]) + r[iy]
        u1 = np.maximum(np.abs(ev[:, 0] + xl*ev[:, 1]), np.abs(ev[:, 0] + xh*ev[:, 1]))
        u2 = yh*np.abs(ev[:, 1])
        h = max(h, float(np.max(np.sqrt(u1*u1 + u2*u2))))
    return h*(1 + 1e-12)

def monomial_tm(alpha):
    """alpha indexes tm.MON: 0 -> 1, 1..4 -> t_i, 5.. -> t_i t_j"""
    m = MON[alpha]
    out = TM.const(1)
    for i in m: out = out*TM.var(i, 0.0, 1.0)
    return out

def block_quad(X, K, Y, absK=None, mrow=None, diag=False):
    """enclosures of X^T Kex Y (X: n x a, Y: n x b) as an a x b list of arb; K = correctly rounded exact matrix"""
    n = X.shape[0]
    if diag:
        Z = K[:, None]*Y; S = X.T @ Z; Sabs = np.abs(X).T @ (np.abs(K)[:, None]*np.abs(Y)); g = gamma(n) + 3*U
    else:
        Z = K @ Y; S = X.T @ Z; Sabs = np.abs(X).T @ (absK @ np.abs(Y)); g = gamma(n) + gamma(mrow) + 2.0000001*U
    R = g*Sabs*SAFE
    return [[ball(S[i, j], R[i, j]) for j in range(S.shape[1])] for i in range(S.shape[0])]

def qform_tm(xs, ys, mats, coefs, diag=False, absmats=None, mrows=None):
    """TM of x(t)^T (sum_s coef_s(t) mats_s) y(t), x(t) = sum_alpha t_alpha xs[alpha] (t_0 = 1)"""
    X = np.array(xs).T; Y = np.array(ys).T
    blocks = [block_quad(X, M, Y, None if diag else absmats[i], None if diag else mrows[i], diag) for i, M in enumerate(mats)]
    tot = TM.const(0)
    for a in range(X.shape[1]):
        for b in range(Y.shape[1]):
            lin = TM.const(0)
            for cs, blk in zip(coefs, blocks): lin = lin + cs*blk[a][b]
            if a == 0 and b == 0: tot = tot + lin
            else: tot = tot + lin*(monomial_tm(a)*monomial_tm(b))
    return tot

class BoxCert:
    def __init__(self, ch, order=2):
        self.ch = ch; self.order = order
        f = ch.crfree
        self.Kt = [ch.cr[(s, nm)][f][:, f].tocsr() for s in 'AC' for nm in ('11', '12', '22')]
        self.Md = [ch.cr[('A', 'M')][f], ch.cr[('C', 'M')][f]]
        g = ch.p2free
        self.Pt = [ch.p2[(s, nm)][g][:, g].tocsr() for s in 'AC' for nm in ('11', '12', '22')]
        self.PM = [ch.p2[(s, 'M')][g][:, g].tocsr() for s in 'AC']
        self.part = ch.crpart
        self.absKt = [abs(K).tocsr() for K in self.Kt]
        self.mKt = [int(np.max(np.diff(K.indptr))) for K in self.Kt]
        self.absPt = [abs(K).tocsr() for K in self.Pt]
        self.mPt = [int(np.max(np.diff(K.indptr))) for K in self.Pt]
        self.absPM = [abs(K).tocsr() for K in self.PM]
        self.mPM = [int(np.max(np.diff(K.indptr))) for K in self.PM]

    def solve(self, p0, r):
        p0 = np.asarray(p0, float); r = np.asarray(r, float)
        gs, Js = g_at(p0); dg, dJ, d2g, d2J = geometry_derivs(p0)
        K = sum(g*Ks for g, Ks in zip(gs, self.Kt)).tocsr(); md = Js[0]*self.Md[0] + Js[1]*self.Md[1]
        w, V = eigsh(K, k=3, M=sp.diags(md), sigma=0.0, which='LM'); o = np.argsort(w); w = w[o]; V = V[:, o]
        x0 = V[:, 0]/np.sqrt(V[:, 0] @ (md*V[:, 0]))
        xs = [x0] + eig_derivs(self.Kt, self.Md, gs, Js, dg, dJ, d2g, d2J, w[0], x0, r, order=self.order)
        K2 = sum(g*Ks for g, Ks in zip(gs, self.Pt)).tocsr(); M2 = (Js[0]*self.PM[0] + Js[1]*self.PM[1]).tocsr()
        w2, V2 = eigsh(K2, k=3, M=M2, sigma=0.0, which='LM'); o = np.argsort(w2); w2 = w2[o]; V2 = V2[:, o]
        ys = []
        for a in range(2):
            y = V2[:, a]/np.sqrt(V2[:, a] @ (M2 @ V2[:, a]))
            # near-degenerate eigenvalues (gap to the neighbours small) make the derivative vectors huge; any fixed vectors are
            # admissible, so drop corrections whose M-norm exceeds 0.25 (a vector in a near-degenerate eigenspace is already good)
            gap = min([abs(w2[a] - w2[b]) for b in range(len(w2)) if b != a])
            if gap < 0.02*w2[a]:
                ys.append([y]); continue
            ders = eig_derivs(self.Pt, self.PM, gs, Js, dg, dJ, d2g, d2J, w2[a], y, r, order=1)
            nrm = max(np.sqrt(abs(d @ (M2 @ d))) for d in ders)
            ys.append([y] + ders if nrm < 0.25 else [y])
        return dict(xs=xs, lamh=w, ys=ys, lamp2=w2, p0=p0, r=r)

    def certify(self, p0, r, sigma2, sol=None):
        p0 = np.asarray(p0, float); r = np.asarray(r, float)
        if sol is None: sol = self.solve(p0, r)
        xs = sol['xs']; n = len(xs[0])
        gt, Jt = geometry_tm(p0, r)
        # ---- (1) Rayleigh quotient ----
        # rho = lam0 + sum_X (1/eta_X) [ x^T P_X(t) x - lam0 eta_X^2 x^T Mhat_X x ] / b,   P_X = sum_s pX_s Khat_X,s
        lam0 = float(sol['lamh'][0])
        Xv = np.array(xs).T
        axP = Poly.var(0, p0[0], r[0]); ayP = Poly.var(1, p0[1], r[1]); cxP = Poly.var(2, p0[2], r[2]); cyP = Poly.var(3, p0[3], r[3])
        eA = ayP; eC = -cyP
        pAP = [eA*eA + axP*axP, -axP, Poly.const(1)]; pCP = [eC*eC + cxP*cxP, -cxP, Poly.const(1)]
        monoP = [Poly.mono(al) for al in range(len(xs))]
        def vecquad_poly(blocks, coefs):
            out = Poly.const(0)
            for a_ in range(len(xs)):
                for b_ in range(len(xs)):
                    lin = Poly.const(0)
                    for cs, blk in zip(coefs, blocks): lin = lin + cs*blk[a_][b_]
                    out = out + lin*(monoP[a_]*monoP[b_])
            return out
        bK = [block_quad(Xv, self.Kt[i], Xv, self.absKt[i], self.mKt[i]) for i in range(6)]
        bM = [block_quad(Xv, self.Md[X], Xv, diag=True) for X in range(2)]
        bracketA = vecquad_poly(bK[0:3], pAP) - vecquad_poly([bM[0]], [eA*eA*lam0])
        bracketC = vecquad_poly(bK[3:6], pCP) - vecquad_poly([bM[1]], [eC*eC*lam0])
        bpoly = vecquad_poly(bM, [eA, eC])
        def poly_to_tm(P):
            T = TM.const(0); c = [arb(0)]*NC; R = arb(0)
            for e, v in P.d.items():
                deg = sum(e)
                if deg <= 2:
                    m = tuple(sorted(sum([[i]*e[i] for i in range(NV)], [])))
                    from tm import IDX
                    c[IDX[m]] = c[IDX[m]] + v
                else: R = R + abs(v)
            return TM(c, R)
        b = poly_to_tm(bpoly)
        ib = b.recip()
        ieA = gt[2]; ieC = gt[5]      # 1/eta_A, 1/eta_C TMs
        rho = TM.const(lam0) + (ieA*poly_to_tm(bracketA) + ieC*poly_to_tm(bracketC))*ib
        # ---- (2) residual, scaled per dof set so that all stiffness coefficients are exact polynomials ----
        # G_X = (1/eta_X) * (eta_X^2 + x_X^2, -x_X, 1) =: pX/eta_X ;  mass M = eta_A Mhat_A + eta_C Mhat_C.
        rho_p, rho_R = Poly.from_tm(rho)
        ax = Poly.var(0, p0[0], r[0]); ay = Poly.var(1, p0[1], r[1]); cx = Poly.var(2, p0[2], r[2]); cy = Poly.var(3, p0[3], r[3])
        etaA = ay; etaC = -cy
        pA = [etaA*etaA + ax*ax, -ax, Poly.const(1)]; pC = [etaC*etaC + cx*cx, -cx, Poly.const(1)]
        etaA_hi = float((arb(p0[1]) + arb(r[1])).upper()); etaC_hi = float((arb(-p0[3]) + arb(r[3])).upper())
        etaA_lo = float((arb(p0[1]) - arb(r[1])).lower()); etaC_lo = float((arb(-p0[3]) - arb(r[3])).lower())
        nal = len(xs)
        colsK = [[self.Kt[s] @ xs[al] for s in range(6)] for al in range(nal)]
        colsM = [[self.Md[X]*xs[al] for X in range(2)] for al in range(nal)]
        errK = []
        for al in range(nal):
            axa = np.abs(xs[al]); row = []
            for s in range(6):
                row.append((gamma(self.mKt[s]) + 2.0000001*U)*(self.absKt[s] @ axa)*SAFE)
            errK.append(row)
        errM = [[3*U*np.abs(self.Md[X])*np.abs(xs[al])*SAFE for X in range(2)] for al in range(nal)]
        mono = [Poly.mono(al) for al in range(nal)]
        def set_norm(idx, entries, w):
            """entries: list of (coef Poly, coef remainder (arb), column (full vec), column error (full vec));
            returns rigorous upper bound of sup_t ||sum coef(t) col||_w over rows idx"""
            keys = sorted({k for (cp, cr, c, e) in entries for k in cp.d})
            kidx = {k: i for i, k in enumerate(keys)}
            ncol = len(entries)
            Cm = np.empty((len(idx), ncol)); Em = np.empty((len(idx), ncol))
            cmid = np.zeros((ncol, len(keys))); crad = np.zeros((ncol, len(keys))); crem = np.zeros(ncol); cabs = np.zeros(ncol)
            for j, (cp, cr, c, e) in enumerate(entries):
                Cm[:, j] = c[idx]; Em[:, j] = e[idx]
                for kk, v in cp.d.items():
                    cmid[j, kidx[kk]] = float(v.mid()); crad[j, kidx[kk]] = float(v.rad())
                crem[j] = float(arb(cr).upper()); cabs[j] = float(cp.absbound().upper())
            Rg = Cm @ cmid
            aC = np.abs(Cm)
            err = (aC @ (crad + gamma(ncol + 1)*np.abs(cmid))).sum(axis=1) + Em @ (cabs) + (aC + Em) @ crem
            err = err*SAFE
            nr = np.sqrt(((Rg**2)*w[:, None]).sum(axis=0))
            ne = np.sqrt(np.sum(err**2*w))
            return (float(np.sum(nr)) + float(ne))*(1 + 1e-9)
        tot = 0.0
        # set A rows: r = (1/etaA) [ sum_s pA_s K_s x - rho etaA^2 MA x ]
        idx = self.part['A']
        if len(idx):
            ent = []
            for al in range(nal):
                for s in range(3): ent.append((pA[s]*mono[al], arb(0), colsK[al][s], errK[al][s]))
                q = rho_p*etaA*etaA*mono[al]; ent.append((-q, rho_R*etaA_hi**2, colsM[al][0], errM[al][0]))
            nrm = set_norm(idx, ent, (1 + 8*U)/self.Md[0][idx])
            tot += nrm**2/etaA_lo**3
        idx = self.part['C']
        if len(idx):
            ent = []
            for al in range(nal):
                for s in range(3): ent.append((pC[s]*mono[al], arb(0), colsK[al][3+s], errK[al][3+s]))
                q = rho_p*etaC*etaC*mono[al]; ent.append((-q, rho_R*etaC_hi**2, colsM[al][1], errM[al][1]))
            nrm = set_norm(idx, ent, (1 + 8*U)/self.Md[1][idx])
            tot += nrm**2/etaC_lo**3
        idx = self.part['BD']
        if len(idx):
            ent = []
            for al in range(nal):
                for s in range(3): ent.append((etaC*pA[s]*mono[al], arb(0), colsK[al][s], errK[al][s]))
                for s in range(3): ent.append((etaA*pC[s]*mono[al], arb(0), colsK[al][3+s], errK[al][3+s]))
                qA = rho_p*etaA*etaC*etaA*mono[al]; ent.append((-qA, rho_R*etaA_hi**2*etaC_hi, colsM[al][0], errM[al][0]))
                qC = rho_p*etaA*etaC*etaC*mono[al]; ent.append((-qC, rho_R*etaA_hi*etaC_hi**2, colsM[al][1], errM[al][1]))
            J0A = float(p0[1]); J0C = float(-p0[3])
            mu = min(etaA_lo/J0A, etaC_lo/J0C)*(1 - 1e-12)
            w = (1 + 8*U)/(J0A*self.Md[0][idx] + J0C*self.Md[1][idx])
            nrm = set_norm(idx, ent, w)
            tot += nrm**2/(mu*(etaA_lo*etaC_lo)**2)
        eps2b = arb(tot*(1 + 1e-9))
        eps2 = ib*eps2b
        sig = arb(sigma2)
        gap = TM.const(sig) - rho
        if not (gap.lower() > 0): return dict(ok=False, why='sigma2 <= rho')
        templ = eps2*gap.recip()
        lam = rho - templ
        # ---- (3) Liu ----
        h = hmax_box(self.ch, p0, r)
        C2 = (LIU*arb(h))**2
        L1 = lam*(TM.const(1) + lam*C2).recip()
        # ---- (4) P2 Rayleigh-Ritz with corrected vectors ----
        ys = sol['ys']
        K2 = [[qform_tm(ys[i], ys[j], self.Pt, gt, absmats=self.absPt, mrows=self.mPt) for j in range(2)] for i in range(2)]
        M2 = [[qform_tm(ys[i], ys[j], self.PM, Jt, absmats=self.absPM, mrows=self.mPM) for j in range(2)] for i in range(2)]
        detM = M2[0][0]*M2[1][1] - M2[0][1]*M2[1][0]
        if not (detM.lower() > 0 and M2[0][0].lower() > 0): return dict(ok=False, why='M2 not PD')
        s = L1*2.5
        N = [[K2[i][j] - s*M2[i][j] for j in range(2)] for i in range(2)]
        # N negative definite  <=>  N11 < 0 and Schur complement N22 - N12^2/N11 < 0;
        # since N11 < 0:  N22 - N12^2/N11 <= N22 + N12^2/|N11|  (upper bounds of each term over the box)
        n11 = N[0][0].upper()
        n12sq = (N[0][1]*N[0][1]).upper() if False else arb(max(abs(float(v.upper())) for v in (N[0][1].range()[0], N[0][1].range()[1])))**2
        if not (n11 < 0): return dict(ok=False, why='N11 >= 0', n11=float(n11.mid()))
        schur = N[1][1].upper() + n12sq/abs(n11)
        detN = -schur
        ok = bool(schur < 0)
        rr = lambda T: tuple(float(v.mid()) for v in T.range())
        return dict(ok=ok, n11=float(n11.mid()), detN=float(detN.mid()), h=h, rho=rr(rho), L1=rr(L1),
                    temple=float(templ.upper().mid()), lamh=float(sol['lamh'][0]), lam2p2=float(sol['lamp2'][1]),
                    xi_float=float(sol['lamp2'][1]/sol['lamh'][0]))
