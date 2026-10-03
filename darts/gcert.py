"""RIGOROUS box certificate on general multi-piece affine charts (chart.py).

Problem data for a box P = p0 + r*[-1,1]^4 (t = (p - p0)/r):
  vfun(Z) -> dict label -> (x, y): vertex positions, POLYNOMIAL in the parameters Z = (Z0..Z3); Z may be floats or poly.Poly.
  cr_chart: Crouzeix-Raviart lower bound for lambda_1 (mixed BCs allowed: Neumann cuts, Lemmas 6.1, 9.1, 9.3)
  p2_chart: conforming P2 Rayleigh-Ritz upper bound for lambda_2 of a domain contained in Q (Dirichlet everywhere)
  tails: list of constants c_i with lambda_1(Q) >= min(mu_1(cr domain), c_i)  (Neumann bracketing tails), or empty.
Per piece q: E_q = [P1-P0, P2-P0], R_q = reference edge matrix; G_q = N_q/d_q with N_q = R adj(E)adj(E)^T R^T/|det R|,
d_q = |det E_q| (sign fixed on the box, checked), J_q = d_q/|det R_q|; K = sum_q (1/d_q) sum_kl N_q,kl Khat_q^kl, M = sum_q J_q Mhat_q.
Certificate (all rigorous, see LEMMAS.md section 8):  Temple + Liu for lambda_1, Rayleigh-Ritz for lambda_2, Schur test of
N = K2 - (5/2) min(L1, tails) M2 < 0 on P."""
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh, splu
from flint import arb
from nanguard import fin, afin, fdown, fup, amax
from tm import TM, NV, MON, NC, IDX, ORDER
from poly import Poly
from ldl import ldl_band
from schur import schur_upper

U = 2.0**-53
def gamma(m): return m*U/(1 - m*U)
LIU = arb('0.1893')
SAFE = 1.01

def ball(v, r): return arb(float(v), float(r)*(1 + 1e-12) + 1e-300)   # RIGOR: exact: float() of float64 is exact; the radius is only inflated (r (1 + 1e-12) + 1e-300 >= r after rounding)

# ---------------------------------------------------------------- polynomial helpers
def pvars(p0, r):
    return [Poly.var(i, p0[i], r[i]) for i in range(4)]   # RIGOR: exact: Poly.var builds the exact variable p0 + r t (not a NumPy variance)
def poly_eval_float(P, t=None):   # RIGOR-FUNCTION: heuristic: float evaluation for the float solve and for free positive normalisers (J0), never a bound
    """value at t (default 0) using float mids"""
    if t is None: return float(sum(float(v.mid()) for e, v in P.d.items() if sum(e) == 0))
    s = 0.0
    for e, v in P.d.items():
        s += float(v.mid())*np.prod([t[i]**e[i] for i in range(NV)])
    return s
def poly_coeff_float(P, e):   # RIGOR-FUNCTION: heuristic: float Taylor coefficients for the float solve (test vectors) only
    v = P.d.get(tuple(e)); return float(v.mid()) if v is not None else 0.0
def poly_to_tm(P):
    c = [arb(0)]*NC; R = arb(0)
    for e, v in P.d.items():
        deg = sum(e)
        if deg <= ORDER:
            m = tuple(sorted(sum([[i]*e[i] for i in range(NV)], []))); c[IDX[m]] = c[IDX[m]] + v
        else: R = R + abs(v)
    return TM(c, R)
def poly_range_abs(P):
    s = arb(0)
    for v in P.d.values(): s = s + abs(v)
    return s

class Problem:
    def __init__(self, cr_chart, p2_chart, vfun, tails_fn=None, name=''):
        self.crc = cr_chart; self.p2c = p2_chart; self.vfun = vfun; self.tails_fn = tails_fn; self.name = name
        self.cr = self._mats(cr_chart, 'cr'); self.p2 = self._mats(p2_chart, 'p2')
    def _mats(self, ch, fam):
        out = dict(K=[], M=[], absK=[], mrow=[], ch=ch)
        if fam == 'cr':
            f = ch.crfree; n = ch.nedge
            for q in range(ch.npiece):
                ks = []
                for nm in ('11', '12', '22'):
                    r_, c_, v_ = ch.cr[(q, nm)]; ks.append(sp.csr_matrix((v_, (r_, c_)), shape=(n, n))[f][:, f].tocsr())
                out['K'].append(ks); out['M'].append(ch.cr[(q, 'M')][f])
            out['free'] = f
            # dof groups by piece sets
            groups = {}
            for jj, j in enumerate(f): groups.setdefault(ch.edge_pieces[j], []).append(jj)
            out['groups'] = {k: np.array(v) for k, v in groups.items()}   # RIGOR: exact: integer DOF indices
        else:
            f = ch.p2free; n = ch.np2
            for q in range(ch.npiece):
                ks = []
                for nm in ('11', '12', '22'):
                    r_, c_, v_ = ch.p2[(q, nm)]; ks.append(sp.csr_matrix((v_, (r_, c_)), shape=(n, n))[f][:, f].tocsr())
                out['K'].append(ks)
                r_, c_, v_ = ch.p2[(q, 'M')]; out['M'].append(sp.csr_matrix((v_, (r_, c_)), shape=(n, n))[f][:, f].tocsr())
            out['free'] = f
        out['absK'] = [[abs(K).tocsr() for K in ks] for ks in out['K']]
        out['mrow'] = [[int(np.max(np.diff(K.indptr))) if K.nnz else 1 for K in ks] for ks in out['K']]   # RIGOR: exact: integer row lengths
        if fam == 'p2':
            out['absM'] = [abs(M).tocsr() for M in out['M']]; out['mrowM'] = [int(np.max(np.diff(M.indptr))) for M in out['M']]   # RIGOR: exact: integer row lengths
        return out
    # ------------------------------------------------------------ geometry as polynomials
    def geometry(self, ch, Z):
        """per piece: (N11, N12, N22, d, J) as polynomials/floats in Z; d made positive"""
        V = self.vfun(Z); out = []
        for q, pc in enumerate(ch_pieces(ch)):
            L0, L1, L2 = pc['labels']; P0, P1, P2 = V[L0], V[L1], V[L2]
            e11 = P1[0] - P0[0]; e21 = P1[1] - P0[1]; e12 = P2[0] - P0[0]; e22 = P2[1] - P0[1]
            (r0, r1, r2) = ch.Rref_exact[q]
            R11 = r1[0] - r0[0]; R21 = r1[1] - r0[1]; R12 = r2[0] - r0[0]; R22 = r2[1] - r0[1]
            detR = R11*R22 - R12*R21
            # adj(E) = [[e22, -e12],[-e21, e11]];  W = R adj(E)
            Ra = [[_c(R11), _c(R12)], [_c(R21), _c(R22)]]; adRa = abs(_c(detR))
            a11 = e22; a12 = -e12; a21 = -e21; a22 = e11
            W11 = Ra[0][0]*a11 + Ra[0][1]*a21; W12 = Ra[0][0]*a12 + Ra[0][1]*a22
            W21 = Ra[1][0]*a11 + Ra[1][1]*a21; W22 = Ra[1][0]*a12 + Ra[1][1]*a22
            N11 = (W11*W11 + W12*W12)*(1/adRa); N12 = (W11*W21 + W12*W22)*(1/adRa); N22 = (W21*W21 + W22*W22)*(1/adRa)
            d = e11*e22 - e12*e21
            out.append([N11, N12, N22, d, adRa, (e11, e12, e21, e22)])
        return out

def _c(x):
    """exact Fraction -> arb (exact for dyadic), float stays float"""
    from fractions import Fraction
    if isinstance(x, Fraction): return arb(x.numerator)/arb(x.denominator)
    return x

def ch_pieces(ch): return ch.pieces_meta

# ------------------------------------------------------------------ main certificate
class GCert:
    def __init__(self, prob, order=2):
        self.pb = prob; self.order = order
    # ---- numerics at p0 (floats) ----
    def _float_geo(self, ch, p0, r):   # RIGOR-FUNCTION: heuristic: float geometry for the float solve (test vectors) only
        """Taylor coefficients (in t) up to degree 2 of g_q,kl = N/d and J_q, as float arrays"""
        Z = pvars(p0, r)
        geo = self.pb.geometry(ch, Z)
        out = []
        for (N11, N12, N22, d, adR, E) in geo:
            d0 = poly_eval_float(d); sgn = 1.0 if d0 > 0 else -1.0
            dd = d*sgn
            # 1/d to second order: (1/d0)(1 - delta + delta^2), delta = (d - d0)/d0
            d0a = arb(abs(d0))
            delta = (dd - Poly.const(abs(d0)))*(1/d0a)
            inv = (Poly.const(1) - delta + delta*delta)*(1/d0a)
            gs = [N*inv for N in (N11, N12, N22)]
            J = dd*(1/adR)
            out.append((gs, J, sgn))
        return out
    def _taylor_mats(self, mats, fgeo, family):   # RIGOR-FUNCTION: heuristic: float matrices for the float solve (test vectors) only
        """K(t), M(t) Taylor coefficient matrices: dict mono -> (K, M) for monos of degree <= 2 (float)"""
        monos = [()] + [(i,) for i in range(4)] + [(i, j) for i in range(4) for j in range(i, 4)]
        def ex(m):
            e = [0]*4
            for i in m: e[i] += 1
            return tuple(e)
        out = {}
        n = mats['K'][0][0].shape[0]
        for m in monos:
            K = sp.csr_matrix((n, n)); M = sp.csr_matrix((n, n)) if family == 'p2' else np.zeros(n)
            for q, (gs, J, sgn) in enumerate(fgeo):
                for kl in range(3):
                    c = poly_coeff_float(gs[kl], ex(m))
                    if c != 0: K = K + c*mats['K'][q][kl]
                c = poly_coeff_float(J, ex(m))
                if c != 0: M = M + c*mats['M'][q]
            out[m] = (K.tocsr(), (sp.diags(M) if family == 'cr' else M).tocsr())
        return out
    def _vectors(self, TK, lam_idx, nev, order):   # RIGOR-FUNCTION: heuristic: approximate eigenvectors and their Taylor transport (free test functions)
        K0, M0 = TK[()]
        w, V = eigsh(K0, k=nev, M=M0, sigma=0.0, which='LM', v0=np.ones(K0.shape[0]))
        o = np.argsort(w); w = w[o]; V = V[:, o]
        res = []
        for a in lam_idx:
            x = V[:, a]/np.sqrt(V[:, a] @ (M0 @ V[:, a])); lam = w[a]
            gap = min([abs(w[a] - w[b]) for b in range(len(w)) if b != a])
            vecs = {(): x}
            if order >= 1 and gap > 0.02*abs(lam):
                Mx = M0 @ x
                Bm = sp.bmat([[(K0 - lam*M0).tocsc(), sp.csc_matrix(Mx[:, None])], [sp.csc_matrix(Mx[None, :]), None]]).tocsc()
                lu = splu(Bm); n = len(x)
                xi = {}; li = {}
                for i in range(4):
                    Ki, Mi = TK[(i,)]
                    li[i] = float(x @ (Ki @ x) - lam*(x @ (Mi @ x)))
                    rhs = -(Ki @ x - li[i]*Mx - lam*(Mi @ x))
                    sol = lu.solve(np.concatenate([rhs, [-0.5*float(x @ (Mi @ x))]])); xi[i] = sol[:n]
                ok = max(np.sqrt(abs(xi[i] @ (M0 @ xi[i]))) for i in range(4)) < (0.5 if order >= 2 else 0.25)
                if not ok: res.append((vecs, lam)); continue
                for i in range(4): vecs[(i,)] = xi[i]
                if order >= 2:
                    for i in range(4):
                        for j in range(i, 4):
                            Kij, Mij = TK[(i, j)]; Ki, Mi = TK[(i,)]; Kj, Mj = TK[(j,)]
                            Ai = lambda v: Ki @ v - li[i]*(M0 @ v) - lam*(Mi @ v)
                            Aj = lambda v: Kj @ v - li[j]*(M0 @ v) - lam*(Mj @ v)
                            if i != j:
                                T_ = Ai(xi[j]) + Aj(xi[i]) + Kij @ x - li[i]*(Mj @ x) - li[j]*(Mi @ x) - lam*(Mij @ x)
                                nrm = -(xi[i] @ (M0 @ xi[j]) + xi[i] @ (Mj @ x) + xi[j] @ (Mi @ x) + 0.5*(x @ (Mij @ x)))
                            else:
                                T_ = Ai(xi[i]) + Kij @ x - li[i]*(Mi @ x) - lam*(Mij @ x)
                                nrm = -(0.5*(xi[i] @ (M0 @ xi[i])) + xi[i] @ (Mi @ x) + 0.5*(x @ (Mij @ x)))
                            lij = float(x @ T_)
                            rhs = -(T_ - lij*Mx)
                            sol = lu.solve(np.concatenate([rhs, [nrm]])); vecs[(i, j)] = sol[:n]
            res.append((vecs, lam))
        return w, res
    def solve(self, p0, r):
        p0 = np.asarray(p0, float); r = np.asarray(r, float)   # RIGOR: exact: float64 arrays of Python floats (no rounding)
        fg_cr = self._float_geo(self.pb.crc, p0, r); fg_p2 = self._float_geo(self.pb.p2c, p0, r)
        TKc = self._taylor_mats(self.pb.cr, fg_cr, 'cr'); TKp = self._taylor_mats(self.pb.p2, fg_p2, 'p2')
        wc, rc = self._vectors(TKc, [0], 3, self.order)
        wp, rp = self._vectors(TKp, [0, 1], 3, 1)
        return dict(p0=p0, r=r, lamh=wc, xvecs=rc[0][0], lam0=rc[0][1], yvecs=[rp[0][0], rp[1][0]], lamp2=wp, TKc=TKc)
    # ---- rigorous part ----
    def certify(self, p0, r, sigma2, sol=None):
        p0 = np.asarray(p0, float); r = np.asarray(r, float)   # RIGOR: exact: float64 arrays of Python floats (no rounding)
        if sol is None: sol = self.solve(p0, r)
        Z = pvars(p0, r)
        pb = self.pb
        # exact geometry polynomials + sign check
        geoC = pb.geometry(pb.crc, Z); geoP = pb.geometry(pb.p2c, Z)
        def prep(geo):
            out = []
            for (N11, N12, N22, d, adR, E) in geo:
                dT = poly_to_tm(d); lo, hi = dT.range()
                if lo > 0: sgn = 1
                elif hi < 0: sgn = -1
                else: raise ValueError('piece degenerates on the box')
                dpos = d*sgn
                dTp = poly_to_tm(dpos); invd = dTp.recip()
                J = dpos*(1/adR); JT = poly_to_tm(J)
                out.append(dict(N=[N11, N12, N22], d=dpos, dT=dTp, invd=invd, J=J, JT=JT, adR=adR, E=[[poly_to_tm(E[0]), poly_to_tm(E[1])], [poly_to_tm(E[2]), poly_to_tm(E[3])]]))
            return out
        GC = prep(geoC); GP = prep(geoP)
        mats = pb.cr
        keys = list(sol['xvecs'].keys()); xs = [sol['xvecs'][k] for k in keys]
        monoP = [mono_poly(k) for k in keys]
        X = np.array(xs).T   # RIGOR: exact: the float test vectors (free choice)
        lam0 = float(sol['lam0'])   # RIGOR: exact: any float shift is valid: rho = lam0 + (xKx - lam0 xMx)/xMx is exact algebra with this float
        # (1) Rayleigh quotient, cancelled form
        def vq(K, absK, mrow, diag=False):
            return block_quad(X, K, X, absK, mrow, diag)
        bM = [vq(mats['M'][q], None, None, True) for q in range(len(GC))]
        bpoly = Poly.const(0); num = TM.const(0)
        for q, g in enumerate(GC):
            bK = [vq(mats['K'][q][kl], mats['absK'][q][kl], mats['mrow'][q][kl]) for kl in range(3)]
            Bq = vecquad(bK, g['N'], monoP) - vecquad([bM[q]], [g['d']*g['J']*lam0], monoP)
            num = num + g['invd']*poly_to_tm(Bq)
            bpoly = bpoly + vecquad([bM[q]], [g['J']], monoP)
        b = poly_to_tm(bpoly); ib = b.recip()
        rho = TM.const(lam0) + num*ib
        rho_p, rho_R = Poly.from_tm(rho)
        # (2) residual by dof groups
        n = X.shape[0]
        colsK = {}; errK = {}
        for q in range(len(GC)):
            for kl in range(3):
                Kq = mats['K'][q][kl]; aK = mats['absK'][q][kl]; m = mats['mrow'][q][kl]
                for ai, xa in enumerate(xs):
                    colsK[(q, kl, ai)] = Kq @ xa
                    errK[(q, kl, ai)] = (gamma(m) + 2.0000001*U)*(aK @ np.abs(xa))*SAFE
        colsM = {(q, ai): mats['M'][q]*xa for q in range(len(GC)) for ai, xa in enumerate(xs)}
        errM = {(q, ai): 3*U*np.abs(mats['M'][q])*np.abs(xa)*SAFE for q in range(len(GC)) for ai, xa in enumerate(xs)}
        tot = arb(0)
        for S, idx in mats['groups'].items():
            prodd = Poly.const(1)
            for q in S: prodd = prodd*GC[q]['d']
            ent = []
            for q in S:
                others = Poly.const(1)
                for q2 in S:
                    if q2 != q: others = others*GC[q2]['d']
                for ai in range(len(xs)):
                    for kl in range(3):
                        ent.append((others*GC[q]['N'][kl]*monoP[ai], arb(0), colsK[(q, kl, ai)], errK[(q, kl, ai)]))
                    cp = rho_p*prodd*GC[q]['J']*monoP[ai]
                    rem = rho_R*poly_range_abs(prodd*GC[q]['J'])
                    ent.append((-cp, rem, colsM[(q, ai)], errM[(q, ai)]))
            # weights
            # weights w >= 1/M(p)_ii (M(p) = sum_q J_q(p) Mhat_q lumped, positive; Mhat_q stored as nearest floats of the exact
            # rationals): at most 3|S| + 1 relative roundings of positive quantities, covered by the factor 1 + (3|S| + 8) u
            if len(S) == 1:
                q = S[0]; Jlo = fdown(GC[q]['JT'].lower())
                if not Jlo > 0: return dict(ok=False, why='J<=0')
                w = (1 + 10*U)/(Jlo*mats['M'][q][idx])
            else:
                J0 = {q: poly_eval_float(GC[q]['J']) for q in S}
                if not all(J0[q] > 0 for q in S): return dict(ok=False, why='J0<=0')
                mu = min(fdown(GC[q]['JT'].lower())/J0[q] for q in S)*(1 - 1e-12)      # RIGOR: accounted: positive quotients (checked), one rounding each, times (1 - 1e-12): mu <= J_q(p)/J0_q; min of floats is exact
                if not mu > 0: return dict(ok=False, why='mu<=0')
                w = (1 + (3*len(S) + 8)*U)/(mu*sum(J0[q]*mats['M'][q][idx] for q in S))
            if not (np.all(np.isfinite(w)) and np.all(w > 0)): return dict(ok=False, why='weights')   # RIGOR: exact: boolean tests
            nrm = set_norm(idx, ent, w)
            pd_lo = fdown(poly_to_tm(prodd).lower())
            if not pd_lo > 0: return dict(ok=False, why='prod d <= 0')
            an = arb(nrm); ap = arb(pd_lo)
            tot = tot + an*an/(ap*ap)                                   # ball arithmetic (no float rounding)
        eps2 = ib*arb(fup(tot))
        gap = TM.const(arb(sigma2)) - rho
        if not (gap.lower() > 0): return dict(ok=False, why='sigma2 <= rho', rho=tm_rng(rho))
        lam = rho - eps2*gap.recip()
        # (3) Liu
        h = hmax_box(pb.crc, GC, geoC, Z)
        lh = LIU*arb(h); C2 = lh*lh
        L1 = lam*(TM.const(1) + lam*C2).recip()
        # (4) P2 Rayleigh-Ritz
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
        detM = M2[0][0]*M2[1][1] - M2[0][1]*M2[1][0]
        if not (detM.lower() > 0 and M2[0][0].lower() > 0): return dict(ok=False, why='M2 not PD')
        out = dict(rho=tm_rng(rho), L1=tm_rng(L1), h=h, temple=fup((eps2*gap.recip()).upper()),
                   lamh=float(sol['lamh'][0]), lam2p2=float(sol['lamp2'][1]))   # RIGOR: info: float eigenvalue estimates, logged only
        sigs = [L1]
        tails = pb.tails_fn(p0, r) if pb.tails_fn else []
        for c in tails: sigs.append(TM.const(arb(c)))
        worsts = []
        for s_ in sigs:
            N = [[K2[i][j] - (s_*2.5)*M2[i][j] for j in range(2)] for i in range(2)]
            ok_, schur, why = schur_upper(N)
            if not ok_: out.update(ok=False, why=why, **({'schur': fup(schur)} if schur is not None else {})); return out
            worsts.append(schur)
        out.update(ok=True, schur=fup(amax(*worsts)), tails=list(tails))
        return out

def tm_rng(T):
    """DIRECTED float enclosure (lo, hi) of the range of a Taylor model: lo <= T <= hi on the box"""
    lo, hi = T.range()
    return (fdown(lo), fup(hi))

def mono_poly(key):
    e = [0]*NV
    for i in key: e[i] += 1
    return Poly({tuple(e): arb(1)})

def block_quad(X, K, Y, absK=None, mrow=None, diag=False):
    n = X.shape[0]
    if diag:
        Z = K[:, None]*Y; S = X.T @ Z; Sabs = np.abs(X).T @ (np.abs(K)[:, None]*np.abs(Y)); g = gamma(n) + 3*U
    else:
        Z = K @ Y; S = X.T @ Z; Sabs = np.abs(X).T @ (absK @ np.abs(Y)); g = gamma(n) + gamma(mrow) + 2.0000001*U
    R = g*Sabs*SAFE
    return [[ball(S[i, j], R[i, j]) for j in range(S.shape[1])] for i in range(S.shape[0])]

def vecquad(blocks, coefs, mono):
    out = Poly.const(0)
    for a in range(len(mono)):
        for b in range(len(mono)):
            lin = Poly.const(0)
            for cs, blk in zip(coefs, blocks): lin = lin + cs*blk[a][b]
            out = out + lin*(mono[a]*mono[b])
    return out

def vecquad2(blocks, coefs, mi, mj):
    out = Poly.const(0)
    for a in range(len(mi)):
        for b in range(len(mj)):
            lin = Poly.const(0)
            for cs, blk in zip(coefs, blocks): lin = lin + cs*blk[a][b]
            out = out + lin*(mi[a]*mj[b])
    return out

def set_norm(idx, entries, w):
    if len(idx) > 10**6: raise ValueError('set_norm: group too large for the 1e-9 float-summation allowance')
    keys = sorted({k for (cp, cr, c, e) in entries for k in cp.d})
    kidx = {k: i for i, k in enumerate(keys)}
    ncol = len(entries)
    Cm = np.empty((len(idx), ncol)); Em = np.empty((len(idx), ncol))
    cmid = np.zeros((ncol, len(keys))); crad = np.zeros((ncol, len(keys))); crem = np.zeros(ncol); cabs = np.zeros(ncol)
    for j, (cp, cr, c, e) in enumerate(entries):
        Cm[:, j] = c[idx]; Em[:, j] = e[idx]
        for kk, v in cp.d.items():
            cmid[j, kidx[kk]] = float(v.mid()); crad[j, kidx[kk]] = fup(v.rad())   # RIGOR: accounted: |v - cmid| <= crad + u|cmid|, the u term is in gamma(ncol + 1) below
        crem[j] = fup(arb(cr)); cabs[j] = fup(cp.absbound())
    Rg = Cm @ cmid
    aC = np.abs(Cm)
    err = (aC @ (crad + gamma(ncol + 1)*np.abs(cmid))).sum(axis=1) + Em @ cabs + (aC + Em) @ crem   # RIGOR: accounted: float sums of nonnegative terms covered by SAFE = 1.01
    err = err*SAFE
    nr = np.sqrt(((Rg*Rg)*w[:, None]).sum(axis=0))   # RIGOR: accounted: covered by the factor 1 + 1e-9 below (group size < 1e6)
    ne = np.sqrt(np.sum(err*err*w))   # RIGOR: accounted: covered by the factor 1 + 1e-9 below (group size < 1e6)
    return (float(np.sum(nr)) + float(ne))*(1 + 1e-9)   # RIGOR: accounted: float sums of < 1e6 positive terms and sqrt, relative error < 1e-9 (group-size guard above)

def hmax_box(ch, G, geo, Z):
    """sup over the box of the max element diameter: |DF e| for reference edges e, DF = E R^{-1} (interval entries).
    R^{-1} is formed EXACTLY (rationals -> arb balls); DF entries are enclosed by directed floats [a_lo, a_hi]; for a reference
    edge e (ev = componentwise nearest float of the exact edge, |ev - e| <= u |ev|) and every corner (a, b) of the DF row box,
    |a e0 + b e1| <= |fl(a ev0 + b ev1)| + 8u (|a ev0| + |b ev1|) (product, sum and representation roundings); the maximum of the
    convex function |a e0 + b e1| over the box is attained at a corner; the final sqrt / max in floats is covered by 1 + 1e-9."""
    from fractions import Fraction as Fr
    h = 0.0
    for q in range(len(G)):
        (r0, r1, r2) = ch.Rref_exact[q]
        R = [[Fr(r1[0]) - Fr(r0[0]), Fr(r2[0]) - Fr(r0[0])], [Fr(r1[1]) - Fr(r0[1]), Fr(r2[1]) - Fr(r0[1])]]
        det = R[0][0]*R[1][1] - R[0][1]*R[1][0]
        if det == 0: raise ValueError('degenerate reference triangle')
        Ri = [[R[1][1]/det, -R[0][1]/det], [-R[1][0]/det, R[0][0]/det]]
        Ria = [[arb(x.numerator)/arb(x.denominator) for x in row] for row in Ri]
        E = G[q]['E']; ev = ch.evecs[q]
        DF = [[None, None], [None, None]]
        for i in range(2):
            for j in range(2):
                lo, hi = (E[i][0]*Ria[0][j] + E[i][1]*Ria[1][j]).range()
                DF[i][j] = (fdown(lo), fup(hi))
        def comp(row):
            (a_lo, a_hi), (b_lo, b_hi) = row
            outs = []
            for a in (a_lo, a_hi):
                for b in (b_lo, b_hi):
                    pa = np.abs(a*ev[:, 0]); pb = np.abs(b*ev[:, 1])
                    outs.append(np.abs(a*ev[:, 0] + b*ev[:, 1]) + 8*U*(pa + pb))
            return np.maximum.reduce(outs)   # RIGOR: exact: elementwise maximum of floats (selection); the rounding of the entries is accounted above
        u1 = comp(DF[0]); u2 = comp(DF[1])
        hq = fin(np.max(np.sqrt(u1*u1 + u2*u2)))   # RIGOR: accounted: the float squares / sqrt are covered by the final factor 1 + 1e-9
        h = max(h, hq)                              # RIGOR: exact: max of floats
    return h*(1 + 1e-9)
