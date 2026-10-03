"""Verified lower bound sigma2 <= lambda_{2,h}(p) on a box, for general charts (gcert.Problem): sigma_box.
 Strain transfer (LEMMAS 1.10, per piece q): v^T G_q(p) v >= mu_q v^T G_q(pc) v and J_q(p) <= j_q J_q(pc) on the box
 (piece_factors, Taylor-model ranges in ball arithmetic, directed conversions), so K(p) >= sum_q mu_q K_q(pc) and
 M(p) <= sum_q j_q M_q(pc) (Loewner); ONE verified LDL^T of the minorant A_low = sum mu_q K_q(pc) - s' sum j_q M_q(pc)
 (ldl.py, LEMMAS 1.9) with exactly one negative pivot gives lambda_2(K(p), M(p)) >= s' - (eA + dA)/lambda_min(M_low).
 Float assembly: K_low, M_low assembled from the stored reference matrices with float coefficients; the assembly error is bounded
 componentwise by rel * (sum |g||Khat| + s' sum J Mhat), rel = 14u + max relative error of the coefficients (ball arithmetic),
 and the float evaluation of the error terms is covered by the factor 1.01 (each term is a sum of at most a few thousand positive
 products: relative rounding error << 1e-10)."""
import numpy as np, scipy.sparse as sp
from flint import arb
from nanguard import fin, afin, fdown, fup, amax
from ldl import ldl_band
from gcert import pvars, poly_to_tm
from poly import Poly
from tm import TM

U = 2.0**-53

def piece_factors(gc, pc, p0, rad):
    """per piece q: (mu_q lower bound, j_q upper bound) over the box with
       v^T G_q(p) v >= mu_q v^T G_q(pc) v,   J_q(p) <= j_q J_q(pc)."""
    pb = gc.pb
    Z = pvars(np.asarray(p0, float), np.asarray(rad, float))   # RIGOR: exact: float64 arrays of Python floats (no rounding)
    Zc = [Poly.const(arb(float(v))) for v in pc]           # RIGOR: exact: float() of float64 box values is exact
    geo = pb.geometry(pb.crc, Z); geoc = pb.geometry(pb.crc, Zc)
    out = []
    for (g, gcn) in zip(geo, geoc):
        E = g[5]; Ec = [list(e.d.values())[0] if e.d else arb(0) for e in gcn[5]]
        e11, e12, e21, e22 = E; c11, c12, c21, c22 = Ec
        detT = poly_to_tm(e11*e22 - e12*e21)
        lo, hi = detT.range()
        if lo > 0: sgn = 1
        elif hi < 0: sgn = -1
        else: out.append((0.0, np.inf)); continue
        idet = (detT*sgn).recip()
        a11 = poly_to_tm(e22); a12 = poly_to_tm(-e12); a21 = poly_to_tm(-e21); a22 = poly_to_tm(e11)
        T11 = (a11*c11 + a21*c12)*idet*sgn; T12 = (a12*c11 + a22*c12)*idet*sgn
        T21 = (a11*c21 + a21*c22)*idet*sgn; T22 = (a12*c21 + a22*c22)*idet*sgn
        S11 = T11*T11 + T12*T12 - 1; S12 = T11*T21 + T12*T22; S22 = T21*T21 + T22*T22 - 1
        r11 = S11.range(); r12 = S12.range(); r22 = S22.range()
        apc_lo = (r11[0] + r22[0])/2
        # balls only (no nearest-float conversion): amc >= max |a - c| over a in r11, c in r22; b >= max |r12|
        amc = amax(abs(afin(r11[0] - r22[1])), abs(afin(r11[1] - r22[0])))/2
        bb = amax(abs(afin(r12[0])), abs(afin(r12[1])))
        lmin = afin(1 + apc_lo - (amc*amc + bb*bb).sqrt())
        dc = c11*c22 - c12*c21; dc = dc if dc > 0 else -dc
        jr = (detT*sgn)*(1/dc); jlo, jhi = jr.range()
        mu = fdown(arb(lmin.lower())*arb(jlo.lower()))       # directed: lower bound rounded down
        out.append((mu if mu > 0.0 else 0.0, fup(jhi)))      # upper bound rounded up; mu <= 0 is rejected by sigma_box
    return out

def sigma_box(gc, sol, p0, rad, frac=0.85):
    """sigma2 <= lambda_{2,h}(p) for all p in the box, via ONE verified LDL^T of the Loewner minorant
    A_low = sum_q mu_q K_q(pc) - s sum_q j_q M_q(pc)  (pc = p0):  K(p) >= sum mu_q K_q(pc), M(p) <= sum j_q M_q(pc)."""
    pb = gc.pb; pc = np.asarray(p0, float)   # RIGOR: exact: float64 arrays of Python floats (no rounding)
    fac = piece_factors(gc, pc, p0, rad)
    if any(m <= 0 for m, j in fac): return None, dict(why='piece factor <= 0')
    lam1, lam2 = sol['lamh'][0], sol['lamh'][1]
    Z = [Poly.const(arb(float(v))) for v in pc]            # RIGOR: exact: float() of float64 box values is exact
    geo = pb.geometry(pb.crc, Z)
    n = pb.cr['K'][0][0].shape[0]
    K = sp.csr_matrix((n, n)); absK = sp.csr_matrix((n, n)); md = np.zeros(n); rel = 12*U
    for q, (N11, N12, N22, d, adR, E) in enumerate(geo):
        mu, jq = fac[q]
        d0 = list(d.d.values())[0]; dv = d0 if d0 > 0 else -d0
        gs = [list(N.d.values())[0]/dv if N.d else arb(0) for N in (N11, N12, N22)]
        J = dv/adR
        for kl, g in enumerate(gs):
            gm = float(g.mid()); gf = gm*mu                  # RIGOR: accounted: |g - gm| <= |g.mid - gm| + g.rad enters rel below
            if gf != 0:
                K = K + gf*pb.cr['K'][q][kl]; absK = absK + abs(gf)*pb.cr['absK'][q][kl]
                rel = max(rel, 14*U + fup((abs(g.mid() - arb(gm)) + g.rad())/abs(arb(gm))))   # RIGOR: exact: max of floats; the .mid() term is exact ball arithmetic of the representation error
        Jm = float(J.mid()); Jf = Jm*jq; md = md + Jf*pb.cr['M'][q]; rel = max(rel, 14*U + fup((abs(J.mid() - arb(Jm)) + J.rad())/abs(arb(Jm))))   # RIGOR: accounted: the representation error of Jm (and of gm above) is bounded in ball arithmetic and enters rel; max of floats is exact
    if not (np.all(np.isfinite(md)) and float(np.min(md)) > 0): return None, dict(why='mass diagonal')   # RIGOR: exact: min of a float array
    for f in (frac, 0.6, 0.35):
        s_prime = lam1 + f*(lam2 - lam1)
        A = (K - s_prime*sp.diags(md)).tocsr()
        eA = rel*(float(np.max(np.asarray(absK.sum(axis=1)).ravel())) + s_prime*float(np.max(np.abs(md))))*1.01   # RIGOR: accounted: float row sums of positive terms, covered by the factor 1.01 (module docstring)
        r = ldl_band(A)
        if r['zero'] or r['neg'] != 1: continue
        X = (eA + r['dA_bound'])/float(np.min(md))*1.01      # RIGOR: accounted: upper bound of (eA + dA)/lambda_min(M); two float ops covered by 1.01
        s = fdown(arb(s_prime) - arb(X))                     # directed: sigma2 <= s' - X
        return s, dict(frac=f, n=r['n'], b=r['b'])
    return None, dict(why='inertia')
