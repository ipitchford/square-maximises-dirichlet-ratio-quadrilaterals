"""Verified lower bound for the second CR eigenvalue lambda_{2,h}(p) on a box (cell):
 (a) at a cell centre pc: LDL^T inertia of A_fl = K_fl(pc) - s' M_fl(pc) has exactly one negative pivot;
     ||A_fl - A_ex||_2 <= assembly bound, ||dA||_2 <= LDL bound  =>  K_ex - s M_ex has <= 1 nonpositive eigenvalue for
     s = s' - (eA + dA)/m_min  =>  lambda_{2,h}(pc) > s.
 (b) strain transfer to p in the cell: for every CR vector v, a_p(v) >= min_X mu_X(p) a_pc(v), m_p(v) <= max_X (J_X(p)/J_X(pc)) m_pc(v),
     mu_X = lambda_min(G_X(pc)^{-1} G_X(p)) = (tau - sqrt(tau^2-4))/2, tau = tr(G_X(pc)^{-1}G_X(p))  (det G_X = 1)."""
import numpy as np, scipy.sparse as sp
from flint import arb
from ldl import ldl_band, gamma, U
from assemble import Gfun

def Garb(X):
    x, y = arb(X[0]), abs(arb(X[1]))
    return [y + x*x/y, -x/y, 1/y], y

def sigma_centre(ch, pc, s_prime):
    A = (pc[0], pc[1]); C = (pc[2], pc[3])
    K, md = ch.cr_mats(A, C)
    Afl = (K - s_prime*sp.diags(md)).tocsr()
    # assembly error bound (componentwise): coefficients g, J computed in float vs exact arb; reference rounding u
    gA, JA = Garb(A); gC, JC = Garb(C)
    f = ch.crfree
    absK = sum(abs(float(gv.mid()))*abs(ch.cr[(s, nm)][f][:, f]) for s, gg in (('A', gA), ('C', gC)) for gv, nm in zip(gg, ('11', '12', '22')))
    dg = max(float(abs(gv - float(gv.mid())).upper()) / max(abs(float(gv.mid())), 1e-300) for gv in gA + gC)
    dJ = max(float(abs(JA - float(JA.mid())).upper())/float(JA.mid()), float(abs(JC - float(JC.mid())).upper())/float(JC.mid()))
    rel = 12*U + dg + dJ
    eA = rel*(np.max(np.asarray(absK.sum(axis=1)).ravel()) + s_prime*np.max(np.abs(md)))*1.01
    r = ldl_band(Afl)
    if r['zero'] or r['neg'] != 1: return None, r
    s = s_prime - (eA + r['dA_bound'])/np.min(md)*1.01
    return s, r

def strain_factor(pc, p0, rad):
    """Rigorous lower bound, over the box p0 +- rad, of   min_X mu_X(p) / max_X (J_X(p)/J_X(pc)).
    For X with F = [[1, x],[0, eta]] (eta = |y|):  v^T G_p v / v^T G_c v = (eta/eta_c) |T^T w|^2/|w|^2 with T = F_c F_p^{-1}
    = [[1, b],[0, c]],  b = (x_c - x)/eta,  c = eta_c/eta;  so mu_X >= (eta/eta_c) * lambda_min(T T^T) and
    lambda_min(T^T T) = 1 + (a - sqrt(a^2 + 4 b^2))/2,  a = b^2 + c^2 - 1  (increasing in a, decreasing in b^2).
    (T T^T and T^T T have the same eigenvalues.)"""
    pc = [arb(v) for v in pc]
    worst_mu = None; worst_J = None
    for (ix, iy, sg) in ((0, 1, 1), (2, 3, -1)):
        x = arb(p0[ix] - rad[ix]).union(arb(p0[ix] + rad[ix]))
        eta = arb(sg*p0[iy] - rad[iy]).union(arb(sg*p0[iy] + rad[iy]))
        if not (eta.lower() > 0): return arb(0)
        xc = pc[ix]; etac = sg*pc[iy]
        b = (xc - x)/eta; c = etac/eta
        B2 = arb((b*b).upper()); a = b*b + c*c - 1
        a_lo = arb(a.lower())
        lmin = 1 + (a_lo - (a_lo*a_lo + 4*B2).sqrt())/2
        mu = arb(((eta/etac)*lmin).lower())
        jr = arb((eta/etac).upper())
        worst_mu = mu if worst_mu is None else (mu if mu.lower() < worst_mu.lower() else worst_mu)
        worst_J = jr if worst_J is None else (jr if jr.upper() > worst_J.upper() else worst_J)
    if not (worst_mu.lower() > 0): return arb(0)
    return arb((worst_mu/worst_J).lower())
