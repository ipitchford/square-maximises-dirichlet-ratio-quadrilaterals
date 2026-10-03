"""Fast FLOAT approximation of bo2.bo_family, used ONLY to choose slicing directions and arm variants (cover heuristic, no
rigour): the certificate is always recomputed with the rigorous bo2.bo_family for the chosen direction, and the replay recomputes
it from the recorded direction.  Mirrors bo_family: same nodes, same per-slab rule (sum of the A and C extents if the slab meets
the projection of BD, else the max, each capped by the per-dart chord bound), same far-arm cuts and tails; float arithmetic,
vectorised slab extents, and the ground state of -v'' + V by vectorised Sturm shooting over a grid of trial energies."""
# RIGOR-MODULE: heuristic: float screening of slicing directions / arm modes for the production driver; never used by a replay
import math, numpy as np
from bo1d import graded_nodes

def _e(b, side): return np.array([math.cos(b), side*math.sin(b)])

def _polar_pts(arm):
    pts = [(0.0, 0.0), (1.0, 0.0)]
    for r in (arm.r1, arm.r2):
        for b in (arm.b1, arm.b2): pts.append(tuple(r*_e(b, arm.side)))
    half = (arm.b2 - arm.b1)/2
    if not half < math.pi/2: return None
    pts.append(tuple(arm.r2/math.cos(half)*_e((arm.b1 + arm.b2)/2, arm.side)))
    return np.array(pts)

def _far_pts(arm, u, xc):
    pts = [(0.0, 0.0), (1.0, 0.0)]
    for b in (arm.b1, arm.b2):
        e = _e(b, arm.side); de = e @ u
        for P0 in (np.zeros(2), np.array([1.0, 0.0])):
            lam = (xc - P0 @ u)/de
            if lam < 0: return None
            pts.append(tuple(P0 + lam*e))
    return np.array(pts)

def _dir_sign(arm, u):
    ds = [_e(b, arm.side) @ u for b in np.linspace(arm.b1, arm.b2, 5)]
    if min(ds) > 0: return 1
    if max(ds) < 0: return -1
    return 0

def _chord(arm, u):
    bs = np.linspace(arm.b1, arm.b2, 9)
    d = np.array([_e(b, arm.side) @ u for b in bs])
    if d.min() <= 0 <= d.max(): return None
    return float(np.max(np.sin(bs)/np.abs(d)))

def _ext(P, u, nodes):
    """approximate y-extent of conv(P) in each slab [nodes[k], nodes[k+1]]"""
    up = np.array([-u[1], u[0]])
    X = P @ u; Y = P @ up
    i, j = np.triu_indices(len(P), 1)
    dx = X[j] - X[i]
    c = np.asarray(nodes)[:, None]
    with np.errstate(divide='ignore', invalid='ignore'):
        t = (c - X[i])/dx
        y = Y[i] + t*(Y[j] - Y[i])
    ok = (dx != 0) & (t >= 0) & (t <= 1)
    ylo = np.where(ok, y, np.inf).min(1); yhi = np.where(ok, y, -np.inf).max(1)
    lo = np.minimum(ylo[:-1], ylo[1:]); hi = np.maximum(yhi[:-1], yhi[1:])
    nd = np.asarray(nodes)
    for xp, yp in zip(X, Y):
        k = np.nonzero((nd[:-1] <= xp) & (xp <= nd[1:]))[0]
        lo[k] = np.minimum(lo[k], yp); hi[k] = np.maximum(hi[k], yp)
    e = hi - lo
    return np.where(np.isfinite(e) & (e > 0), e, 0.0)

def _mu(nodes, H, kL, kR):
    """approximate bottom of the spectrum of -v'' + pi^2/H^2 (Dirichlet ends; exponential tails approximated by appended
    regions of length 8/sqrt(k) with V = k): finite volumes on the nodes, lowest eigenvalue of a tridiagonal matrix"""
    from scipy.linalg import eigh_tridiagonal
    x0 = np.asarray(nodes, float); V0 = np.minimum((math.pi/np.maximum(H, 1e-6))**2, 1e12)
    # merge nodes closer than 1e-7 of the length (FD stiffness 1/h would blow up); a merged cell takes the max potential
    tol = 1e-7*(x0[-1] - x0[0]); x = [x0[0]]; V = []; vmax = 0.0
    for k in range(len(V0)):
        vmax = max(vmax, V0[k])
        if x0[k + 1] - x[-1] > tol or k == len(V0) - 1:
            if x0[k + 1] - x[-1] > tol: x.append(x0[k + 1]); V.append(vmax)
            else: V[-1] = max(V[-1], vmax) if V else vmax; x[-1] = x0[k + 1] if not V else x[-1]
            vmax = 0.0
    if kL is not None:
        Lt = 8/math.sqrt(kL); h = (x[1] - x[0]); ext = []
        while sum(ext) < Lt: ext.append(min(h, Lt/20)); h *= 1.3
        xs = [x[0] - c for c in np.cumsum(ext)][::-1]; x = xs + x; V = [kL]*len(ext) + V
    if kR is not None:
        Lt = 8/math.sqrt(kR); h = (x[-1] - x[-2]); ext = []
        while sum(ext) < Lt: ext.append(min(h, Lt/20)); h *= 1.3
        x = x + [x[-1] + c for c in np.cumsum(ext)]; V = V + [kR]*len(ext)
    x = np.asarray(x); V = np.asarray(V); h = np.diff(x)
    hl, hr = h[:-1], h[1:]; Vl, Vr = V[:-1], V[1:]
    m = (hl + hr)/2
    d = (1/hl + 1/hr + (Vl*hl + Vr*hr)/2)/m
    off = -(1/hr[:-1])/np.sqrt(m[:-1]*m[1:])
    w = eigh_tridiagonal(d, off, select='i', select_range=(0, 0), eigvals_only=True)
    return float(max(w[0], 0.0))

def bo_family_f(armA, armC, theta, K=160):
    """float approximation of bo2.bo_family(armA, armC, theta, K)[0] (0.0 when bo_family would give up)"""
    u = np.array([math.cos(theta), math.sin(theta)])
    arms = [a for a in (armA, armC) if a is not None]
    bd_lo, bd_hi = min(0.0, u[0]), max(0.0, u[0])
    cpts = [np.zeros(2), np.array([1.0, 0.0])]; signs = []; ppts = {}
    for a in arms:
        if a.kind == 'polar':
            P = _polar_pts(a)
            if P is None: return 0.0
            ppts[id(a)] = P; cpts += list(P)
        else:
            s = _dir_sign(a, u)
            if s == 0: return 0.0
            signs.append((a, s))
    prc = np.array([p @ u for p in cpts])
    xc_lo, xc_hi = prc.min(), prc.max(); span = max(xc_hi - xc_lo, 1e-6)
    x_cut_R = xc_hi + 0.05*span; x_cut_L = xc_lo - 0.05*span
    tails = {+1: [], -1: []}
    for a, s in signs:
        P = _far_pts(a, u, x_cut_R if s > 0 else x_cut_L)
        cb = _chord(a, u)
        if P is None or cb is None: return 0.0
        ppts[id(a)] = P; tails[s].append(cb)
    x_lo = x_cut_L if tails[-1] else xc_lo; x_hi = x_cut_R if tails[+1] else xc_hi
    nodes = graded_nodes(x_lo, x_hi, list(prc) + [x_cut_L, x_cut_R], K)
    mids_lo = np.asarray(nodes[:-1]); mids_hi = np.asarray(nodes[1:])
    meets = (mids_hi >= bd_lo) & (mids_lo <= bd_hi)
    exts = []
    for a in (armA, armC):
        if a is None: exts.append(np.zeros(len(nodes) - 1)); continue
        e = _ext(ppts[id(a)], u, nodes)
        cb = _chord(a, u)
        if cb is not None: e = np.where(meets, e, np.minimum(e, cb))
        exts.append(e)
    H = np.where(meets, exts[0] + exts[1], np.maximum(exts[0], exts[1]))
    kL = (math.pi/max(tails[-1]))**2 if tails[-1] else None
    kR = (math.pi/max(tails[+1]))**2 if tails[+1] else None
    return _mu(nodes, H, kL, kR)
