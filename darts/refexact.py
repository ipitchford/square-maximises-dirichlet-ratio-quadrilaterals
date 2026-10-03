"""EXACT reference matrices for the BD-normalised dart chart (rigorous input data).

Reference triangle Rhat = {(0,0)=B, (1,0)=D, (0,1)=X}; mesh = k uniform refinements + radial grading toward D with
integer exponent kappa (s = 1 - lambda_D -> s^kappa); all node coordinates are exact dyadic rationals.
For both triangles X in {A, C} the SAME reference mesh is used; the global dart mesh glues the two copies along
BD (xhat2 = 0).  Physical matrices:  K(p) = sum_X sum_{kl} G_X[kl](p) Khat^{kl}_X,  M(p) = sum_X J_X(p) Mhat_X,
G_X = |det F_X| F_X^{-1}F_X^{-T}, J_X = |det F_X|, F_X = [[1, x],[0, y]] for apex X = (x, y).
Khat^{11}: int d1u d1v,  Khat^{12}: int (d1u d2v + d2u d1v),  Khat^{22}: int d2u d2v  (reference coordinates).
Every global entry is computed exactly (fractions.Fraction) and rounded once to the nearest double, so
|stored - exact| <= u |exact|, u = 2^-53 (we use the bound 1.2e-16*|stored|, which dominates).
Element families: 'cr' (Crouzeix-Raviart, dofs = edges, mass diagonal) and 'p2' (Lagrange P2, dofs = nodes + edges)."""
import numpy as np, pickle, os, time
from fractions import Fraction as Fr
from skfem import MeshTri

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'refcache')

def ref_mesh_exact(k, kappa):
    m = MeshTri.init_refdom().refined(k)
    n = 2**k
    P = np.rint(m.p*n).astype(np.int64)   # RIGOR: exact: the refined reference nodes are multiples of 1/n (dyadic), rint of exact integers
    pts = []
    for (i, j) in P.T:
        lD = Fr(int(i), n); lX = Fr(int(j), n); lB = 1 - lD - lX; s = lB + lX
        if s == 0: pts.append((Fr(1), Fr(0)))
        else:
            sk = s**kappa; pts.append((1 - sk, lX*sk/s))   # RIGOR: exact: Fraction ** int
    T = m.t.copy()
    # orient every reference triangle counter-clockwise (exact test)
    for e in range(T.shape[1]):
        a, b, c = (pts[T[i, e]] for i in range(3))
        d = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
        assert d != 0
        if d < 0: T[1, e], T[2, e] = T[2, e], T[1, e]
    return pts, T

def grad_bary(a, b, c):
    """exact gradients of barycentric coordinates for ccw triangle a,b,c; returns ([g0,g1,g2], area)"""
    D = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
    g0 = ((b[1]-c[1])/D, (c[0]-b[0])/D)
    g1 = ((c[1]-a[1])/D, (a[0]-c[0])/D)
    g2 = ((a[1]-b[1])/D, (b[0]-a[0])/D)
    return [g0, g1, g2], D/2

# ---- universal P2 integrals (divided by area) ----
def _mono_int(a, b, c):
    from math import factorial as f
    return Fr(2*f(a)*f(b)*f(c), f(a+b+c+2))
def _polymul(p, q):
    r = {}
    for ka, va in p.items():
        for kb, vb in q.items():
            k = (ka[0]+kb[0], ka[1]+kb[1], ka[2]+kb[2]); r[k] = r.get(k, 0) + va*vb
    return r
def _int(p): return sum(v*_mono_int(*k) for k, v in p.items())
def _e(i):
    t = [0, 0, 0]; t[i] = 1; return tuple(t)
def p2_universal():
    one = {(0, 0, 0): Fr(1)}
    lam = [{_e(i): Fr(1)} for i in range(3)]
    # basis: vertices 0,1,2 then edges (1,2),(2,0),(0,1)
    basis = []
    for i in range(3): basis.append(_polymul(lam[i], {_e(i): Fr(2), (0, 0, 0): Fr(-1)}))
    edges = [(1, 2), (2, 0), (0, 1)]
    for (i, j) in edges: basis.append({tuple(np.add(_e(i), _e(j))): Fr(4)})
    # gradient coefficients: grad phi_a = sum_m P[a][m] grad lambda_m, P linear polys
    P = [[{} for _ in range(3)] for _ in range(6)]
    for i in range(3): P[i][i] = {_e(i): Fr(4), (0, 0, 0): Fr(-1)}
    for t, (i, j) in enumerate(edges):
        P[3+t][i] = {_e(j): Fr(4)}; P[3+t][j] = {_e(i): Fr(4)}
    S = np.empty((6, 6, 3, 3), dtype=object)
    for a in range(6):
        for b in range(6):
            for m in range(3):
                for n in range(3):
                    S[a, b, m, n] = _int(_polymul(P[a][m], P[b][n])) if (P[a][m] and P[b][n]) else Fr(0)
    Mu = np.empty((6, 6), dtype=object)
    for a in range(6):
        for b in range(6): Mu[a, b] = _int(_polymul(basis[a], basis[b]))
    return S, Mu, edges

def build(k, kappa=2):   # RIGOR-FUNCTION: heuristic: legacy float assembly of the first session; not used by any replay (chart.py uses only grad_bary and p2_universal)
    os.makedirs(CACHE, exist_ok=True)
    fn = os.path.join(CACHE, f'ref_k{k}_kap{kappa}.pkl')
    if os.path.exists(fn):
        with open(fn, 'rb') as f: return pickle.load(f)
    t0 = time.time()
    pts, T = ref_mesh_exact(k, kappa)
    npt = len(pts); ne = T.shape[1]
    onBD = np.array([p[1] == 0 for p in pts])
    # global node ids: A copy = 0..npt-1; C copy: BD nodes shared, others new
    cmap = np.empty(npt, dtype=np.int64); nxt = npt
    for i in range(npt):
        if onBD[i]: cmap[i] = i
        else: cmap[i] = nxt; nxt += 1
    nnode = nxt
    # global edges
    loc = [(1, 2), (2, 0), (0, 1)]
    def edges_for(Tg):
        out = np.empty((3, ne), dtype=object)
        for e in range(ne):
            for i, (a, b) in enumerate(loc):
                u, v = int(Tg[a, e]), int(Tg[b, e]); out[i, e] = (u, v) if u < v else (v, u)
        return out
    TA = T; TC = cmap[T]
    eA = edges_for(TA); eC = edges_for(TC)
    eid = {}
    for arr in (eA, eC):
        for key in arr.ravel():
            if key not in eid: eid[key] = len(eid)
    nedge = len(eid)
    t2eA = np.vectorize(lambda kk: eid[kk])(eA).astype(np.int64)
    t2eC = np.vectorize(lambda kk: eid[kk])(eC).astype(np.int64)
    cnt = np.bincount(np.concatenate([t2eA.ravel(), t2eC.ravel()]), minlength=nedge)
    bnd_edge = cnt == 1
    edge_nodes = np.empty((2, nedge), dtype=np.int64)
    for key, i in eid.items(): edge_nodes[:, i] = key
    # node coordinates as floats (for diameters; exact values in pts / C copy mirrors)
    S2, Mu2, p2edges = p2_universal()
    geo = [grad_bary(*(pts[T[i, e]] for i in range(3))) for e in range(ne)]
    # ---------------- CR ----------------
    cr = {}
    for side, t2e in (('A', t2eA), ('C', t2eC)):
        acc = {('11'): {}, ('12'): {}, ('22'): {}}; mdiag = {}
        for e in range(ne):
            g, ar = geo[e]
            for i in range(3):
                mdiag[t2e[i, e]] = mdiag.get(t2e[i, e], 0) + ar/3
                for j in range(3):
                    key = (int(t2e[i, e]), int(t2e[j, e]))
                    # CR phi_i = 1 - 2 lambda_i  -> grad = -2 grad lambda_i
                    v11 = 4*ar*g[i][0]*g[j][0]; v22 = 4*ar*g[i][1]*g[j][1]
                    v12 = 4*ar*(g[i][0]*g[j][1] + g[i][1]*g[j][0])
                    for nm, v in (('11', v11), ('12', v12), ('22', v22)):
                        d = acc[nm]; d[key] = d.get(key, 0) + v
        for nm in acc:
            items = acc[nm]; r = np.array([kk[0] for kk in items]); c = np.array([kk[1] for kk in items])
            v = np.array([float(x) for x in items.values()])
            cr[(side, nm)] = (r, c, v)
        md = np.zeros(nedge); 
        for i, v in mdiag.items(): md[i] = float(v)
        cr[(side, 'M')] = md
    # ---------------- P2 ----------------
    # global P2 dofs: nodes 0..nnode-1, then edges nnode + edge id
    p2 = {}
    for side, Tg, t2e in (('A', TA, t2eA), ('C', TC, t2eC)):
        acc = {'11': {}, '12': {}, '22': {}, 'M': {}}
        for e in range(ne):
            g, ar = geo[e]
            # local dofs: vertices 0,1,2 then edges (1,2),(2,0),(0,1) = local edge index 0,1,2 (opposite vertex)
            dofs = [int(Tg[0, e]), int(Tg[1, e]), int(Tg[2, e])] + [nnode + int(t2e[i, e]) for i in range(3)]
            for a in range(6):
                for b in range(6):
                    v11 = 0; v12 = 0; v22 = 0
                    for m in range(3):
                        for n in range(3):
                            s = S2[a, b, m, n]
                            if s == 0: continue
                            v11 += s*g[m][0]*g[n][0]; v22 += s*g[m][1]*g[n][1]
                            v12 += s*(g[m][0]*g[n][1] + g[m][1]*g[n][0])
                    key = (dofs[a], dofs[b])
                    for nm, v in (('11', v11*ar), ('12', v12*ar), ('22', v22*ar), ('M', Mu2[a, b]*ar)):
                        d = acc[nm]; d[key] = d.get(key, 0) + v
        for nm in acc:
            items = acc[nm]; r = np.array([kk[0] for kk in items]); c = np.array([kk[1] for kk in items])
            v = np.array([float(x) for x in items.values()])
            p2[(side, nm)] = (r, c, v)
    # P2 boundary dofs: nodes on the dart boundary and boundary edges
    nodes_on_bnd = np.zeros(nnode, bool)
    for i in range(nedge):
        if bnd_edge[i]: nodes_on_bnd[edge_nodes[:, i]] = True
    p2_bnd = np.concatenate([nodes_on_bnd, bnd_edge])
    # partition of CR dofs: edges with both endpoints on BD are the shared ones
    xy = np.array([[float(p[0]), float(p[1])] for p in pts]).T
    xyC_ids = cmap  # C copy node i -> global id cmap[i]
    node_side = np.zeros(nnode, dtype=np.int8)   # 0 = BD, 1 = A-only, 2 = C-only
    node_side[:npt][~onBD] = 1
    node_side[cmap[~onBD]] = 2
    edge_side = np.empty(nedge, dtype=np.int8)
    for i in range(nedge):
        s0, s1 = node_side[edge_nodes[0, i]], node_side[edge_nodes[1, i]]
        edge_side[i] = 0 if (s0 == 0 and s1 == 0) else (1 if 1 in (s0, s1) else 2)
    # reference diameters per element (exact squared edge lengths -> float upper bound) and edge vectors for h(p)
    evecs = []
    for e in range(ne):
        a, b, c = (pts[T[i, e]] for i in range(3))
        for (u, v) in ((a, b), (b, c), (c, a)): evecs.append((float(v[0]-u[0]), float(v[1]-u[1])))
    evecs = np.array(evecs)
    out = dict(k=k, kappa=kappa, npt=npt, ne=ne, nnode=nnode, nedge=nedge, T=T, cmap=cmap, onBD=onBD,
               xy_ref=xy, edge_nodes=edge_nodes, bnd_edge=bnd_edge, edge_side=edge_side, node_side=node_side,
               t2eA=t2eA, t2eC=t2eC, cr=cr, p2=p2, p2_bnd=p2_bnd, evecs=evecs, build_time=time.time()-t0)
    with open(fn, 'wb') as f: pickle.dump(out, f)
    return out

if __name__ == '__main__':
    import sys
    for k in map(int, sys.argv[1:]):
        R = build(k); print(f"k={k}: nedge={R['nedge']} nnode={R['nnode']} build {R['build_time']:.1f}s", flush=True)
