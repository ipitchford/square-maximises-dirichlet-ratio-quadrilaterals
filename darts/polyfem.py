"""Non-rigorous P2/P3 FEM on unions of convex polygons (fan-triangulated from a 'centre' vertex, graded toward it).
Used for float pilots of inclusion-box bounds.  pieces: list of (centre, [v1, v2, ..., vk]) where the convex polygon
is centre, v1, ..., vk (ccw or cw) -- fan triangles (centre, v_i, v_{i+1})."""
import numpy as np
from scipy.sparse.linalg import eigsh
from skfem import MeshTri, ElementTriP2, ElementTriP3, ElementTriP1, Basis, asm, condense
from skfem.models.poisson import laplace, mass
from fem_ref import graded_tri

def poly_mesh(pieces, k=4, kappa=2.5):
    Xs = []; Ts = []; off = 0
    for (c, vs) in pieces:
        c = np.asarray(c, float)
        for i in range(len(vs) - 1):
            P1 = np.asarray(vs[i], float); P2 = np.asarray(vs[i+1], float)
            ar = (P1[0]-c[0])*(P2[1]-c[1]) - (P1[1]-c[1])*(P2[0]-c[0])
            if abs(ar) < 1e-14: continue
            X, T = graded_tri(c, P1, P2, k, kappa)
            Xs.append(X); Ts.append(T + off); off += X.shape[1]
    X = np.hstack(Xs); T = np.hstack(Ts)
    key = np.round(X.T*1e9).astype(np.int64)
    _, idx, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    P = X[:, idx]; T = inv.ravel()[T]
    a = P[:, T[0]]; b = P[:, T[1]]; cc = P[:, T[2]]
    sgn = (b[0]-a[0])*(cc[1]-a[1]) - (b[1]-a[1])*(cc[0]-a[0])
    T = T.copy(); sw = sgn < 0; T[1, sw], T[2, sw] = T[2, sw], T[1, sw].copy()
    return MeshTri(P, T)

def poly_eigs(pieces, k=4, kappa=2.5, deg=2, nev=3):
    m = poly_mesh(pieces, k, kappa)
    E = {1: ElementTriP1, 2: ElementTriP2, 3: ElementTriP3}[deg]()
    basis = Basis(m, E)
    K = asm(laplace, basis); M = asm(mass, basis)
    Kc, Mc = condense(K, M, D=basis.get_dofs(), expand=False)
    w = eigsh(Kc, k=nev, M=Mc, sigma=0.0, which='LM', return_eigenvectors=False)
    return np.sort(w)

def hull_ccw(pts):
    pts = sorted(set(map(tuple, np.round(np.asarray(pts, float), 15))))
    def cr(o, a, b): return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])
    lo = []; up = []
    for p in pts:
        while len(lo) >= 2 and cr(lo[-2], lo[-1], p) <= 0: lo.pop()
        lo.append(p)
    for p in reversed(pts):
        while len(up) >= 2 and cr(up[-2], up[-1], p) <= 0: up.pop()
        up.append(p)
    return [np.array(p) for p in lo[:-1] + up[:-1]]

def fan_from(c, hull):
    """hull: ccw convex polygon containing vertex c; return vertex list starting after c going around back to before c"""
    hull = [np.asarray(h, float) for h in hull]
    i = int(np.argmin([np.linalg.norm(h - c) for h in hull]))
    assert np.linalg.norm(hull[i] - c) < 1e-12
    n = len(hull); return [hull[(i + j) % n] for j in range(1, n)]
