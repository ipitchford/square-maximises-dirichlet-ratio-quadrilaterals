"""Reference meshes and reference matrices for darts in the BD-normalised frame B=(0,0), D=(1,0).
Triangle X in {A, C}: affine map F_X(xhat) = xhat1*(D-B) + xhat2*(X-B) from the reference triangle
(0,0)=B, (1,0)=D, (0,1)=X.  Mesh: k uniform refinements, then radial grading toward D=(1,0):
 s = 1 - lambda_D, s' = s^kappa (integer kappa), ratios lambda_B:lambda_X kept.  All reference coordinates are dyadic
rationals, so every reference matrix entry is an exact rational; we store floats plus a rigorous bound on the rounding.
Physical stiffness on triangle X:  K = sum_{kl} G_X[k,l] Khat_X^{(kl)},  G_X = |det F_X| F_X^{-1} F_X^{-T};
mass M = |det F_X| Mhat_X."""
import numpy as np
from fractions import Fraction as Fr
from skfem import MeshTri

def ref_points(k, kappa):
    m = MeshTri.init_refdom().refined(k)
    n = 2**k
    # exact dyadic coordinates: round to the grid 1/n
    P = np.rint(m.p * n).astype(np.int64)          # integer coordinates on grid
    T = m.t.copy()
    pts = []
    for (i, j) in P.T:
        lD = Fr(int(i), n); lX = Fr(int(j), n); lB = 1 - lD - lX
        s = lB + lX
        if s == 0:
            pts.append((Fr(1), Fr(0)))
        else:
            s2 = s**kappa
            pts.append((1 - s2, lX*s2/s))
    return pts, T

def mesh_both(k, kappa):
    """Return exact reference points for A-side and C-side (identical reference meshes), triangles, and the shared
    BD node identification.  Global node set: union, nodes with xhat2 == 0 (on BD) shared."""
    pts, T = ref_points(k, kappa)
    npt = len(pts)
    onBD = [i for i, p in enumerate(pts) if p[1] == 0]
    # global numbering: A-side nodes 0..npt-1; C-side nodes: BD nodes map to A-side ids, others new ids
    cmap = np.empty(npt, dtype=np.int64); nxt = npt
    bdset = set(onBD)
    for i in range(npt):
        if i in bdset: cmap[i] = i
        else: cmap[i] = nxt; nxt += 1
    return pts, T, cmap, nxt, onBD

def edges_of(T):
    """unique edges of triangle list T (3 x ne); returns edge array (2 x nE) and element->edge (3 x ne),
    local edge i is opposite local vertex i."""
    ne = T.shape[1]
    loc = [(1, 2), (2, 0), (0, 1)]
    E = {}
    t2e = np.empty((3, ne), dtype=np.int64)
    elist = []
    for e in range(ne):
        for i, (a, b) in enumerate(loc):
            u, v = int(T[a, e]), int(T[b, e])
            key = (u, v) if u < v else (v, u)
            if key not in E:
                E[key] = len(elist); elist.append(key)
            t2e[i, e] = E[key]
    return np.array(elist, dtype=np.int64).T, t2e
