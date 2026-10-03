"""Reference (NON-rigorous) FEM for Dirichlet eigenvalues of darts: P_k Lagrange elements on a mesh of
ABD u BCD graded radially toward the reflex vertex D (r -> r^kappa).  Used only for validation/pilots."""
import numpy as np
from scipy.sparse.linalg import eigsh
from skfem import MeshTri, ElementTriP2, ElementTriP3, ElementTriP1, Basis, asm, condense
from skfem.models.poisson import laplace, mass

_REF = {}
def ref_mesh(k):
    if k not in _REF:
        m = MeshTri.init_refdom().refined(k); _REF[k] = (m.p.copy(), m.t.copy())
    return _REF[k]

def graded_tri(D, P1, P2, k, kappa):
    p, t = ref_mesh(k); u, v = p[0], p[1]; s = u + v
    with np.errstate(invalid='ignore', divide='ignore'):
        f = np.where(s > 0, s**(kappa-1.0), 0.0)
    u2, v2 = u*f, v*f
    X = D[:, None] + np.outer(P1 - D, u2) + np.outer(P2 - D, v2)
    return X, t

def dart_mesh(A, B, C, D, k=5, kappa=2.0):
    X1, T1 = graded_tri(D, B, A, k, kappa)
    X2, T2 = graded_tri(D, B, C, k, kappa)
    X = np.hstack([X1, X2]); T = np.hstack([T1, T2 + X1.shape[1]])
    key = np.round(X.T*1e9).astype(np.int64)
    _, idx, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    P = X[:, idx]; T = inv.ravel()[T]
    # orient elements ccw
    a = P[:, T[0]]; b = P[:, T[1]]; c = P[:, T[2]]
    sgn = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
    T = T.copy(); sw = sgn < 0; T[1, sw], T[2, sw] = T[2, sw], T[1, sw].copy()
    return MeshTri(P, T)

def eigs(A, B, C, D, k=5, kappa=2.0, deg=2, nev=3):
    m = dart_mesh(A, B, C, D, k, kappa)
    E = {1: ElementTriP1, 2: ElementTriP2, 3: ElementTriP3}[deg]()
    basis = Basis(m, E)
    K = asm(laplace, basis); M = asm(mass, basis)
    Kc, Mc = condense(K, M, D=basis.get_dofs(), expand=False)
    w = eigsh(Kc, k=nev, M=Mc, sigma=0.0, which='LM', return_eigenvectors=False)
    return np.sort(w)

def xi(A, B, C, D, **kw):
    w = eigs(A, B, C, D, **kw); return w[1]/w[0], w
