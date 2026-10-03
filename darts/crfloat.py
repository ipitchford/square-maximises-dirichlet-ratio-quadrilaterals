"""FLOAT (non-rigorous) CR / P2 assembly on the BD-normalised dart via reference matrices (for pilots and as the
numerical core of the rigorous certificate; rigour is added in cert.py)."""
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from refmesh import mesh_both, edges_of
from skfem import MeshTri, ElementTriP2, Basis, BilinearForm, asm

class RefDart:
    def __init__(self, k=5, kappa=2):
        pts, T, cmap, nglob, onBD = mesh_both(k, kappa)
        self.k, self.kappa = k, kappa
        self.Pref = np.array([[float(x), float(y)] for x, y in pts]).T     # 2 x npt (reference coords)
        self.T = T; self.cmap = cmap; self.nglob = nglob; self.npt = len(pts)
        # global triangles: A-side uses ids as is, C-side via cmap
        self.TA = T; self.TC = cmap[T]
        TT = np.hstack([self.TA, self.TC])
        E, t2e = edges_of(TT)
        self.E = E; self.t2e = t2e; self.ne = T.shape[1]
        nE = E.shape[1]
        # boundary edges: edges belonging to exactly one element
        cnt = np.bincount(t2e.ravel(), minlength=nE); self.bnd = np.where(cnt == 1)[0]
        self.free = np.setdiff1d(np.arange(nE), self.bnd)
        # reference element geometry (same reference mesh for both sides)
        P = self.Pref
        x = P[:, T]   # 2 x 3 x ne
        # edge vectors opposite vertex i: e_i = x_{i+2} - x_{i+1}
        ev = np.stack([x[:, (i+2) % 3, :] - x[:, (i+1) % 3, :] for i in range(3)], axis=1)   # 2 x 3 x ne
        det = (x[0, 1] - x[0, 0])*(x[1, 2] - x[1, 0]) - (x[1, 1] - x[1, 0])*(x[0, 2] - x[0, 0])
        self.aref = np.abs(det)/2           # reference element areas
        # grad lambda_i = R(e_i)/(det) with R(x,y) = (-y, x)?  For ccw: grad lambda_i = (-e_i_y, e_i_x)/det  (sign det)
        g = np.stack([-ev[1], ev[0]], axis=0)/det        # 2 x 3 x ne  (works for either orientation)
        self.gl = g                                      # reference grad of barycentrics
        self.diamref = np.max(np.linalg.norm(ev, axis=0), axis=0)
        self._build_cr(); self._build_p2()
    def _build_cr(self):
        ne = self.ne; g = self.gl; a = self.aref; nE = self.E.shape[1]
        # CR basis phi_i = 1 - 2 lambda_i on edge opposite vertex i: grad = -2 grad lambda_i
        mats = {}
        for side, off in (('A', 0), ('C', ne)):
            rows = self.t2e[:, off:off+ne]
            for (k, l) in ((0, 0), (0, 1), (1, 1)):
                vals = np.empty((3, 3, ne))
                for i in range(3):
                    for j in range(3):
                        if (k, l) == (0, 1):
                            vals[i, j] = 4*a*(g[0, i]*g[1, j] + g[1, i]*g[0, j])   # symmetric part (G12 multiplies both)
                        else:
                            vals[i, j] = 4*a*g[k, i]*g[l, j]
                R = np.repeat(rows[:, None, :], 3, axis=1); C = np.repeat(rows[None, :, :], 3, axis=0)
                mats[(side, k, l)] = sp.csr_matrix((vals.ravel(), (R.ravel(), C.ravel())), shape=(nE, nE))
            md = np.zeros(nE); np.add.at(md, rows.ravel(), np.tile(a/3, 3))
            mats[(side, 'M')] = md
        self.cr = mats
    def _build_p2(self):
        # P2 on the global node mesh: use skfem on the reference mesh twice (A and C side) with custom forms
        ne = self.ne
        m = MeshTri(self.Pref, self.T)
        basis = Basis(m, ElementTriP2())
        forms = {(0, 0): BilinearForm(lambda u, v, w: u.grad[0]*v.grad[0]),
                 (0, 1): BilinearForm(lambda u, v, w: u.grad[0]*v.grad[1] + u.grad[1]*v.grad[0]),
                 (1, 1): BilinearForm(lambda u, v, w: u.grad[1]*v.grad[1])}
        Kh = {kl: asm(f, basis) for kl, f in forms.items()}
        Mh = asm(BilinearForm(lambda u, v, w: u*v), basis)
        # dof coordinates to merge A/C sides on BD: P2 dofs = nodes + edge midpoints
        X = basis.doflocs                     # 2 x ndof (reference)
        nd = X.shape[1]
        onBD = np.abs(X[1]) < 1e-14
        cmap = np.empty(nd, dtype=np.int64); nxt = nd
        for i in range(nd):
            if onBD[i]: cmap[i] = i
            else: cmap[i] = nxt; nxt += 1
        Pm = sp.csr_matrix((np.ones(nd), (cmap, np.arange(nd))), shape=(nxt, nd))
        PA = sp.csr_matrix((np.ones(nd), (np.arange(nd), np.arange(nd))), shape=(nxt, nd))
        self.p2 = {}
        for kl in Kh:
            self.p2[('A',) + kl] = (PA @ Kh[kl] @ PA.T).tocsr(); self.p2[('C',) + kl] = (Pm @ Kh[kl] @ Pm.T).tocsr()
        self.p2[('A', 'M')] = (PA @ Mh @ PA.T).tocsr(); self.p2[('C', 'M')] = (Pm @ Mh @ Pm.T).tocsr()
        # boundary dofs: reference boundary of the triangle except the BD edge (x2 = 0) for each side
        bref = basis.get_dofs().all()
        bd = np.zeros(nd, bool); bd[bref] = True
        bA = bd & ~onBD; bC = bd & ~onBD
        bglob = np.zeros(nxt, bool); bglob[np.arange(nd)[bA]] = True; bglob[cmap[bC]] = True
        # B and D endpoints are on the boundary of the dart too (vertices of BD): x2==0 and x1 in {0,1}
        endpts = onBD & ((np.abs(X[0]) < 1e-14) | (np.abs(X[0]-1) < 1e-14)); bglob[np.arange(nd)[endpts]] = True
        self.p2free = np.where(~bglob)[0]; self.np2 = nxt

def Gmat(X):
    """X = apex (x, y); F = [[1, x],[0, y]]; G = |det F| F^{-1} F^{-T}, J = |det F|"""
    x, y = X; J = abs(y)
    Fi = np.array([[1.0, -x/y], [0.0, 1.0/y]])
    return J*Fi @ Fi.T, J

def assemble(ref, A, C, kind='cr'):
    GA, JA = Gmat(A); GC, JC = Gmat(C)
    if kind == 'cr':
        d = ref.cr
        K = (GA[0, 0]*d[('A', 0, 0)] + GA[0, 1]*d[('A', 0, 1)] + GA[1, 1]*d[('A', 1, 1)]
             + GC[0, 0]*d[('C', 0, 0)] + GC[0, 1]*d[('C', 0, 1)] + GC[1, 1]*d[('C', 1, 1)])
        M = JA*d[('A', 'M')] + JC*d[('C', 'M')]
        f = ref.free
        return K[f][:, f], sp.diags(M[f])
    else:
        d = ref.p2
        K = (GA[0, 0]*d[('A', 0, 0)] + GA[0, 1]*d[('A', 0, 1)] + GA[1, 1]*d[('A', 1, 1)]
             + GC[0, 0]*d[('C', 0, 0)] + GC[0, 1]*d[('C', 0, 1)] + GC[1, 1]*d[('C', 1, 1)])
        M = JA*d[('A', 'M')] + JC*d[('C', 'M')]
        f = ref.p2free
        return K[f][:, f], M[f][:, f]

def hmax(ref, A, C):
    out = 0.0
    for X in (A, C):
        F = np.array([[1.0, X[0]], [0.0, X[1]]])
        P = F @ ref.Pref; x = P[:, ref.T]
        ev = np.stack([x[:, (i+2) % 3, :] - x[:, (i+1) % 3, :] for i in range(3)], axis=1)
        out = max(out, np.max(np.linalg.norm(ev, axis=0)))
    return out

def eig2(K, M, nev=3):
    w, V = eigsh(K, k=nev, M=M, sigma=0.0, which='LM')
    o = np.argsort(w); return w[o], V[:, o]
