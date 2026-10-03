"""Float assembly of the BD-chart matrices from the exact reference data (refexact.build)."""
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh

def Gfun(X):
    x, y = float(X[0]), abs(float(X[1]))
    return np.array([y + x*x/y, -x/y, 1.0/y]), y          # (G11, G12, G22), J

class Chart:
    """Reference matrices as scipy CSR; 'free' dofs; dof partitions."""
    def __init__(self, R):
        self.R = R
        n = R['nedge']; self.ncr = n
        self.cr = {}
        for side in 'AC':
            for nm in ('11', '12', '22'):
                r, c, v = R['cr'][(side, nm)]; self.cr[(side, nm)] = sp.csr_matrix((v, (r, c)), shape=(n, n))
            self.cr[(side, 'M')] = R['cr'][(side, 'M')]
        self.crfree = np.where(~R['bnd_edge'])[0]
        n2 = R['nnode'] + R['nedge']; self.np2 = n2
        self.p2 = {}
        for side in 'AC':
            for nm in ('11', '12', '22', 'M'):
                r, c, v = R['p2'][(side, nm)]; self.p2[(side, nm)] = sp.csr_matrix((v, (r, c)), shape=(n2, n2))
        self.p2free = np.where(~R['p2_bnd'])[0]
        es = R['edge_side'][self.crfree]
        self.crpart = {'BD': np.where(es == 0)[0], 'A': np.where(es == 1)[0], 'C': np.where(es == 2)[0]}
    def cr_mats(self, A, C, free=True):
        gA, JA = Gfun(A); gC, JC = Gfun(C)
        K = sum(g*self.cr[(s, nm)] for s, gg in (('A', gA), ('C', gC)) for g, nm in zip(gg, ('11', '12', '22')))
        M = JA*self.cr[('A', 'M')] + JC*self.cr[('C', 'M')]
        if free:
            f = self.crfree; return K[f][:, f].tocsr(), M[f]
        return K, M
    def p2_mats(self, A, C):
        gA, JA = Gfun(A); gC, JC = Gfun(C)
        K = sum(g*self.p2[(s, nm)] for s, gg in (('A', gA), ('C', gC)) for g, nm in zip(gg, ('11', '12', '22')))
        M = JA*self.p2[('A', 'M')] + JC*self.p2[('C', 'M')]
        f = self.p2free; return K[f][:, f].tocsr(), M[f][:, f].tocsr()
    def hmax(self, A, C):
        ev = self.R['evecs']; h = 0.0
        for X in (A, C):
            F = np.array([[1.0, X[0]], [0.0, X[1]]]); h = max(h, np.max(np.linalg.norm(ev @ F.T, axis=1)))
        return h

def solve_cr(ch, A, C, nev=3):
    K, md = ch.cr_mats(A, C)
    w, V = eigsh(K, k=nev, M=sp.diags(md), sigma=0.0, which='LM'); o = np.argsort(w); return w[o], V[:, o], K, md
def solve_p2(ch, A, C, nev=3):
    K, M = ch.p2_mats(A, C)
    w, V = eigsh(K, k=nev, M=M, sigma=0.0, which='LM'); o = np.argsort(w); return w[o], V[:, o], K, M
