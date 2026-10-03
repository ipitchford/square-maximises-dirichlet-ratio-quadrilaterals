"""Verified inertia test for a symmetric sparse matrix via banded LDL^T without pivoting (float64) + componentwise
backward error bound (Higham, ASNA 2nd ed., Thm 9.3 and the remark on banded matrices; for symmetric elimination
U = D L^T):  L D L^T = A_fl + dA,  |dA| <= gamma_{b+1} |L||D||L^T|   (b = bandwidth, gamma_m = m u/(1 - m u)).
Returns the number of negative pivots and a rigorous upper bound for ||dA||_2 (<= max row sum of the bound matrix).
If no pivot is zero, Sylvester's law of inertia gives In(A_fl + dA) = In(D)."""
import numpy as np, scipy.sparse as sp
from scipy.sparse.csgraph import reverse_cuthill_mckee

U = 2.0**-53
def gamma(m): return m*U/(1 - m*U)

def ldl_band(A):
    """A: scipy sparse symmetric (float).  Returns dict(neg, zero, dA_bound, b, n, minabsd)."""
    A = sp.csr_matrix(A)
    perm = reverse_cuthill_mckee(A, symmetric_mode=True)
    Ap = A[perm][:, perm].tocoo()
    n = Ap.shape[0]
    b = int(np.max(np.abs(Ap.row - Ap.col)))   # RIGOR: exact: integer bandwidth
    # lower band storage: AB[d, j] = Ap[j+d, j]
    AB = np.zeros((b+1, n))
    m = Ap.row >= Ap.col
    AB[(Ap.row - Ap.col)[m], Ap.col[m]] = Ap.data[m]
    Lb = np.zeros((b+1, n)); dvec = np.zeros(n)
    W = np.zeros((b+1, b+1))
    # initial window: rows/cols 0..b
    def col_full(j):            # column j of the (updated later) matrix, entries j..j+b, from AB
        return AB[:, j]
    # fill W with the original entries of block [0..b]x[0..b]
    for jj in range(min(b+1, n)):   # RIGOR: exact: integers
        v = AB[:b+1-jj, jj]        # entries (jj+d, jj), d=0..b-jj
        W[jj:jj+len(v), jj] = v; W[jj, jj:jj+len(v)] = v
    for j in range(n):
        dj = W[0, 0]; dvec[j] = dj
        if dj == 0.0: return dict(neg=-1, zero=True)
        l = W[1:, 0]/dj
        Lb[1:, j] = l
        W[1:, 1:] -= dj*np.outer(l, l)   # RIGOR: accounted: every float operation of the factorisation is covered by the backward error bound |dA| <= gamma_{b+1}|L||D||L^T| (Higham, Thm 9.3)
        # shift window by one: new window covers j+1 .. j+1+b
        W[:-1, :-1] = W[1:, 1:]
        newi = j + 1 + b
        W[-1, :] = 0.0; W[:, -1] = 0.0
        if newi < n:
            # new row/col: original entries A[newi, newi-b .. newi] -> positions 0..b of window
            # A[newi, c] for c in [j+1, newi] = AB[newi - c, c]
            cs = np.arange(j+1, newi+1)
            vals = AB[newi - cs, cs]
            W[-1, :] = vals; W[:, -1] = vals
    neg = int(np.sum(dvec < 0))   # RIGOR: exact: count of negative pivots (the inertia)
    # bound matrix |L||D||L^T| row sums: row i: sum_k |L_ik| |d_k| colsum_k, colsum_k = sum_j |L_jk| (incl. diag 1)
    absL = np.abs(Lb); absL[0, :] = 1.0
    colsum = absL.sum(axis=0)                      # column k: sum over d of |L_{k+d,k}|   # RIGOR: accounted: float sums of nonnegative terms, covered by the factor 1 + 1e-6 of dA
    w = np.abs(dvec)*colsum                        # per column k
    rows = np.zeros(n)
    for d in range(b+1):
        # contributes |L_{k+d,k}| * w_k to row k+d
        rows[d:] += absL[d, :n-d]*w[:n-d]
    # each entry of the bound involves at most b+1 products; safety factor for the float evaluation of the bound
    dA = gamma(b+2)*np.max(rows)*(1 + 1e-6)   # RIGOR: accounted: max of computed row sums; the float error of the bound is covered by 1 + 1e-6
    return dict(neg=neg, zero=False, dA_bound=dA, b=b, n=n, minabsd=float(np.min(np.abs(dvec))), dvec=dvec)   # RIGOR: info: diagnostics

if __name__ == '__main__':
    # self-test on a random banded SPD-minus-shift matrix with a known inertia (dense check)
    rng = np.random.default_rng(0)
    n = 400
    B = sp.random(n, n, density=0.02, random_state=1); A = (B + B.T)/2 + sp.eye(n)*3
    A = A.tocsr(); w = np.linalg.eigvalsh(A.toarray())
    for s in (w[0] - 0.1, (w[0] + w[1])/2, (w[2] + w[3])/2):
        r = ldl_band(A - s*sp.eye(n)); print('shift between eigs: neg pivots', r['neg'], 'true', int(np.sum(w < s)), 'b', r['b'], 'dA', r['dA_bound'])   # RIGOR: info: self-test output
