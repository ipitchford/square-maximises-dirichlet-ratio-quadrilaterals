"""Float (NON-rigorous) Born-Oppenheimer / 1-D lower bound for lambda_1 of a simple polygon, used for planning only.

lambda_1(Omega) >= mu_1( -d^2/dx^2 + pi^2/l(x)^2 , Dirichlet on the projection interval ), l(x) = length of the LONGEST
component of the slice Omega cap {x = const} (x = coordinate along the unit vector u).  Here: l(x) on a grid by exact
line/polygon intersection, mu_1 by a piecewise-constant MINORANT potential on cells (max of l over each cell, from the cell ends
and the kinks inside) and Sturm shooting (bisection on E) -- the same scheme as the rigorous version, in floats."""
import numpy as np

def slice_components(P, x, u):
    """lengths of the components of P cap {<., u> = x} for a simple polygon P (array (n,2)), u unit"""
    v = np.array([-u[1], u[0]])
    X = P @ u; Y = P @ v
    ys = []
    n = len(P)
    for i in range(n):
        x0, x1 = X[i], X[(i + 1) % n]
        if (x0 - x)*(x1 - x) < 0:
            s = (x - x0)/(x1 - x0); ys.append(Y[i] + s*(Y[(i + 1) % n] - Y[i]))
        elif x0 == x and x1 != x:
            # vertex exactly on the line: count it once (by the half-open rule on the edge starting there)
            if (X[(i - 1) % n] - x)*(x1 - x) < 0: ys.append(Y[i])
    ys.sort()
    return [ys[2*j + 1] - ys[2*j] for j in range(len(ys)//2)]

def longest(P, x, u):
    c = slice_components(P, x, u)
    return max(c) if c else 0.0

def kinks(P, u):
    return np.sort(P @ u)

def cell_heights(P, u, K=400, grade=0.15):
    """cells over the projection interval, geometrically graded towards both ends (where l -> 0); H_k = max of l over the cell
    (l is piecewise linear between vertex projections: max at cell ends or kinks inside)"""
    X = P @ u; x0, x1 = X.min(), X.max(); L = x1 - x0
    # graded nodes: uniform in a stretched variable, plus all kinks
    s = np.linspace(0, 1, K + 1)
    g = 0.5*(1 - np.cos(np.pi*s))                      # cosine grading (dense at both ends)
    nodes = np.unique(np.concatenate([x0 + L*g, np.sort(X)]))
    H = np.zeros(len(nodes) - 1)
    ks = np.sort(X)
    for k in range(len(nodes) - 1):
        a, b = nodes[k], nodes[k + 1]
        pts = [a + 1e-12*L, b - 1e-12*L] + [z for z in ks if a < z < b]
        H[k] = max(longest(P, z, u) for z in pts)
    return nodes, H

def no_zero(nodes, H, E, tail=None):
    """Sturm: solution of -v'' + (V - E) v = 0, v(x0)=0, v'(x0)=1 stays positive on (x0, xK] (and, with tail = V_inf, on the
    constant-potential half line beyond xK: v'(xK) + sqrt(V_inf - E) v(xK) > 0).  Floats."""
    v, dv = 0.0, 1.0
    for k in range(len(H)):
        L = nodes[k + 1] - nodes[k]
        Vk = np.pi**2/H[k]**2 if H[k] > 0 else 1e30
        q = Vk - E
        if q > 0:
            kap = np.sqrt(q); th = np.tanh(kap*L)          # (v, dv) scaled by 1/cosh(kap L) > 0 (sign-preserving)
            vn = v + dv*th/kap; dvn = v*kap*th + dv
        else:
            om = np.sqrt(-q); c = np.cos(om*L); s = np.sin(om*L)
            # zero inside the cell?
            if v == 0:
                if om*L >= np.pi: return False
            else:
                ph = np.arctan2(dv/om, v)
                if np.pi/2 + ph - om*L <= 0: return False
            vn = v*c + dv*s/om; dvn = -v*om*s + dv*c
        if vn <= 0: return False
        v, dv = vn, dvn
        m = max(abs(v), abs(dv)); v /= m; dv /= m
    if tail is not None:
        if tail <= E: return False
        return dv + np.sqrt(tail - E)*v > 0
    return True

def mu1(P, u, K=400, tail=None):
    nodes, H = cell_heights(P, u, K)
    lo = 0.0
    if not no_zero(nodes, H, lo, tail): return 0.0
    hi = np.pi**2/max(H.max(), 1e-12)**2 + 1.0
    while no_zero(nodes, H, hi, tail) and hi < 1e12: hi *= 2
    for _ in range(50):
        mid = 0.5*(lo + hi)
        if no_zero(nodes, H, mid, tail): lo = mid
        else: hi = mid
    return lo

def dirs_for(P):
    """candidate slicing axes: all edge and diagonal directions of the polygon"""
    out = []
    n = len(P)
    for i in range(n):
        for j in range(i + 1, n):
            d = P[j] - P[i]; nd = np.hypot(*d)
            if nd > 1e-12: out.append(d/nd)
    return out

def bo_lower(P, K=300, extra_dirs=()):
    best = 0.0; bu = None
    for u in list(dirs_for(P)) + list(extra_dirs):
        m = mu1(P, u, K)
        if m > best: best, bu = m, u
    return best, bu
