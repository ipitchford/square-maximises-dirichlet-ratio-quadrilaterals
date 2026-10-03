"""Float prototype of the Lemma B' chart cover (thin band 1/10 < w/D <= 2/5, convex quadrilaterals).
Chart I: diameter = diagonal P0P2 = (0,0)-(1,0); P1=(x1,-y1) below, P3=(x3,y3) above.
Box bound: lambda2(Q) <= R2(Q_in) (Ritz, weighted 1-D form on inner profiles a_in<y<b_in),
           lambda1(Q) >= mu1(h_out)  (1-D Born-Oppenheimer lower bound on the outer slice length).
Here both are float P1 FEM approximations (lumped mass); this only estimates box counts."""
import numpy as np, sys, time
from scipy.linalg import eigh_tridiagonal, eigh
PI2 = np.pi**2
C6 = 1/6 + 1/(4*PI2)

def tent_max(x, l, u):  # max over apex abscissa in [l,u] of tent(x; apex), height 1
    return np.minimum(1.0, np.minimum(x/l, (1-x)/(1-u)))
def tent_min(x, l, u):  # min over apex abscissa in [l,u]
    return np.minimum(x/u, (1-x)/(1-l))

def mesh(kinks, N):
    xs = np.sort(np.concatenate([np.linspace(0, 1, N+1), np.clip(kinks, 0, 1)]))
    keep = np.concatenate([[True], np.diff(xs) > 1e-6]); xs = xs[keep]
    # re-insert exact kinks so profiles stay linear per element (replace the nearest node)
    for k in kinks:
        if 0 < k < 1: xs[np.argmin(np.abs(xs-k))] = k
    return np.sort(xs)

def mu1_lower_proxy(hfun, kinks, N=400):
    xs = mesh(kinks, N); xm = 0.5*(xs[1:]+xs[:-1]); dx = np.diff(xs)
    # 2-point Gauss on each element for potential with P1 hats; lumped mass
    g = np.array([-1, 1])/np.sqrt(3)
    Kd = np.zeros(len(xs)); Ko = np.zeros(len(xs)-1); Md = np.zeros(len(xs))
    for gp in g:
        xq = xm + 0.5*dx*gp; V = PI2/np.maximum(hfun(xq), 1e-300)**2; w = 0.5*dx
        phiL = 0.5*(1-gp); phiR = 0.5*(1+gp)
        Kd[:-1] += w*V*phiL**2; Kd[1:] += w*V*phiR**2; Ko += w*V*phiL*phiR
    Kd[:-1] += 1/dx; Kd[1:] += 1/dx; Ko += -1/dx
    Md[:-1] += dx/2; Md[1:] += dx/2
    d = Kd[1:-1]/Md[1:-1]; e = Ko[1:-1]/np.sqrt(Md[1:-2]*Md[2:-1])
    return eigh_tridiagonal(d, e, select='i', select_range=(0, 0), eigvals_only=True)[0]

def mu2_weighted(afun, bfun, kinks, N=400):
    xs = mesh(kinks, N); xm = 0.5*(xs[1:]+xs[:-1]); dx = np.diff(xs)
    a0, a1 = afun(xs[:-1]), afun(xs[1:]); b0, b1 = bfun(xs[:-1]), bfun(xs[1:])
    ap = (a1-a0)/dx; bp = (b1-b0)/dx; hp = bp-ap
    g = np.array([-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053])
    wg = np.array([0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454])
    Kd = np.zeros(len(xs)); Ko = np.zeros(len(xs)-1); Md = np.zeros(len(xs)); Mo = np.zeros(len(xs)-1)
    coef = 0.5 + 0.5*ap*bp + hp**2*C6   # a'^2/2 + a'h'/2 = a'b'/2
    for gp, wq in zip(g, wg):
        xq = xm + 0.5*dx*gp; h = bfun(xq) - afun(xq); w = 0.5*dx*wq
        pL = 0.5*(1-gp); pR = 0.5*(1+gp); dL = -1/dx; dR = 1/dx
        # f'^2 h/2 + f f' h'/2 + f^2 (pi^2/h) coef
        def form(p, dp, q, dq): return dp*dq*h/2 + 0.5*(p*dq + q*dp)*hp/2 + p*q*PI2/h*coef
        Kd[:-1] += w*form(pL, dL, pL, dL); Kd[1:] += w*form(pR, dR, pR, dR); Ko += w*form(pL, dL, pR, dR)
        Md[:-1] += w*pL**2*h/2; Md[1:] += w*pR**2*h/2; Mo += w*pL*pR*h/2
    n = len(xs)-2
    K = np.diag(Kd[1:-1]) + np.diag(Ko[1:-1], 1) + np.diag(Ko[1:-1], -1)
    M = np.diag(Md[1:-1]) + np.diag(Mo[1:-1], 1) + np.diag(Mo[1:-1], -1)
    return eigh(K, M, subset_by_index=[0, 1], eigvals_only=True)[1]

def chartI_box_ratio(x1l, x1u, y1l, y1u, x3l, x3u, y3l, y3u, N=400):
    hout = lambda x: y1u*tent_max(x, x1l, x1u) + y3u*tent_max(x, x3l, x3u)
    ain = lambda x: -y1l*tent_min(x, x1l, x1u)
    bin_ = lambda x: y3l*tent_min(x, x3l, x3u)
    kinks = np.array([x1l, x1u, x3l, x3u, x1u/(1+x1u-x1l), x3u/(1+x3u-x3l)])
    return mu2_weighted(ain, bin_, kinks, N) / mu1_lower_proxy(hout, kinks, N)

def width_chartI(x1, y1, x3, y3):
    V = np.array([[0, 0], [x1, -y1], [1, 0], [x3, y3]])
    w = np.inf
    for i in range(4):
        e = V[(i+1) % 4]-V[i]; n = np.array([-e[1], e[0]])/np.hypot(*e)
        w = min(w, np.max(np.abs((V-V[i]) @ n)))
    return w

if __name__ == "__main__":
    # point checks against the thin study (vertices given there are diameter-normalised)
    for p in [(0.4, 0.1, 0.5, 0.3), (0.3, 0.2, 0.45, 0.25), (0.6, 0.15, 0.3, 0.35)]:
        x1, y1, x3, y3 = p; t = time.time()
        r0 = chartI_box_ratio(x1, x1, y1, y1, x3, x3, y3, y3)
        print("point", p, "w", round(width_chartI(*p), 4), "1-D ratio", round(r0, 4), "%.1f ms" % (1e3*(time.time()-t)))
        for rel in [0.005, 0.01, 0.02, 0.04]:
            dy1, dy3 = rel*y1, rel*y3; dx = rel
            r = chartI_box_ratio(x1-dx, x1+dx, y1-dy1, y1+dy1, x3-dx, x3+dx, y3-dy3, y3+dy3)
            print("   rel half-width %.3f: box ratio %.4f (loss %.3f)" % (rel, r, r/r0-1))

# ---------------- chart II: diameter = side P0P1 = (0,0)-(1,0); P2=(x2,y2), P3=(x3,y3) above, 0<x3<x2<1
def upper_hull_profile(points):
    """upper concave envelope over [0,1] of a point set containing (0,0),(1,0); returns callable + kinks"""
    P = sorted(set(map(tuple, points)))
    H = []
    for p in P:
        while len(H) >= 2 and (H[-1][0]-H[-2][0])*(p[1]-H[-2][1]) - (H[-1][1]-H[-2][1])*(p[0]-H[-2][0]) >= 0: H.pop()
        H.append(p)
    hx = np.array([q[0] for q in H]); hy = np.array([q[1] for q in H])
    return (lambda x: np.interp(x, hx, hy)), hx[1:-1]

def chartII_box_ratio(x2l, x2u, y2l, y2u, x3l, x3u, y3l, y3u, N=400):
    corners = [(0.0, 0.0), (1.0, 0.0)] + [(x, y) for x in (x2l, x2u) for y in (y2l, y2u)] + [(x, y) for x in (x3l, x3u) for y in (y3l, y3u)]
    bout, k1 = upper_hull_profile(np.array(corners))
    def apex(xl, xu, yl):
        xs = xu/(1+xu-xl); return (xs, yl*xs/xu)
    bin_, k2 = upper_hull_profile(np.array([(0.0, 0.0), (1.0, 0.0), apex(x2l, x2u, y2l), apex(x3l, x3u, y3l)]))
    zero = lambda x: 0.0*x
    kinks = np.concatenate([k1, k2])
    return mu2_weighted(zero, bin_, kinks, N) / mu1_lower_proxy(bout, kinks, N)
