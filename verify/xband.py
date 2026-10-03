"""Independent verifier for Lemma B' box certificates (own code; python-flint arb; exact rational geometry).
For a chart box B (BAND_SPEC.md sections 1-3):
  lambda1(Q) >= E   via Sturm shooting for -v'' + (V - E) v = 0 with the cellwise-constant minorant V = pi^2/max_cell(h_out)^2,
  lambda2(Q) <= R2  via the weighted 1-D Ritz form on the inner domain Q_in = {a_in < y < b_in}, exact element integrals,
and checks R2 < (5/2) E.  Floats are used only to choose E and the Ritz test vectors."""
import math, sys
from fractions import Fraction as Fr
import numpy as np
from scipy.linalg import eigh
from flint import arb, ctx
ctx.prec = 200
PI = arb.pi(); PI2 = PI*PI; C6 = arb(1)/6 + 1/(4*PI2)

def A(x):                      # exact rational -> arb
    x = Fr(x); return arb(x.numerator)/arb(x.denominator)

# ---------------- polylines with exact rational breakpoints -----------------
def collapse(pts, keep):
    """merge points with equal abscissa: keep='max' for OUTER (upper) profiles, 'min' for INNER profiles."""
    d = {}
    for x, y in pts:
        x, y = Fr(x), Fr(y)
        d[x] = y if x not in d else (max(d[x], y) if keep == "max" else min(d[x], y))
    return sorted(d.items())

class Poly:
    def __init__(self, pts):     # list of (x, y) rationals, STRICTLY increasing x, covering [0,1]
        self.p = [(Fr(x), Fr(y)) for x, y in pts]
        assert all(a[0] < b[0] for a, b in zip(self.p, self.p[1:])), "polyline abscissae must be strictly increasing"
        assert self.p[0][0] == 0 and self.p[-1][0] == 1, "polyline must cover [0,1]"
    def __call__(self, x):
        x = Fr(x); P = self.p
        for (x0, y0), (x1, y1) in zip(P, P[1:]):
            if x0 <= x <= x1: return y0 if x1 == x0 else y0 + (y1-y0)*(x-x0)/(x1-x0)
        raise ValueError(x)
    def kinks(self): return [x for x, _ in self.p]

def upper_hull(points):
    P = collapse(points, "max"); H = []          # equal abscissae: only the highest point can be on the upper hull
    for q in P:
        while len(H) >= 2 and (H[-1][0]-H[-2][0])*(q[1]-H[-2][1]) - (H[-1][1]-H[-2][1])*(q[0]-H[-2][0]) >= 0: H.pop()
        H.append(q)
    return Poly(H)

def tent_max_pts(l, u, y):       # y * min(1, x/l, (1-x)/(1-u)) on [0,1]; plateau reaches x=0 if l=0 and x=1 if u=1
    return collapse([(0, 0), (l, y), (u, y), (1, 0)], "max")
def tent_min_pts(l, u, y):       # y * min(x/u, (1-x)/(1-l)), apex at u/(1+u-l)
    assert u > 0 and l < 1
    xs = u/(1+u-l); return collapse([(0, 0), (xs, y*xs/u), (1, 0)], "min")

def add_polys(P1, P2, sign=1):
    xs = sorted(set(P1.kinks()) | set(P2.kinks()))   # both inputs continuous with strictly increasing kinks
    return Poly([(x, P1(x) + sign*P2(x)) for x in xs])

def profiles(chart, box):
    """box: dict of exact rational intervals. Returns (h_out, a_in, b_in) as Poly."""
    if chart == "I":
        (y3l, y3u), (tl, tu), (x1l, x1u), (x3l, x3u) = box["y3"], box["t"], box["x1"], box["x3"]
        y1l, y1u = tl*y3l, tu*y3u
        hout = add_polys(Poly(tent_max_pts(x1l, x1u, y1u)), Poly(tent_max_pts(x3l, x3u, y3u)))
        ain = Poly([(x, -y) for x, y in tent_min_pts(x1l, x1u, y1l)])
        bin_ = Poly(tent_min_pts(x3l, x3u, y3l))
        return hout, ain, bin_
    else:
        (y2l, y2u), (tl, tu), (x2l, x2u), (x3l, x3u) = box["y2"], box["t"], box["x2"], box["x3"]
        y3l, y3u = tl*y2l, tu*y2u
        corners = [(0, 0), (1, 0)] + [(x, y) for x in (x2l, x2u) for y in (y2l, y2u)] + [(x, y) for x in (x3l, x3u) for y in (y3l, y3u)]
        hout = upper_hull(corners)
        def apex(l, u, yl): xs = u/(1+u-l); return (xs, yl*xs/u)
        bin_ = upper_hull([(0, 0), (1, 0), apex(x2l, x2u, y2l), apex(x3l, x3u, y3l)])
        return hout, Poly([(0, 0), (1, 0)]), bin_

def graded_mesh(kinks, n=240, grade=Fr(21, 20), first=Fr(1, 2000)):
    """rational mesh: all kinks, geometric grading towards 0 and 1, uniform interior."""
    pts = set(Fr(k) for k in kinks) | {Fr(0), Fr(1)}
    x = first
    while x < Fr(1, 10):
        pts.add(x); pts.add(1-x); x = x*grade
        x = Fr(x).limit_denominator(10**9)
    for i in range(1, n): pts.add(Fr(i, n))
    xs = sorted(p for p in pts if 0 <= p <= 1)
    out = [xs[0]]
    for v in xs[1:]:
        if v - out[-1] > Fr(1, 10**7) or v in kinks: out.append(v)
    return out

# ---------------- lower bound: Sturm shooting on the minorant -----------------
def minorant(hout, xs):
    return [max(hout(x0), hout(x1)) for x0, x1 in zip(xs, xs[1:])]      # exact rationals (h linear on cells)

def shoot_float(H, xs, E):
    f, df = 0.0, 1.0
    for Hk, x0, x1 in zip(H, xs, xs[1:]):
        d = float(x1-x0); q = math.pi**2/float(Hk)**2 - E
        if q > 0:
            k = math.sqrt(q); c, s = math.cosh(k*d), math.sinh(k*d)
            f, df = c*f + s/k*df, k*s*f + c*df
            m = max(abs(f), abs(df))
            if m > 1e200: f, df = f/m, df/m
        else:
            k = math.sqrt(-q)
            if k*d >= math.pi: return False
            c, s = math.cos(k*d), math.sin(k*d)
            f, df = c*f + (s/k if k > 0 else d)*df, -k*s*f + c*df
        if f <= 0: return False
    return True

def mu1_float(H, xs):
    lo, hi = 0.0, math.pi**2/float(min(H))**2 * 4 + 1e3
    while not shoot_float(H, xs, lo): lo /= 2
    for _ in range(60):
        mid = (lo+hi)/2
        if shoot_float(H, xs, mid): lo = mid
        else: hi = mid
    return lo

def sturm_certify(H, xs, E):
    """rigorous: True only if the solution with energy E has no zero in (0,1] (ball arithmetic)."""
    f, df = arb(0), arb(1); Ea = A(E)
    for Hk, x0, x1 in zip(H, xs, xs[1:]):
        d = A(x1-x0); q = PI2/A(Hk)**2 - Ea
        if q > 0:
            k = q.sqrt(); kd = k*d; c, s = kd.cosh(), kd.sinh()
            f, df = c*f + s/k*df, k*s*f + c*df
        elif q < 0:
            k = (-q).sqrt(); kd = k*d
            if not (kd < PI): return False
            c, s = kd.cos(), kd.sin()
            f, df = c*f + s/k*df, -k*s*f + c*df
        else:
            # q straddles 0: power series of cosh(sqrt(q) d), sinh(sqrt(q) d)/sqrt(q) with rigorous tail
            z = q*d*d; C = arb(0); S = arb(0); term_c = arb(1); term_s = d
            for n in range(0, 30):
                C += term_c; S += term_s
                term_c = term_c*z/((2*n+1)*(2*n+2)); term_s = term_s*z/((2*n+2)*(2*n+3))
            # rigorous tail (referee audit 2026-10-02: the former float tail 2*zr**30 could underflow): for |z| < 1 the omitted
            # terms are bounded by |z|^30/(60! (1-|z|)) (cosh series) and d |z|^30/(61! (1-|z|)) (sinh(sqrt(q)d)/sqrt(q) series),
            # computed in ball arithmetic and widened to their rigorous upper bounds; fail closed unless |z| < 1 is proved.
            zr = arb(abs(z).upper())
            if not (zr < 1): return False
            fac60 = arb.fac_ui(60)
            C += arb(0, (zr**30/fac60/(1 - zr)).upper()); S += arb(0, (d*zr**30/(fac60*61)/(1 - zr)).upper())
            dS = q*S          # derivative of C is q*S
            f, df = C*f + S*df, dS*f + C*df
        if not (f > 0): return False
        m = max(abs(f).upper(), abs(df).upper())
        if m > arb(2)**600: f, df = f/m, df/m
    return True

# ---------------- upper bound: weighted 1-D Ritz on the inner domain -----------------
def ritz_mesh_values(ain, bin_, xs):
    a = [ain(x) for x in xs]; b = [bin_(x) for x in xs]; return a, b

def assemble_float(a, b, xs):
    n = len(xs); K = np.zeros((n, n)); M = np.zeros((n, n))
    for e in range(n-1):
        d = float(xs[e+1]-xs[e]); a0, a1 = float(a[e]), float(a[e+1]); b0, b1 = float(b[e]), float(b[e+1])
        h0, h1 = b0-a0, b1-a1; ap, bp = (a1-a0)/d, (b1-b0)/d; hp = bp-ap
        coef = 0.5 + 0.5*ap*bp + hp*hp*float(1/6 + 1/(4*math.pi**2))
        gp = np.array([-0.906179845938664, -0.538469310105683, 0, 0.538469310105683, 0.906179845938664])
        gw = np.array([0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189])
        for g, w in zip(gp, gw):
            t = 0.5*(1+g); h = h0 + (h1-h0)*t; W = 0.5*d*w
            phi = np.array([1-t, t]); dphi = np.array([-1/d, 1/d])
            for i in range(2):
                for j in range(2):
                    K[e+i, e+j] += W*(dphi[i]*dphi[j]*h/2 + (phi[i]*dphi[j]+dphi[i]*phi[j])*hp/4 + phi[i]*phi[j]*math.pi**2/h*coef)
                    M[e+i, e+j] += W*phi[i]*phi[j]*h/2
    I = slice(1, n-1)
    w, V = eigh(K[I, I], M[I, I], subset_by_index=[0, 1])
    F = np.zeros((n, 2)); F[1:n-1, :] = V
    return w, F

def elem_integrals(f, g, a0, a1, b0, b1, d):
    """exact (ball) element integrals of N and D for P1 f=(f0,f1), g=(g0,g1); a,b linear; d = element length (arb)."""
    f0, f1 = f; g0, g1 = g
    h0, h1 = b0-a0, b1-a1
    ap, bp = (a1-a0)/d, (b1-b0)/d; hp = bp-ap
    coef = arb(1)/2 + ap*bp/2 + hp*hp*C6
    fp, gp_ = (f1-f0)/d, (g1-g0)/d
    I_h = (h0+h1)*d/2                                   # int h
    I_f = (f0+f1)*d/2; I_g = (g0+g1)*d/2
    N = fp*gp_*I_h/2 + hp/4*(gp_*I_f + fp*I_g)
    # int f g h (three linear factors) and int f g / h
    fgh = d*((f0*g0*h0 + f1*g1*h1)/4 + (f1*g0*h0 + f0*g1*h0 + f0*g0*h1)/12 + (f1*g1*h0 + f1*g0*h1 + f0*g1*h1)/12)
    Dm = fgh/2
    return N, Dm, coef, (f0, f1, g0, g1, h0, h1)

def int_fg_over_h(f0, f1, g0, g1, h0, h1, d, h0_exact, h1_exact):
    if h0_exact == 0 or h1_exact == 0:
        # f, g vanish where h does (Dirichlet ends): f/h constant on the element
        if h0_exact == 0:  return (f1/h1)*(g0+g1)*d/2 if True else None
        else:              return (f0/h0)*(g0+g1)*d/2
    if h0_exact == h1_exact:
        return d/h0*(f0*g0/3 + (f0*g1+f1*g0)/6 + f1*g1/3)
    # t = h(x) linear; f = af + bf t, g = ag + bg t; dx = dt/h', h' = (h1-h0)/d
    hp = (h1-h0)/d
    bf = (f1-f0)/(h1-h0); af = f0 - bf*h0
    bg = (g1-g0)/(h1-h0); ag = g0 - bg*h0
    F = lambda t: af*ag*t.log() + (af*bg+ag*bf)*t + bf*bg*t*t/2
    return (F(h1) - F(h0))/hp

def ritz_certify(ain, bin_, xs, F):
    a = [A(ain(x)) for x in xs]; b = [A(bin_(x)) for x in xs]
    hex_ = [bin_(x)-ain(x) for x in xs]
    Fv = [[arb(float(F[k, j])) for j in range(2)] for k in range(len(xs))]
    Nm = [[arb(0)]*2 for _ in range(2)]; Dm = [[arb(0)]*2 for _ in range(2)]
    Nm = [[arb(0), arb(0)], [arb(0), arb(0)]]; Dm = [[arb(0), arb(0)], [arb(0), arb(0)]]
    for e in range(len(xs)-1):
        d = A(xs[e+1]-xs[e])
        for i in range(2):
            for j in range(i, 2):
                f = (Fv[e][i], Fv[e+1][i]); g = (Fv[e][j], Fv[e+1][j])
                N, D, coef, (f0, f1, g0, g1, h0, h1) = elem_integrals(f, g, a[e], a[e+1], b[e], b[e+1], d)
                N = N + PI2*coef*int_fg_over_h(f0, f1, g0, g1, h0, h1, d, hex_[e], hex_[e+1])
                Nm[i][j] += N; Dm[i][j] += D
    Nm[1][0] = Nm[0][1]; Dm[1][0] = Dm[0][1]
    detD = Dm[0][0]*Dm[1][1] - Dm[0][1]**2
    if not (Dm[0][0] > 0 and detD > 0): return None
    B = Nm[0][0]*Dm[1][1] + Nm[1][1]*Dm[0][0] - 2*Nm[0][1]*Dm[0][1]
    detN = Nm[0][0]*Nm[1][1] - Nm[0][1]**2
    disc = B*B - 4*detD*detN
    if not (disc > 0): return None
    return (B + disc.sqrt())/(2*detD)

def certify_box(chart, box, n=240, verbose=False):
    hout, ain, bin_ = profiles(chart, box)
    kinks = sorted(set(hout.kinks()) | set(ain.kinks()) | set(bin_.kinks()))
    xs = graded_mesh(kinks, n)
    H = minorant(hout, xs)
    mu = mu1_float(H, xs); E = Fr(mu*(1 - 2e-6)).limit_denominator(10**12)
    okE = sturm_certify(H, xs, E)
    a, b = ritz_mesh_values(ain, bin_, xs)
    w, F = assemble_float(a, b, xs)
    R2 = ritz_certify(ain, bin_, xs, F)
    ok = okE and R2 is not None and (R2 < A(Fr(5, 2))*A(E))
    if verbose: print(f"  mesh {len(xs)}, mu1_float {mu:.6f}, E ok {okE}, Ritz float {w[1]:.6f}, R2 {R2}, ratio {float(R2.mid())/float(E) if R2 is not None else None}")
    return ok, E, R2

if __name__ == "__main__":
    # point boxes (degenerate intervals) of the prototype's test shapes, chart I: (x1,y1,x3,y3)
    for (x1, y1, x3, y3) in [(Fr(4, 10), Fr(1, 10), Fr(5, 10), Fr(3, 10)), (Fr(3, 10), Fr(2, 10), Fr(45, 100), Fr(25, 100))]:
        box = {"y3": (y3, y3), "t": (y1/y3, y1/y3), "x1": (x1, x1), "x3": (x3, x3)}
        print("chart I point", (float(x1), float(y1), float(x3), float(y3)), certify_box("I", box, verbose=True)[0])
