"""Rigorous one-dimensional (Born-Oppenheimer) LOWER bound for lambda_1 of every dart of a box (LEMMAS.md, Lemma 11.1).

Lemma 11.1 (1-D reduction with the longest slice component).  Let Omega be open and bounded or not, u a unit vector,
x = <., u>, and l(x) the length of the LONGEST component of the slice Omega cap {x = const} (l = 0 if empty).  Then
lambda_1(Omega) >= mu_1(-d^2/dx^2 + V) on the projection interval (Dirichlet), for every V with 0 <= V <= pi^2/l^2.
(For phi in C_c^oo(Omega): v(x) = ||phi(x,.)|| has |v'| <= ||phi_x(x,.)|| (compact support: no boundary terms), and on every
slice component J (an interval with phi = 0 near its ends) int_J |phi_y|^2 >= (pi/|J|)^2 int_J phi^2 >= (pi/l)^2 int_J phi^2.)

Outer domains.  For a box of darts we use, for every slab S_k = {x_k <= x <= x_{k+1}}, an upper bound H_k of l over S_k valid for
every dart of the box:  a dart is Q = int(ABD u CBD u BD) with ABD subset conv(SA), CBD subset conv(SC) (SA = {B, D, corners of
the A-box}, SC likewise).  The y-extent of conv(S) cap S_k is computed from the points of S in the slab and the intersections of
ALL segments between points of S with the two slab lines (a superset of the vertices of conv(S) cap S_k), in ball arithmetic.
If the slab meets the projection [x_B, x_D] of the segment BD:  H_k = ext_A + ext_C;  otherwise slices do not meet the open
segment BD, the A- and C-parts of a slice are different components and H_k = max(ext_A, ext_C).
Spikes (vertical slicing, A = a (1, t), t <= t_hi, any a > 0): the A-part has slice length <= t_hi min(x, 1) for x >= 0
(the triangle BDA lies below the ray BA and, beyond x = 1, between DA and BA, at vertical distance <= t_hi), so ext_A(S_k) =
t_hi min(x_{k+1}, 1); beyond x_1 = max(1, x-range of conv(SC)) the only slices are spike slices: V >= V_oo := pi^2/t_hi^2 on
[x_1, oo) (tail).

Sturm certificate (ball arithmetic): the solution of -v'' + (V_k - E) v = 0, v(x_0) = 0, v'(x_0) = 1, has no zero in (x_0, x_K]
(hyperbolic cells: positivity at both cell ends; trigonometric cells: phase condition), and with a tail additionally
v'(x_K) + sqrt(V_oo - E) v(x_K) > 0 (then v stays positive on [x_K, oo) where it is v cosh + v'/kappa sinh).  By Sturm
comparison this gives mu_1 > E, hence lambda_1(Q) > E for every dart of the box."""
import math, numpy as np
from flint import arb, ctx
from nanguard import fin, afin, fdown, fup, amax, amin, upper_min, lower_max
ctx.prec = 106
PI = arb.pi()

def A(x): return x if isinstance(x, arb) else arb(x)

# ----------------------------------------------------------------------------------------- slab extents (rigorous)
def slab_extent(S, u, a, b):
    """upper bound of the y-extent (y = <p, u_perp>) of conv(S) cap {a <= <p,u> <= b}; S = list of (x, y) floats (exact points);
    u = (ux, uy) arb; a, b floats.  Returns arb upper bound (0 if the slab certainly misses conv(S))."""
    ux, uy = u
    X = [afin(A(p[0])*ux + A(p[1])*uy) for p in S]
    Y = [afin(-A(p[0])*uy + A(p[1])*ux) for p in S]
    aa, bb = A(a), A(b)
    lo = None; hi = None
    def take(y):
        nonlocal lo, hi
        afin(y)                                                     # a NaN would be dropped by the comparisons below
        l, h = y.lower(), y.upper()
        lo = l if lo is None or l < lo else lo
        hi = h if hi is None or h > hi else hi
    n = len(S)
    for i in range(n):
        if X[i].upper() >= a and X[i].lower() <= b: take(Y[i])
    for i in range(n):
        for j in range(i + 1, n):
            xlo = amin(X[i], X[j]); xhi = amax(X[i], X[j])                    # exact arb points
            dx = X[j] - X[i]
            for xi in (aa, bb):
                if not (xlo <= xi.upper() and xi.lower() <= xhi): continue
                if dx.contains(0):
                    take(Y[i].union(Y[j])); continue
                t = (xi - X[i])/dx
                tl = lower_max(t, arb(0)); th = upper_min(t, arb(1))              # t in [0, 1] on the segment
                if tl > th: continue
                tt = arb(tl).union(arb(th))
                take(Y[i] + tt*(Y[j] - Y[i]))
    if lo is None: return arb(0)
    return arb(hi) - arb(lo)

# ----------------------------------------------------------------------------------------- Sturm shooting
def no_zero_float(nodes, H, E, tail=None):   # RIGOR-FUNCTION: heuristic: float search for a candidate E (verified by no_zero_arb)
    v, dv = 0.0, 1.0
    for k in range(len(H)):
        L = nodes[k + 1] - nodes[k]
        q = (math.pi**2/H[k]**2 if H[k] > 0 else 1e300) - E
        if q > 0:
            kap = math.sqrt(q); th = math.tanh(kap*L)
            vn = v + dv*th/kap; dvn = v*kap*th + dv
        elif q == 0:
            vn = v + dv*L; dvn = dv
        else:
            om = math.sqrt(-q); c = math.cos(om*L); s = math.sin(om*L)
            if v == 0.0:
                if om*L >= math.pi: return False
            else:
                if math.pi/2 + math.atan2(dv/om, v) - om*L <= 0: return False
            vn = v*c + dv*s/om; dvn = -v*om*s + dv*c
        if not vn > 0: return False
        m = max(abs(vn), abs(dvn)); v, dv = vn/m, dvn/m
    if tail is not None:
        if tail <= E: return False
        return dv + math.sqrt(tail - E)*v > 0
    return True

def no_zero_arb(nodes, Hup, E, tail=None):
    """rigorous version: Hup = list of arb upper bounds of the slice lengths (V_k = pi^2/Hup_k^2 is a minorant)"""
    E = A(E); v = arb(0); dv = arb(1); started = False
    for k in range(len(Hup)):
        L = A(nodes[k + 1]) - A(nodes[k])
        h = Hup[k]
        if not (h.upper() > 0):
            Vk = arb('1e300')                                   # certainly empty slab (no slice): any minorant value is allowed
        else:
            hu = arb(h.upper()); Vk = PI*PI/(hu*hu)
        q = Vk - E
        if q.lower() > 0:
            kap = q.sqrt(); th = (kap*L).tanh()
            vn = v + dv*th/kap; dvn = v*kap*th + dv              # scaled by 1/cosh > 0
            if not (vn.lower() > 0): return False
            if started and not (v.lower() > 0): return False
            if not started and not (dv.lower() > 0): return False
        elif q.upper() < 0:
            om = (-q).sqrt(); c = (om*L).cos(); s = (om*L).sin()
            if not started:
                if not (dv.lower() > 0): return False
                if not ((PI - om*L).lower() > 0): return False
            else:
                if not (v.lower() > 0): return False
                ph = arb.atan2(dv/om, v)
                if not ((PI/2 + ph - om*L).lower() > 0): return False
            vn = v*c + dv*s/om; dvn = -v*om*s + dv*c
            if not (vn.lower() > 0): return False
        else:
            return False
        m = max(abs(float(vn.mid())), abs(float(dvn.mid())), 1e-300)   # RIGOR: heuristic: free positive rescaling of (v, v'): every test is scale invariant
        v = vn/m; dv = dvn/m; started = True
    if tail is not None:
        T = A(tail) - E
        if not (T.lower() > 0): return False
        return bool((dv + T.sqrt()*v).lower() > 0)
    return True

def mu_search(nodes, Hf, tail=None, iters=60):
    """float bisection for the largest E with no_zero_float True"""
    lo = 0.0
    if not no_zero_float(nodes, Hf, lo, tail): return 0.0
    hi = 1.0
    while no_zero_float(nodes, Hf, hi, tail) and hi < 1e15: lo = hi; hi *= 2
    for _ in range(iters):
        mid = 0.5*(lo + hi)
        if no_zero_float(nodes, Hf, mid, tail): lo = mid
        else: hi = mid
    return lo

def certify_E(nodes, Hup, tail=None):
    """largest certified E (rigorous): float search on the minorant, then ball-arithmetic verification with back-off"""
    Hf = [float(h.upper()) if h.upper() > 0 else 0.0 for h in Hup]   # RIGOR: heuristic: float search data; the returned E is verified by no_zero_arb
    Ef = mu_search(nodes, Hf, tail)
    for back in (1e-9, 1e-6, 1e-4, 1e-3, 1e-2, 3e-2):
        E = Ef*(1 - back)
        if E <= 0: break
        if no_zero_arb(nodes, Hup, E, tail): return E
    return 0.0

# ----------------------------------------------------------------------------------------- nodes
def graded_nodes(x0, x1, kinks, K=160, rmin=2e-4, ratio=1.25):   # RIGOR-FUNCTION: heuristic: node placement (x0 and x1 are kept exactly)
    """nodes on [x0, x1]: uniform interior spacing L/K, geometric grading (ratio) towards both ends down to rmin*L, plus kinks"""
    L = x1 - x0
    pts = set([x0, x1])
    h = L/K
    for k in range(1, K): pts.add(x0 + k*h)
    for end, sgn in ((x0, 1), (x1, -1)):
        d = h
        while d > rmin*L:
            d /= ratio; pts.add(end + sgn*d)
    for z in kinks:
        if x0 < z < x1: pts.add(float(z))
    return sorted(pts)

# ----------------------------------------------------------------------------------------- dart boxes
def corners2(lo, hi):
    return [(lo[0], lo[1]), (hi[0], lo[1]), (lo[0], hi[1]), (hi[0], hi[1])]

def bo_box(Abox, Cbox, theta, K=160):
    """rigorous lower bound for lambda_1 over all darts Q(A, C), A in Abox = ((ax_lo, ay_lo), (ax_hi, ay_hi)), C in Cbox,
    slicing axis u = (cos theta, sin theta) (theta float; u is an arb unit vector).  Returns (E, info)."""
    th = A(theta); u = (th.cos(), th.sin())
    B = (0.0, 0.0); D = (1.0, 0.0)
    SA = [B, D] + corners2(*Abox); SC = [B, D] + corners2(*Cbox)
    pr = lambda p: (A(p[0])*u[0] + A(p[1])*u[1])
    xs = [pr(p) for p in SA + SC]
    x0 = fdown(amin(*xs)); x1 = fup(amax(*xs))               # directed projection range (left down, right up)
    xB, xD = pr(B), pr(D)
    bd_lo = amin(xB, xD); bd_hi = amax(xB, xD)              # exact arb points
    kinks = [float(x.mid()) for x in xs]   # RIGOR: heuristic: node placement
    nodes = graded_nodes(x0, x1, kinks, K)
    Hup = []
    for k in range(len(nodes) - 1):
        a, b = nodes[k], nodes[k + 1]
        eA = slab_extent(SA, u, a, b); eC = slab_extent(SC, u, a, b)
        if b >= bd_lo and a <= bd_hi: Hup.append(eA + eC)
        else: Hup.append(amax(eA, eC))
    E = certify_E(nodes, Hup)
    return E, dict(theta=theta, n=len(nodes))

def bo_spike(t_hi, Cbox, K=160):
    """rigorous lower bound for lambda_1 over all darts with A = a (1, t), 0 < t <= t_hi, any a > 0 (spikes, collapses), C in Cbox;
    vertical slicing (u = (1, 0)) with the spike tail.  Returns (E, info)."""
    u = (arb(1), arb(0))
    SC = [(0.0, 0.0), (1.0, 0.0)] + corners2(*Cbox)
    xs = [A(p[0]) for p in SC]
    x0 = fdown(amin(arb(0), *xs)); x1 = fup(amax(arb(1), *xs))   # directed; the spike part lives over [0, 1]
    nodes = graded_nodes(x0, x1, [p[0] for p in SC] + [0.0, 1.0], K)
    th = A(t_hi)
    Hup = []
    for k in range(len(nodes) - 1):
        a, b = nodes[k], nodes[k + 1]
        eC = slab_extent(SC, u, a, b)
        eA = th*min(b, 1.0) if b > 0 else arb(0)   # RIGOR: exact: min of two floats
        if b >= 0 and a <= 1: Hup.append(eC + eA)           # the slab meets the projection [0, 1] of BD
        else: Hup.append(amax(eC, eA))
    tail = fdown(PI*PI/(th*th))                              # directed: the tail potential is a lower bound
    E = certify_E(nodes, Hup, tail=tail)
    return E, dict(n=len(nodes), tail=tail)
