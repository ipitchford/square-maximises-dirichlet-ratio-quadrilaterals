"""General rigorous 1-D (Born-Oppenheimer) lower bound for lambda_1 over a FAMILY of darts described by arm specs, with tails
for unbounded arms (LEMMAS.md, Lemma 11.1 and Lemma 11.2).  Dart: Q = int(T_A u T_C u BD), T_A = tri(B, D, A), T_C = tri(B, D, C),
B = (0,0), D = (1,0), A in the upper and C in the lower half plane.

Arm specs (side = +1 for A, -1 for C; direction e(b) = (cos b, side*sin b), b in (0, pi) measured from BD):
  ('polar', b1, b2, r1, r2) : X = rho e(b), b in [b1, b2], rho in [r1, r2]            (compact)
  ('far',   b1, b2, R)      : X = rho e(b), b in [b1, b2], rho >= R  (R >= 0; any length) (unbounded)
Facts used (proved in LEMMAS 11.2):
  (S) T_X is contained in the closed half-strip S(b) = {B + s e(b) + tau (D - B): s >= 0, 0 <= tau <= 1}.
  (C) A slicing line {<p, u> = x} that does not meet the segment BD meets T_X in a segment between the sides BX and DX, of length
      <= w(b)/|<e(b), u>| with w(b) = sin b (the width of S(b)), whenever <e(b), u> != 0.
  (U) For every slab S_k the union over the family of T_X cap S_k lies in conv(P) cap S_k for the point set P below; so the
      y-extent of conv(P) cap S_k bounds every slice of every T_X in the slab (bo1d.slab_extent):
        polar: P = {B, D, r1 e(b1), r1 e(b2), r2 e(b1), r2 e(b2), (r2/cos((b2-b1)/2)) e((b1+b2)/2)};
        far  : only in the near region {<p,u> <= x_cut} (arm going to +x) resp. {>= x_cut} (to -x):
               P = {B, D, the intersections of the boundary rays B + s e(b_i), D + s e(b_i) (i = 1, 2) with the cut line}.
  Per slab: H_k = ext_A + ext_C if the slab meets the projection of BD (slices are connected through BD), else max(ext_A, ext_C);
  ext_X = min(union extent, per-dart chord bound (C) if the slab misses the projection of BD).
  Tails: beyond x_cut a far arm going to +x contributes only slices of length <= its chord bound: V >= pi^2/c^2 there.
1-D problem on [x_L, x_R] (Dirichlet at a finite end, exponential tail condition at an unbounded end), certified by Sturm
shooting in ball arithmetic.  Result: lambda_1(Q) >= E for every dart of the family."""
import math, numpy as np
from flint import arb, ctx
from nanguard import fin, afin, fdown, fup, amax, amin, upper_min, cuts
from bo1d import slab_extent, graded_nodes, A
ctx.prec = 106
PI = arb.pi()

def e_arb(b, side):
    b = A(b); return (b.cos(), side*b.sin())
def b_ball(b1, b2):
    return arb(b1).union(arb(b2))

def far_of(arm):
    """the far arm {b in [b1, b2], rho >= r1} of a polar arm: its family CONTAINS the polar one, so a lower bound of
    lambda_1 over a family with this arm (Lemma 11.3) is a lower bound for the polar family too (LEMMAS 11.3, remark (4))"""
    assert arm.kind == 'polar'
    return Arm(('far', arm.b1, arm.b2, arm.r1), arm.side)

def arm_variant(armA, armC, mode):
    """arms used by a BO2 record with key 'arms' = mode (absent: 'polar')"""
    if mode in (None, 'polar'): return armA, armC
    if mode == 'farA': return far_of(armA), armC
    if mode == 'farC': return armA, far_of(armC)
    if mode == 'farAC': return far_of(armA), far_of(armC)
    raise ValueError('unknown arm mode ' + str(mode))

class Arm:
    def __init__(self, spec, side):
        self.kind = spec[0]; self.side = side
        self.b1, self.b2 = float(spec[1]), float(spec[2])   # RIGOR: exact: arm specs are floats (float() of a float is exact)
        if self.kind == 'polar': self.r1, self.r2 = float(spec[3]), float(spec[4])   # RIGOR: exact: arm specs are floats
        elif self.kind == 'far': self.R = float(spec[3])   # RIGOR: exact: arm specs are floats
        else: raise ValueError(spec)
    # direction sign along u: +1 (arm goes to +x), -1 (to -x), 0 (undetermined)
    def dir_sign(self, u):
        bb = b_ball(self.b1, self.b2)
        d = bb.cos()*u[0] + self.side*bb.sin()*u[1]
        if d.lower() > 0: return 1
        if d.upper() < 0: return -1
        return 0
    def chord_bound(self, u):
        """per-dart chord bound (C), arb upper bound over b in [b1, b2]; None if <e, u> may vanish"""
        bb = b_ball(self.b1, self.b2)
        d = bb.cos()*u[0] + self.side*bb.sin()*u[1]
        if d.contains(0): return None
        # w(b) = sin b; bound sin(b)/|d(b)| by subdividing the angle interval (the ball over the whole interval is coarse)
        n = 8; cs = []
        xs = cuts(self.b1, self.b2, n)                       # cells cover [b1, b2] exactly
        for i in range(n):
            bi = b_ball(xs[i], xs[i + 1])
            di = bi.cos()*u[0] + self.side*bi.sin()*u[1]
            if di.contains(0): return None
            cs.append(afin(bi.sin()/abs(di)))
        return amax(*cs)
    def outer_points(self, u=None, x_cut=None):
        B = (arb(0), arb(0)); D = (arb(1), arb(0))
        if self.kind == 'polar':
            pts = [B, D]
            for r in (self.r1, self.r2):
                for b in (self.b1, self.b2):
                    e = e_arb(b, self.side); pts.append((arb(r)*e[0], arb(r)*e[1]))
            half = (A(self.b2) - A(self.b1))/2
            if not (half.upper() < PI/2): raise ValueError('angle interval too wide')
            e = e_arb((A(self.b1) + A(self.b2))/2, self.side); rr = A(self.r2)/half.cos()      # exact mid-angle (ball)
            pts.append((rr*e[0], rr*e[1]))
            return pts
        # far: near-region points at the cut line <p, u> = x_cut
        pts = [B, D]
        X = A(x_cut)
        for b in (self.b1, self.b2):
            e = e_arb(b, self.side); de = e[0]*u[0] + e[1]*u[1]
            for P0 in (B, D):
                lam = (X - (P0[0]*u[0] + P0[1]*u[1]))/de
                if not (lam.lower() >= 0): raise ValueError('cut line behind the arm base')
                pts.append((P0[0] + lam*e[0], P0[1] + lam*e[1]))
        return pts
    def inner_vertex(self):
        """a float point X_in with tri(B,D,X_in) inside tri(B,D,X) for every X of the spec (angles at B and D not larger), or None"""
        bmin = A(self.b1)
        rho = A(self.r1) if self.kind == 'polar' else A(self.R)
        if rho.upper() <= 0: return None
        bb = b_ball(self.b1, self.b2)
        # angle at D: atan2(rho sin b, 1 - rho cos b), increasing in rho (LEMMAS 11.2); lower bound over b in [b1, b2] at rho_min
        n = 8; dmin = None
        xs = cuts(self.b1, self.b2, n)                       # cells cover [b1, b2] exactly
        for i in range(n):
            bi = b_ball(xs[i], xs[i + 1])
            ang = afin(arb.atan2(rho*bi.sin(), 1 - rho*bi.cos()))
            dmin = ang.lower() if dmin is None or ang.lower() < dmin else dmin     # exact arb points
        dmin = arb(dmin)
        if not (dmin.lower() > 0 and bmin.lower() > 0 and (bmin + dmin).upper() < PI): return None
        # intersection of the ray from B at angle bmin and the ray from D at angle dmin (measured at D from DB)
        # X = t (cos bmin, sin bmin) with tan(dmin) = t sin bmin / (1 - t cos bmin)  ->  t = sin dmin / sin(bmin + dmin)
        t = dmin.sin()/(bmin + dmin).sin()
        xf = float((t*bmin.cos()).mid()); yf = float((t*bmin.sin()).mid())   # RIGOR: heuristic: a candidate inner vertex, verified rigorously below
        # shrink towards the midpoint of BD and verify the angles rigorously
        for shrink in (1e-12, 1e-9, 1e-6, 1e-4):
            px = 0.5 + (xf - 0.5)*(1 - shrink); py = yf*(1 - shrink)
            aB = arb.atan2(arb(py), arb(px)); aD = arb.atan2(arb(py), 1 - arb(px))
            if aB.upper() <= bmin.lower() and aD.upper() <= dmin.lower() and py > 0:
                return (px, self.side*py)
        return None

def _no_zero_float(nodes, H, E, kL=None, kR=None):   # RIGOR-FUNCTION: heuristic: float search for a candidate E (verified by _no_zero_arb)
    if kL is None: v, dv = 0.0, 1.0
    else:
        if kL <= E: return False
        v, dv = 1.0, math.sqrt(kL - E)
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
    if kR is not None:
        if kR <= E: return False
        return dv + math.sqrt(kR - E)*v > 0
    return True

def _no_zero_arb(nodes, Hup, E, kL=None, kR=None):
    E = A(E)
    if kL is None: v = arb(0); dv = arb(1); started = False
    else:
        T = A(kL) - E
        if not (T.lower() > 0): return False
        v = arb(1); dv = T.sqrt(); started = True
    for k in range(len(Hup)):
        L = A(nodes[k + 1]) - A(nodes[k])
        h = Hup[k]
        if not (h.upper() > 0): Vk = arb('1e300')
        else:
            hu = arb(h.upper()); Vk = PI*PI/(hu*hu)
        q = Vk - E
        if q.lower() > 0:
            kap = q.sqrt(); th = (kap*L).tanh()
            if started and not (v.lower() > 0): return False
            if not started and not (dv.lower() > 0): return False
            vn = v + dv*th/kap; dvn = v*kap*th + dv
            if not (vn.lower() > 0): return False
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
    if kR is not None:
        T = A(kR) - E
        if not (T.lower() > 0): return False
        return bool((dv + T.sqrt()*v).lower() > 0)
    return True

def certify_E2(nodes, Hup, kL=None, kR=None):   # RIGOR-FUNCTION: heuristic: float search; the returned E is verified by _no_zero_arb (rigorous)
    Hf = [float(h.upper()) if h.upper() > 0 else 0.0 for h in Hup]
    kLf = float(kL.lower()) if kL is not None else None; kRf = float(kR.lower()) if kR is not None else None
    lo = 0.0
    if not _no_zero_float(nodes, Hf, lo, kLf, kRf): return 0.0
    hi = 1.0
    while _no_zero_float(nodes, Hf, hi, kLf, kRf) and hi < 1e15: lo = hi; hi *= 2
    for _ in range(60):
        mid = 0.5*(lo + hi)
        if _no_zero_float(nodes, Hf, mid, kLf, kRf): lo = mid
        else: hi = mid
    for back in (1e-9, 1e-6, 1e-4, 1e-3, 1e-2, 3e-2):
        E = lo*(1 - back)
        if E <= 0: break
        if _no_zero_arb(nodes, Hup, E, kL, kR): return E
    return 0.0

def proj_range(pts, u):
    """floats (lo, hi) with lo <= <p, u> <= hi for every point p of the (ball) point list, and the projection balls: directed
    (the left extremum rounded down, the right one up)"""
    prc = [afin(p[0]*u[0] + p[1]*u[1]) for p in pts]
    return fdown(amin(*prc)), fup(amax(*prc)), prc

def bo_family(armA, armC, theta, K=160):
    """rigorous lower bound E for lambda_1 of every dart (A, C) with A in armA (an Arm, side +1, or None: no A-triangle, e.g. the
    collar's CBD-only P2 domain is irrelevant here) and C in armC (side -1); slicing axis u = (cos theta, sin theta).
    Returns (E, info); E = 0.0 if nothing could be certified (e.g. an unbounded arm not transversal to the slices)."""
    th = A(theta); u = (th.cos(), th.sin())
    arms = [a for a in (armA, armC) if a is not None]
    prB = arb(0); prD = u[0]                                  # projections of B = (0,0), D = (1,0)
    bd_lo = amin(prB, prD); bd_hi = amax(prB, prD)                    # exact arb points
    # compact extent and far-arm directions
    cpts = [(arb(0), arb(0)), (arb(1), arb(0))]
    signs = []
    for a in arms:
        if a.kind == 'polar': cpts += a.outer_points()
        else:
            s = a.dir_sign(u)
            if s == 0: return 0.0, dict(why='far arm not transversal to the slices')
            signs.append((a, s))
    xc_lo, xc_hi, prc = proj_range(cpts, u)
    # cut positions: just outside the compact projection range (exact floats)
    span = max(xc_hi - xc_lo, 1e-6)   # RIGOR: heuristic: placement of the tail cuts (any position beyond the compact range is valid)
    x_cut_R = xc_hi + 0.05*span; x_cut_L = xc_lo - 0.05*span
    near = {}; tails = {+1: [], -1: []}
    for a, s in signs:
        xc = x_cut_R if s > 0 else x_cut_L
        try: near[id(a)] = a.outer_points(u, xc)
        except ValueError as e: return 0.0, dict(why=str(e))
        cb = a.chord_bound(u)
        if cb is None: return 0.0, dict(why='chord bound undefined')
        tails[s].append(cb)
    x_lo = x_cut_L if tails[-1] else xc_lo
    x_hi = x_cut_R if tails[+1] else xc_hi
    kinks = [float(x.mid()) for x in prc] + [x_cut_L, x_cut_R]   # RIGOR: heuristic: node placement (any nodes are valid; the cuts are exact floats)
    nodes = graded_nodes(x_lo, x_hi, kinks, K)
    # per-arm point sets for the union extents
    ptsA = None; ptsC = None
    def pts_for(a):
        if a is None: return None
        return a.outer_points() if a.kind == 'polar' else near[id(a)]
    ptsA = pts_for(armA); ptsC = pts_for(armC)
    cbA = armA.chord_bound(u) if armA is not None else None
    cbC = armC.chord_bound(u) if armC is not None else None
    Hup = []
    for k in range(len(nodes) - 1):
        a_, b_ = nodes[k], nodes[k + 1]
        meets_bd = (b_ >= bd_lo and a_ <= bd_hi)
        exts = []
        for arm, pts, cb in ((armA, ptsA, cbA), (armC, ptsC, cbC)):
            if arm is None: exts.append(arb(0)); continue
            e = slab_extent(pts, u, a_, b_)
            if not meets_bd and cb is not None: e = upper_min(e, cb)         # both bound the slice: the smaller upper end
            exts.append(e)
        afin(exts[0]); afin(exts[1])
        if meets_bd: Hup.append(exts[0] + exts[1])
        else: Hup.append(amax(exts[0], exts[1]))
    cL = amax(*tails[-1]) if tails[-1] else None          # upper bounds of the chord bounds (balls, no float conversion)
    cR = amax(*tails[+1]) if tails[+1] else None
    kL = (PI*PI/(cL*cL)) if cL is not None else None
    kR = (PI*PI/(cR*cR)) if cR is not None else None
    E = certify_E2(nodes, Hup, kL, kR)
    return E, dict(theta=theta, n=len(nodes), kL=float(kL.mid()) if kL is not None else None, kR=float(kR.mid()) if kR is not None else None)   # RIGOR: info: logged
