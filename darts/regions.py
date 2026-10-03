"""Region bookkeeping for the dart cover (BD chart, WLOG b_A <= b_C by the mirror symmetry in BD).
All predicates on boxes are RIGOROUS (interval arithmetic in arb): a box is discarded only if the predicate holds on the
whole box.  Regions (R = 3, t1 = 0.02, t2 = 0.3, X0 = 1.5):
  MAIN  : |A| <= R, |C| <= R, tan b_A >= t2 (A-triangle not spiky)                          -> BD chart, 4 parameters
  SHORT : |A| <= R, |C| <= R, t1 <= tan b_A <= t2, a_x <= X0                                 -> BD chart, 4 parameters
  SPIKE : t1 <= tan b_A <= t2, a_x >= X0 (any |A|), |C| <= R                                 -> SPIKE2 chart (t, cx, cy)
  COLLAR: tan b_A <= t1 (any |A|), |C| <= R                                                   -> SPIKE chart at t = t1 (cx, cy)
Lemmas: Prop 5.1 (beta <= 0.0339), Lemma 3.1 (collapse collars), Lemma 5.1 (tan b_A + tan b_C <= 0.1).
NOT covered here: |C| > R (or |A| > R with tan b_A >= t2): long arms (Lemma 9.3 charts designed, not run)."""
import numpy as np
from flint import arb
from nanguard import fin, afin, fdown, fup, exact_box, dec, decq, amin, amax

R_MAX = 3.0; T1 = 0.02; T2 = 0.3; X0 = 1.5

def iv(lo, hi): return arb(lo).union(arb(hi))
def corners(p0, r):
    import itertools
    for s in itertools.product((-1, 1), repeat=4):
        yield [p0[i] + s[i]*r[i] for i in range(4)]

def fD_bounds(p0, r):
    """min/max of f_D = (1 - cx) ay + cy (ax - 1) over the box (multilinear -> extremes at vertices), rigorous: exact arb points
    lo <= f_D <= hi on the box (no float conversion)"""
    exact_box(p0, r)
    lo = None; hi = None
    for c in corners(p0, r):
        v = afin((1 - arb(c[2]))*arb(c[1]) + arb(c[3])*(arb(c[0]) - 1))
        l, h = v.lower(), v.upper()
        lo = l if lo is None or l < lo else lo
        hi = h if hi is None or h > hi else hi
    return lo, hi

def beta_max(p0, r):
    """an upper bound (exact arb point, or None if undecided) of beta over the box"""
    exact_box(p0, r)
    ax, ay, cx, cy = [iv(p0[i] - r[i], p0[i] + r[i]) for i in range(4)]
    n1 = ay - cy; n2 = cx - ax
    num = abs(n1*(1 - ax) + n2*(0 - ay)); den = abs(n1*(0 - ax) + n2*(0 - ay))
    if not (den.lower() > 0): return None
    return afin(num/den).upper()

def tanb_bounds(x_lo, x_hi, y_lo, y_hi):
    """bounds (exact arb points) of tan b = |y|/x for x > 0 (y_lo >= 0); None if x_lo <= 0"""
    if x_lo <= 0: return None
    lo = afin(arb(y_lo)/arb(x_hi)).lower(); hi = afin(arb(y_hi)/arb(x_lo)).upper()
    return lo, hi

def main_discard(p0, r):
    """for BD-chart boxes (MAIN and SHORT); returns reason or None"""
    if P_not_dart(p0, r): return 'not a dart (f_D >= 0 on box)'
    if P_prop51(p0, r): return 'Prop 5.1 (beta <= 0.0339)'
    if P_lemma51(p0, r): return 'Lemma 5.1 (tan bA + tan bC <= 0.1)'
    if P_mirror_tan(p0, r): return 'mirror (b_A > b_C)'
    return None

# ---- named predicates (each returns True only if it holds on the WHOLE box, rigorously) ----
def P_not_dart(p0, r): return bool(fD_bounds(p0, r)[0] >= 0)
def P_prop51(p0, r):
    b = beta_max(p0, r)
    return b is not None and bool(b <= dec('0.0339').lower())
def _tanAC(p0, r):
    exact_box(p0, r)
    if not (p0[1] - r[1] >= 0 and p0[3] + r[3] <= 0): return None, None      # A in the closed upper, C in the closed lower half plane
    tA = tanb_bounds(p0[0] - r[0], p0[0] + r[0], p0[1] - r[1], p0[1] + r[1])
    tC = tanb_bounds(p0[2] - r[2], p0[2] + r[2], -(p0[3] + r[3]), -(p0[3] - r[3]))
    return tA, tC
def P_lemma51(p0, r):
    """tan b_A + tan b_C <= 0.1 on the whole box, with both angles acute (a_x > 0, c_x > 0: LEMMAS 5.1 hypotheses)"""
    tA, tC = _tanAC(p0, r)
    if tA is None or tC is None: return False
    return bool(tA[1] + tC[1] <= dec('0.1').lower())
def P_mirror_tan(p0, r):
    """b_A > b_C on the whole box (both acute): tan b_A > tan b_C"""
    tA, tC = _tanAC(p0, r)
    if tA is None or tC is None: return False
    return bool(tA[0] > tC[1])
def P_tanbA_le(p0, r, t):
    exact_box(p0, r)
    if p0[0] - r[0] <= 0 or p0[1] - r[1] < 0: return False
    return bool(tanb_bounds(p0[0] - r[0], p0[0] + r[0], p0[1] - r[1], p0[1] + r[1])[1] <= dec(t).lower())
def P_ay_ge(p0, r, y):
    exact_box(p0, r)
    return bool(arb(p0[1] - r[1]) >= dec(y).upper())
import re
def check_reason(reason, p0, r):
    """re-evaluate a logged discard reason rigorously; returns True iff it holds on the whole box"""
    if reason == 'not a dart (f_D >= 0 on box)': return P_not_dart(p0, r)
    if reason in ('Prop 5.1', 'Prop 5.1 (beta <= 0.0339)'): return P_prop51(p0, r)      # the old runs wrote 'Prop 5.1'
    if reason == 'Lemma 5.1 (tan bA + tan bC <= 0.1)': return P_lemma51(p0, r)
    if reason == 'mirror (b_A > b_C)': return P_mirror_tan(p0, r)
    m = re.match(r'tan b_A <= ([0-9.]+) \(COLLAR/SPIKE certificates, C in CBOX\)$', reason)
    if m: return P_tanbA_le(p0, r, m.group(1)) and in_cbox(p0, r)
    m = re.match(r'a_y >= ([0-9.e+-]+) \(outside this run\)$', reason)
    if m: return P_ay_ge(p0, r, m.group(1))
    return False
CBOX = ((0.40, 0.60), (-0.966, -0.766))
def cbox_snapped():
    """the C-region actually covered by the COLLAR/SPIKE runs: CBOX snapped outward to the dyadic grid 2^-16 (cover2.dyadic_root)"""
    import math
    sc = 2.0**16
    return ((math.floor(CBOX[0][0]*sc)/sc, math.ceil(CBOX[0][1]*sc)/sc), (math.floor(CBOX[1][0]*sc)/sc, math.ceil(CBOX[1][1]*sc)/sc))
def in_cbox(p0, r):
    exact_box(p0, r)
    (x0, x1), (y0, y1) = cbox_snapped()
    return p0[2] - r[2] >= x0 and p0[2] + r[2] <= x1 and p0[3] - r[3] >= y0 and p0[3] + r[3] <= y1

# ---------------------------------------------------------------------------------------------------------------------------
# C-domain predicates for the general spike/collar runs (run_spike_all.py).  K1(alpha0, R) = {C : c_y < 0, angle CBD >= alpha0,
# angle CDB >= alpha0, |C| <= R}.  Each predicate returns True only if it holds on the WHOLE box (coordinates 2, 3 = c_x, c_y).
def _tan_deg(deg):
    return (dec(deg)*arb.pi()/180).tan()
def P_bC_lt(p0, r, deg):
    """angle CBD < deg on the whole box:  c_x > 0 and |c_y| < tan(deg) c_x  for all points"""
    exact_box(p0, r)
    cx_lo = arb(p0[2] - r[2]); cy_lo = arb(p0[3] - r[3])
    if not (cx_lo > 0 and p0[3] + r[3] <= 0): return False
    return bool(((-cy_lo) - _tan_deg(deg)*cx_lo).upper() < 0)
def P_aD_lt(p0, r, deg):
    """angle CDB < deg on the whole box:  1 - c_x > 0 and |c_y| < tan(deg) (1 - c_x)  for all points"""
    exact_box(p0, r)
    cx_hi = arb(p0[2] + r[2]); cy_lo = arb(p0[3] - r[3])
    if not ((1 - cx_hi) > 0 and p0[3] + r[3] <= 0): return False
    return bool(((-cy_lo) - _tan_deg(deg)*(1 - cx_hi)).upper() < 0)
def P_C_far(p0, r, R):
    """|C| > R on the whole box"""
    exact_box(p0, r)
    lo = [p0[2] - r[2], p0[3] - r[3]]; hi = [p0[2] + r[2], p0[3] + r[3]]
    near = [min(max(0.0, lo[i]), hi[i]) for i in range(2)]          # RIGOR: exact: closest point of the box to the origin (selection among exact floats)
    R_ = dec(R)
    return bool(afin(arb(near[0])*arb(near[0]) + arb(near[1])*arb(near[1]) - R_*R_).lower() > 0)
def P_cy_nonneg(p0, r):
    """c_y >= 0 somewhere is excluded by the charts; a box with c_y_hi >= 0 must not be certified -- returned as a discard only if
    the box lies in c_y >= 0 entirely (never for our roots)"""
    exact_box(p0, r)
    return p0[3] - r[3] >= 0
def k1_discard(alpha0, R):
    def disc(p0, r):
        if P_bC_lt(p0, r, alpha0): return f'C outside K1: angle CBD < {alpha0} deg'
        if P_aD_lt(p0, r, alpha0): return f'C outside K1: angle CDB < {alpha0} deg'
        if P_C_far(p0, r, R): return f'C outside K1: |C| > {R}'
        return None
    return disc
_check_reason_v1 = check_reason
def check_reason(reason, p0, r):
    m = re.match(r'C outside K1: angle CBD < ([0-9.]+) deg$', reason)
    if m: return P_bC_lt(p0, r, m.group(1))
    m = re.match(r'C outside K1: angle CDB < ([0-9.]+) deg$', reason)
    if m: return P_aD_lt(p0, r, m.group(1))
    m = re.match(r'C outside K1: \|C\| > ([0-9.]+)$', reason)
    if m: return P_C_far(p0, r, m.group(1))
    return _check_reason_v1(reason, p0, r)

# ---------------------------------------------------------------------------------------------------------------------------
# Polar C-sector predicates: S(bmin, bmax, rmin, rmax) = {C = rho (cos b, -sin b): bmin <= b <= bmax, rmin <= rho <= rmax}.
# b > th  <=>  g_th(C) := (-c_y) cos th - c_x sin th > 0   (b, th in (0, pi)); g_th is linear, so its range over a box is attained
# at corners.  Each predicate holds on the WHOLE box (rigorous, arb).
def _corners2(p0, r):
    exact_box(p0, r)
    return [(arb(p0[2] + s*r[2]), arb(p0[3] + q*r[3])) for s in (-1, 1) for q in (-1, 1)]
def _g(th_deg, cx, cy):
    th = dec(th_deg)*arb.pi()/180
    return (-cy)*th.cos() - cx*th.sin()
def P_bC_gt(p0, r, deg):
    """b_C > deg on the whole box (box in the closed lower half plane; the linear criterion is valid for c_y <= 0, C != 0)"""
    if not (p0[3] + r[3] <= 0): return False
    return all(_g(deg, cx, cy).lower() > 0 for cx, cy in _corners2(p0, r))
def P_bC_lt2(p0, r, deg):
    """b_C < deg on the whole box (box in the closed lower half plane)"""
    if not (p0[3] + r[3] <= 0): return False
    return all(_g(deg, cx, cy).upper() < 0 for cx, cy in _corners2(p0, r))
def P_C_near(p0, r, R):
    """|C| < R on the whole box (farthest corner)"""
    R_ = dec(R)
    return all(afin(cx*cx + cy*cy - R_*R_).upper() < 0 for cx, cy in _corners2(p0, r))
def sector_discard(bmin, bmax, rmin, rmax):
    def disc(p0, r):
        if P_bC_lt2(p0, r, bmin): return f'C outside sector: b_C < {bmin} deg'
        if P_bC_gt(p0, r, bmax): return f'C outside sector: b_C > {bmax} deg'
        if P_C_near(p0, r, rmin): return f'C outside sector: |C| < {rmin}'
        if P_C_far(p0, r, rmax): return f'C outside sector: |C| > {rmax}'
        return None
    return disc
_check_reason_v2 = check_reason
def check_reason(reason, p0, r):
    m = re.match(r'C outside sector: b_C < ([0-9.]+) deg$', reason)
    if m: return P_bC_lt2(p0, r, m.group(1))
    m = re.match(r'C outside sector: b_C > ([0-9.]+) deg$', reason)
    if m: return P_bC_gt(p0, r, m.group(1))
    m = re.match(r'C outside sector: \|C\| < ([0-9.]+)$', reason)
    if m: return P_C_near(p0, r, m.group(1))
    m = re.match(r'C outside sector: \|C\| > ([0-9.]+)$', reason)
    if m: return P_C_far(p0, r, m.group(1))
    return _check_reason_v2(reason, p0, r)

def P_C_in_sector(p0, r, bmin, bmax, rmin, rmax):
    """C-box entirely inside the closed sector {bmin <= b_C <= bmax, rmin <= |C| <= rmax} (rigorous)"""
    if not (p0[3] + r[3] <= 0): return False
    cs = _corners2(p0, r)
    if not all(_g(bmin, cx, cy).lower() >= 0 for cx, cy in cs): return False         # b >= bmin at all corners (linear)
    if not all(_g(bmax, cx, cy).upper() <= 0 for cx, cy in cs): return False         # b <= bmax
    rmax_ = dec(rmax); rmin_ = dec(rmin)
    if not all(afin(cx*cx + cy*cy - rmax_*rmax_).upper() <= 0 for cx, cy in cs): return False   # |C| <= rmax (convex)
    # |C| >= rmin on the whole box: the closest point of the box to the origin
    lo = [p0[2] - r[2], p0[3] - r[3]]; hi = [p0[2] + r[2], p0[3] + r[3]]
    near = [min(max(0.0, lo[i]), hi[i]) for i in range(2)]   # RIGOR: exact: closest point of the box to the origin (selection among exact floats)
    return bool(afin(arb(near[0])*arb(near[0]) + arb(near[1])*arb(near[1]) - rmin_*rmin_).lower() >= 0)
_check_reason_v3 = check_reason
def check_reason(reason, p0, r):
    m = re.match(r'C inside sector \[([0-9.]+), ([0-9.]+)\] deg x \[([0-9.]+), ([0-9.]+)\] \(covered by run ([A-Za-z0-9_]+)\)$', reason)
    if m: return P_C_in_sector(p0, r, m.group(1), m.group(2), m.group(3), m.group(4))
    return _check_reason_v3(reason, p0, r)

# ---------------------------------------------------------------------------------------------------------------------------
# Angle predicates for 4-D BD-chart boxes (A = coords 0,1 with a_y > 0; C = coords 2,3 with c_y < 0).  Angles b_A = arg A,
# b_C = -arg C in (0, pi).  For a box not containing the origin, the extreme angles are attained at corners.
def _ang_range(xlo, xhi, ylo, yhi, lower):
    """[min, max] (arb) of the angle atan2(|y|, x) over the box (box in the closed upper (lower=False) or lower half plane,
    not containing the origin)"""
    vals = []
    for x in (xlo, xhi):
        for y in (ylo, yhi):
            vals.append(afin(arb.atan2(arb(abs(y)), arb(x))))
    return amin(*vals), amax(*vals)
def angA(p0, r):
    exact_box(p0, r)
    return _ang_range(p0[0] - r[0], p0[0] + r[0], p0[1] - r[1], p0[1] + r[1], False)
def angC(p0, r):
    exact_box(p0, r)
    return _ang_range(p0[2] - r[2], p0[2] + r[2], -(p0[3] + r[3]), -(p0[3] - r[3]), True)
def _box_has_origin(xlo, xhi, ylo, yhi): return xlo <= 0 <= xhi and ylo <= 0 <= yhi
def P_mirror(p0, r):
    """b_C < b_A on the whole box (the mirror image (C bar, A bar) has b_A' = b_C < b_C' = b_A: covered by the b_A <= b_C domain)"""
    exact_box(p0, r)
    if p0[1] - r[1] <= 0 or p0[3] + r[3] >= 0: return False
    aA = angA(p0, r); aC = angC(p0, r)
    return bool(aC[1] < aA[0])
def P_tanbA_lt(p0, r, t):
    """tan b_A < t (b_A < atan t) on the whole box: covered by the spike families (any |A|)"""
    exact_box(p0, r)
    if p0[1] - r[1] < 0: return False
    aA = angA(p0, r)
    return bool(aA[1] < dec(t).atan())
def P_bC_ge(p0, r, deg):
    """b_C >= deg on the whole box"""
    exact_box(p0, r)
    if p0[3] + r[3] > 0: return False
    aC = angC(p0, r)
    return bool(aC[0] >= dec(deg)*arb.pi()/180)
def P_thetaB_le(p0, r, rad):
    """theta_B = b_A + b_C <= rad on the whole box (Lemma 5.3 for rad = 0.1772)"""
    exact_box(p0, r)
    if p0[1] - r[1] <= 0 or p0[3] + r[3] >= 0: return False
    return bool((angA(p0, r)[1] + angC(p0, r)[1]) <= dec(rad))
_check_reason_v4 = check_reason
def check_reason(reason, p0, r):
    if reason == 'mirror (b_C < b_A on the box)': return P_mirror(p0, r)
    m = re.match(r'tan b_A < ([0-9.]+) on the box \(spike families\)$', reason)
    if m: return P_tanbA_lt(p0, r, m.group(1))
    m = re.match(r'Lemma 5.2 \(b_C >= ([0-9.]+) deg, theta_B >= b_C\)$', reason)
    if m: return decq(m.group(1)) >= decq('174.29') and P_bC_ge(p0, r, m.group(1))
    m = re.match(r'Lemma 5.3 \(theta_B <= ([0-9.]+)\)$', reason)
    if m: return decq(m.group(1)) <= decq('0.1772') and P_thetaB_le(p0, r, m.group(1))
    return _check_reason_v4(reason, p0, r)
