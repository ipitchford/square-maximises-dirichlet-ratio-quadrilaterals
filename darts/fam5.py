"""Parameter families for cover5 (README, section 2: the global decomposition).  A family maps a box p0 +- r in its 4 parameters
to (i) the arm specs of bo2 (for the 1-D lower bound), (ii) candidate INNER darts contained in every dart of the box (for the P2
upper bound of lambda_2), (iii) optionally a Cartesian BD-chart box containing every dart of the box (for the FEM certificate),
(iv) rigorous discard predicates (each returns a reason string re-checked by family.check(reason, p0, r) in replay).

Coordinates (angles in radians, l = ln rho):
  PP : (bA, lA, bC, lC)    A = e^lA (cos bA, sin bA), C = e^lC (cos bC, -sin bC)
  SP : (t,  0,  bC, lC)    A = a (1, t) with a >= AMIN (spike: AMIN = 1) or any a > 0 (collar: AMIN = 0), C polar
  PF : (bA, lA, bC, 0)     A polar, C = rho (cos bC, -sin bC) with rho >= RF (far)
  FP : (bA, 0,  bC, lC)    A far (rho >= RF), C polar
  FF : (bA, 0,  bC, 0)     both far
  SF : (t,  0,  bC, 0)     spike A, far C
All boxes are dyadic in these coordinates; a box with r_i = 0 in a coordinate does not use it."""
import math, numpy as np
from flint import arb
from nanguard import fin, afin, NaNError, InexactBox, fdown, fup, amax, amin, cuts, enclose, exact_box, dec
def _holds(f, p0, r):
    """a predicate that meets a non-finite value or an inexact box is NOT established (fail closed)"""
    exact_box(p0, r)
    try: return bool(f(p0, r))
    except NaNError: return False
from bo2 import Arm

LN = math.log
def exp_lo(l): return fdown(arb(l).exp())        # directed: a float <= e^l (l a float or an exact ball)
def exp_hi(l): return fup(arb(l).exp())          # directed: a float >= e^l
def pos(x): return x if x > 0.0 else 0.0

FEM_RMAX = 8.0
LONG_L_SCALE = 4.0       # split scale of ln(rho) for arms with rho >= 4 on the box (cover heuristic)
_BD_CACHE = {}
def _bd_problem_cached(k):
    from charts import bd_problem
    if k not in _BD_CACHE: _BD_CACHE[k] = bd_problem(k)
    return _BD_CACHE[k]
def ang_scale(b): return max(min(abs(b), abs(math.pi - b)), 0.02)   # RIGOR: heuristic: split scaling of the cover (no effect on rigour)
class Family:
    name = '?'
    fem = False
    def __init__(self, **kw):
        self.kw = kw
        for k, v in kw.items(): setattr(self, k, v)
    # ---- to be provided
    def arms(self, p0, r): raise NotImplementedError
    def inners(self, p0, r): raise NotImplementedError
    def cart_box(self, p0, r): return None
    def discard(self, p0, r): return None
    def check(self, reason, p0, r): return False
    def dirs(self, p0, r): raise NotImplementedError
    def scale(self, p, mode='bo', r=None): return np.ones(4)
    def fem_plausible(self, p0, r): return True          # cover heuristic only (no effect on rigour)

# ------------------------------------------------------------------------------------------------ shared rigorous predicates
PI = arb.pi()
def ival(lo, hi): return arb(lo).union(arb(hi))
def beta_cond_lower(bA, rhoA, bC, rhoC):
    """lower bound (arb) of g = sin bC / rhoA + sin bA / rhoC - 0.9661 sin(bA + bC) over the ranges (arb balls).
    g >= 0 everywhere  <=>  no dart with beta > 0.0339 (the dart condition is g_0 < 0 with 0.9661 replaced by 1, and
    1 - beta = (sin bC/rhoA + sin bA/rhoC)/sin(theta_B)); used for 'no dart with beta > 0.0339 (Prop 5.1 or not a dart)'."""
    return (bC.sin()/rhoA + bA.sin()/rhoC - arb('0.9661')*(bA + bC).sin())

def polar_ranges(p0, r, iA, iC):
    bA = ival(p0[iA[0]] - r[iA[0]], p0[iA[0]] + r[iA[0]]); bC = ival(p0[iC[0]] - r[iC[0]], p0[iC[0]] + r[iC[0]])
    return bA, bC

def subdiv_min(fun, lohi, n=4):
    """min over a product of intervals of a function of arb balls, by uniform subdivision (n per axis): rigorous lower bound"""
    import itertools
    axes = []
    for lo, hi in lohi:
        if hi > lo:
            xs = cuts(lo, hi, n); axes.append([(xs[i], xs[i + 1]) for i in range(n)])
        else: axes.append([(lo, hi)])
    return amin(*[afin(fun(*[ival(a, b) for a, b in cell])) for cell in itertools.product(*axes)])

# ------------------------------------------------------------------------------------------------ PP: polar x polar (compact)
class PP(Family):
    name = 'PP'; fem = True; inc = True; fem_cheap = True
    def arms(self, p0, r):
        exact_box(p0, r)
        bA1, bA2 = p0[0] - r[0], p0[0] + r[0]; bC1, bC2 = p0[2] - r[2], p0[2] + r[2]
        A_ = Arm(('polar', bA1, bA2, exp_lo(p0[1] - r[1]), exp_hi(p0[1] + r[1])), +1)
        C_ = Arm(('polar', bC1, bC2, exp_lo(p0[3] - r[3]), exp_hi(p0[3] + r[3])), -1)
        return A_, C_
    def inners(self, p0, r):
        A_, C_ = self.arms(p0, r)
        Ai = A_.inner_vertex(); Ci = C_.inner_vertex()
        out = []
        if Ai is not None and Ci is not None: out.append((Ai, Ci))
        if Ci is not None: out.append((None, Ci))
        return out
    def cart_box(self, p0, r):
        """Cartesian box containing {rho e(b)} over the polar box (arb enclosure, outward floats); None (no FEM) for arms longer
        than FEM_RMAX (the Crouzeix-Raviart / Liu bound on the affine reference meshes is useless for very long triangles)"""
        exact_box(p0, r)
        if not (p0[1] + r[1] <= math.log(FEM_RMAX) and p0[3] + r[3] <= math.log(FEM_RMAX)): return None   # routing only
        balls = []
        for (ib, il, side) in ((0, 1, 1), (2, 3, -1)):
            b = ival(p0[ib] - r[ib], p0[ib] + r[ib])
            rho = ival(exp_lo(p0[il] - r[il]), exp_hi(p0[il] + r[il]))
            balls += [afin(rho*b.cos()), afin(side*rho*b.sin())]
        if not (balls[1].lower() > 0 and balls[3].upper() < 0): return None
        enc = [enclose(x.lower(), x.upper()) for x in balls]        # exact-endpoint float boxes containing the enclosures
        return [c for c, h in enc], [h for c, h in enc]
    def dirs(self, p0, r):
        bA, bC = p0[0], p0[2]
        rA, rC = math.exp(p0[1]), math.exp(p0[3])
        ax, ay = rA*math.cos(bA), rA*math.sin(bA); cx, cy = rC*math.cos(bC), -rC*math.sin(bC)
        sp = [0.0, math.pi/2, math.atan2(cy - ay, cx - ax), math.atan2(cy - ay, cx - ax) + math.pi/2, bA, -bC,
              math.atan2(ay, ax - 1), math.atan2(cy, cx - 1), 0.5*(bA - bC)]
        return sp + [i*math.pi/8 - math.pi/2 for i in range(8)]
    def scale(self, p, mode='bo', r=None):
        # a long arm (rho >= 4 on the box) is treated by every method as the far arm rho >= rho_min (BO2 far variant, FARL chart),
        # so its ln(rho) coordinate is split only late (scale LONG_L_SCALE: fat long arms still need it, e.g. PP_2_10_2_11)
        rr_ = r if r is not None else np.zeros(4)
        longA = math.exp(p[1] - rr_[1]) >= 4.0; longC = math.exp(p[3] - rr_[3]) >= 4.0
        if mode == 'fem':
            return np.array([1.0, LONG_L_SCALE if longA else 1.0, 1.0, LONG_L_SCALE if longC else 1.0])   # RIGOR: heuristic: split scaling of the cover
        # BO2: the transverse spread of a polar arm over an angular box is rho * width, relative to the arm thickness ~ sin b;
        # capped at rho = 4, beyond which the far-arm variant (spread-free chord and tail bounds) is used (PP_1_10_1_7: 364 -> 281 s)
        rA = math.exp(p[1] + rr_[1]); rC = math.exp(p[3] + rr_[3])
        return np.array([ang_scale(p[0])/min(max(1.0, rA), 4.0), LONG_L_SCALE if longA else 1.0,   # RIGOR: heuristic: split scaling of the cover
                         ang_scale(p[2])/min(max(1.0, rC), 4.0), LONG_L_SCALE if longC else 1.0])   # RIGOR: heuristic: split scaling of the cover
    def fem_plausible(self, p0, r):
        """cover heuristic only (no effect on rigour): the BD-chart FEM certificate fails for long arms (Liu's bound on long
        elements, strain of thin pieces): observed failing throughout rho_C in [4.5, 7.4] (PP_0_4_0_11), useful up to 4.5; for
        one long arm (rho >= 4 on the box) the FARL chart is used instead (fem_setup)"""
        return (p0[1] + r[1] <= math.log(5.0) and p0[3] + r[3] <= math.log(5.0)) or self._long(p0, r) in ('A', 'C')
    def _long(self, p0, r):
        # long: rho >= 4 on the whole box and rho > 5 somewhere (boxes with rho <= 5 keep the BD chart)
        lA = exp_lo(p0[1] - r[1]) >= 4.0 and p0[1] + r[1] > math.log(5.0) and p0[0] + r[0] <= math.pi/2
        lC = exp_lo(p0[3] - r[3]) >= 4.0 and p0[3] + r[3] > math.log(5.0) and p0[2] + r[2] <= math.pi/2
        return 'A' if lA and not lC else ('C' if lC and not lA else None)
    def fem_setup(self, p0, r, k):
        """FEM chart for the box: the BD chart on the Cartesian enclosure (cart_box) when both arms are at most 5 long; the FARL
        chart of the far families (charts.farl_problem, LEMMAS 9.3/§10) when exactly one arm has rho >= 4 on the box: that arm is
        treated as the far arm rho >= RF := its smallest radius on the box (a family containing the box's darts), the other apex
        as a Cartesian box (mirrored for a long C arm).  Returns (problem, pc, rc) or None."""
        exact_box(p0, r)
        lng = self._long(p0, r)
        if lng == 'A':
            Cbox = _polar_cart_lower(p0[2] - r[2], p0[2] + r[2], exp_lo(p0[3] - r[3]), exp_hi(p0[3] + r[3]))
            if Cbox[1][1] >= 0: return None
            return farl_setup(p0[0] - r[0], p0[0] + r[0], Cbox, k, RFv=exp_lo(p0[1] - r[1]))
        if lng == 'C':
            Abox = _polar_cart_lower(p0[0] - r[0], p0[0] + r[0], exp_lo(p0[1] - r[1]), exp_hi(p0[1] + r[1]))
            if Abox[1][1] >= 0: return None
            return farl_setup(p0[2] - r[2], p0[2] + r[2], Abox, k, RFv=exp_lo(p0[3] - r[3]))
        cb = self.cart_box(p0, r)
        if cb is None: return None
        return _bd_problem_cached(k), cb[0], cb[1]
    # --- discards: ordered table reason -> predicate (each predicate holds on the WHOLE box, rigorously)
    def preds(self):
        def bb(p0, r): return ival(p0[0] - r[0], p0[0] + r[0]), ival(p0[2] - r[2], p0[2] + r[2])
        def not_dart_pi(p0, r):
            bA, bC = bb(p0, r); return (bA + bC).lower() >= PI.upper()
        def mirror(p0, r):
            # exact (box corners are dyadic): boxes touching the diagonal b_A = b_C must be decidable (the union-based ival
            # rounds outward and would leave them undecided forever)
            return bool(((arb(p0[2]) + arb(r[2])) - (arb(p0[0]) - arb(r[0]))).upper() <= 0)
        def l52(p0, r):
            bA, bC = bb(p0, r); return bool((bA + bC).lower() >= (arb('174.29')*PI/180).upper())
        def l53(p0, r):
            bA, bC = bb(p0, r); return bool((bA + bC).upper() <= arb('0.1772').lower())
        def collar(p0, r):
            bA, bC = bb(p0, r); return bool(bA.upper() < arb('0.02').atan().lower())
        def spike(p0, r):
            bA, bC = bb(p0, r); rA_lo = exp_lo(p0[1] - r[1])
            return bool(bA.upper() < arb('0.3').atan().lower() and (arb(rA_lo)*bA.cos()).lower() >= 1)
        def nodart(p0, r):
            # g = sin bC/rhoA + sin bA/rhoC - 0.9661 sin(bA + bC) is non-increasing in rhoA and in rhoC (sin b >= 0 for b in
            # [0, pi]): its minimum over the box is attained at the largest radii (outward-rounded floats rho_max, which can only
            # lower g); the two angles are subdivided into 6 x 6 cells (exact float cell endpoints, cells cover the box)
            if not (p0[0] - r[0] >= 0 and p0[2] - r[2] >= 0 and p0[0] + r[0] <= math.pi and p0[2] + r[2] <= math.pi): return False
            rA = arb(exp_hi(p0[1] + r[1])); rC = arb(exp_hi(p0[3] + r[3]))
            n = 6
            xa = cuts(p0[0] - r[0], p0[0] + r[0], n); xc = cuts(p0[2] - r[2], p0[2] + r[2], n)
            for i in range(n):
                bA = ival(xa[i], xa[i + 1])
                for j in range(n):
                    if not (beta_cond_lower(bA, rA, ival(xc[j], xc[j + 1]), rC).lower() >= 0): return False
            return True
        return [('theta_B >= pi on the box (not a dart)', not_dart_pi), ('mirror (b_C <= b_A on the box)', mirror),
                ('mirror (b_C < b_A on the box)', mirror),
                ('Lemma 5.2 (theta_B >= 174.29 deg)', l52), ('Lemma 5.3 (theta_B <= 0.1772)', l53),
                ('collar family (tan b_A < 0.02 on the box)', collar),
                ('spike family (tan b_A < 0.3 and a_x >= 1 on the box)', spike),
                ('no dart with beta > 0.0339 (Prop 5.1 or not a dart)', nodart)]
    def discard(self, p0, r):
        for why, f in self.preds():
            if _holds(f, p0, r): return why
        return None
    def check(self, reason, p0, r):
        for why, f in self.preds():
            if why == reason: return _holds(f, p0, r)
        return False

FAMILIES = {'PP': PP}

# ------------------------------------------------------------------------------------------------ common predicate builders
KF_SECTOR = ('38', '100', '0.5', '1.8')     # degrees, radii (exact decimals): the FEM runs certs2/kf_* (COLLAR / SPIKE)
# t-ranges covered by the KF runs (verified from their manifests in accounting.py): COLLAR t <= KF_COL_T_HI (the collar chart is
# evaluated at t_hi = float(0.02)); SPIKE: the dyadic snapping (cover2.dyadic_root) of the t slabs [0.02, 0.3]
KF_COL_T_HI = 0.02
COLLAR_T1 = 0.02             # the t of the collar FEM chart (charts.spike_problem via run_fam.fem_problem('SP_col'))
KF_SPK_T = (1310/65536, 19661/65536)
def _deg(x): return arb(x)*PI/180
def c_polar_in_sector(bC, rho_lo, rho_hi, sec):
    """C polar box (angle ball bC, radius range [rho_lo, rho_hi] floats) inside the closed sector sec = (bmin, bmax, rmin, rmax),
    the sector constants read as exact decimals"""
    return bool(bC.lower() >= _deg(dec(sec[0])).upper() and bC.upper() <= _deg(dec(sec[1])).lower()
                and arb(rho_lo) >= dec(sec[2]).upper() and arb(rho_hi) <= dec(sec[3]).lower())

# ------------------------------------------------------------------------------------------------ SP: spike A x polar C
class SP(Family):
    """A = a (1, t); variant 'col': t in (0, 0.02], any a > 0;  variant 'spk': t in [0.02, 0.3], a >= 1.  C polar.
    Coordinates (t, 0, bC, lC)."""
    name = 'SP'; fem = True
    def __init__(self, variant='spk', **kw):
        super().__init__(variant=variant, **kw)
        self.name = 'SP_' + variant
    def arms(self, p0, r):
        exact_box(p0, r)
        t1, t2 = pos(p0[0] - r[0]), p0[0] + r[0]
        A_ = Arm(('far', fdown(arb(t1).atan()) if t1 > 0 else 0.0, fup(arb(t2).atan()), 0.0), +1)   # directed atan
        C_ = Arm(('polar', p0[2] - r[2], p0[2] + r[2], exp_lo(p0[3] - r[3]), exp_hi(p0[3] + r[3])), -1)
        return A_, C_
    def inners(self, p0, r):
        A_, C_ = self.arms(p0, r)
        Ci = C_.inner_vertex()
        out = []
        t1 = p0[0] - r[0]                                    # exact (exact box)
        if self.variant == 'spk' and t1 > 0 and Ci is not None:
            # tri(B, D, (1, t1)) is inside tri(B, D, a (1, t)) for a >= 1, t >= t1 (LEMMAS 9.4; closed containment)
            out.append(((1.0, t1), Ci))
        if Ci is not None: out.append((None, Ci))
        return out
    def cart_box(self, p0, r):
        exact_box(p0, r)
        b = ival(p0[2] - r[2], p0[2] + r[2]); rho = ival(exp_lo(p0[3] - r[3]), exp_hi(p0[3] + r[3]))
        x = afin(rho*b.cos()); y = afin(-rho*b.sin())
        if not y.upper() < 0: return None
        cx, cy = enclose(x.lower(), x.upper()), enclose(y.lower(), y.upper())
        if self.variant == 'col':
            # the collar chart is evaluated at t1 = float(0.02) and covers every t <= t1 (Lemma 9.2): only boxes with t <= t1
            # (the snapped collar roots reach 1311/65536 > t1; their leaves above t1 are certified by BO2 or split)
            if not (p0[0] + r[0] <= COLLAR_T1): return None
            t = (COLLAR_T1, 0.0)
        else: t = (p0[0], r[0])                               # exact box: t-range [p0 - r, p0 + r] exactly
        return [t[0], 0.0, cx[0], cy[0]], [t[1], 0.0, cx[1], cy[1]]
    def dirs(self, p0, r):
        bC = p0[2]; rC = math.exp(p0[3]); cx, cy = rC*math.cos(bC), -rC*math.sin(bC)
        sp = [0.0, math.atan2(cy, cx - 1), -bC, math.atan2(cy, cx), 0.5*(math.atan(p0[0]) - bC), math.pi/2, -math.pi/12, math.pi/12]
        return sp + [i*math.pi/12 - math.pi/2 for i in range(12)]
    def scale(self, p, mode='bo', r=None):
        if mode == 'fem': return np.array([max(p[0], 0.02), 1.0, 1.0, 1.0])   # RIGOR: heuristic: split scaling of the cover
        return np.array([max(p[0], 0.02), 1.0, ang_scale(p[2]), 1.0])   # RIGOR: heuristic: split scaling of the cover
    def preds(self):
        def tb(p0, r):
            t1, t2 = pos(p0[0] - r[0]), p0[0] + r[0]
            bA = ival(fdown(arb(t1).atan()) if t1 > 0 else 0.0, fup(arb(t2).atan()))
            return bA, ival(p0[2] - r[2], p0[2] + r[2])
        def mirror(p0, r):
            bA, bC = tb(p0, r); return bool(bC.upper() <= bA.lower())
        def l53(p0, r):
            bA, bC = tb(p0, r); return bool((bA + bC).upper() <= arb('0.1772').lower())
        def l52(p0, r):
            bA, bC = tb(p0, r); return bool((bA + bC).lower() >= (arb('174.29')*PI/180).upper())
        def kf(p0, r):
            if self.variant == 'col' and not (p0[0] + r[0] <= KF_COL_T_HI): return False
            if self.variant == 'spk' and not (p0[0] - r[0] >= KF_SPK_T[0] and p0[0] + r[0] <= KF_SPK_T[1]): return False
            bC = ival(p0[2] - r[2], p0[2] + r[2])
            return c_polar_in_sector(bC, exp_lo(p0[3] - r[3]), exp_hi(p0[3] + r[3]), KF_SECTOR)
        def nodart(p0, r):
            # Lemma 13.1(a) for spikes / collars A = a (1, t) (rho_A unbounded): g = sin bC/rhoA + sin bA/rhoC - 0.9661 sin(bA + bC)
            # >= sin bA/rhoC - 0.9661 sin(bA + bC) (the dropped term is >= 0 for bC in [0, pi]), and that bound decreases in rhoC:
            # evaluate it at the largest radius over 6 x 6 cells of (t, bC) (bA = atan t on each t-cell, outward)
            if not (p0[2] - r[2] >= 0 and p0[2] + r[2] <= math.pi): return False
            rC = arb(exp_hi(p0[3] + r[3]))
            n = 6
            ts = cuts(pos(p0[0] - r[0]), p0[0] + r[0], n); cs = cuts(p0[2] - r[2], p0[2] + r[2], n)
            for i in range(n):
                bA = ival(fdown(arb(ts[i]).atan()) if ts[i] > 0 else 0.0, fup(arb(ts[i + 1]).atan()))
                for j in range(n):
                    bC = ival(cs[j], cs[j + 1])
                    if not ((bA.sin()/rC - arb('0.9661')*(bA + bC).sin()).lower() >= 0): return False
            return True
        return [('mirror (b_C <= b_A on the box)', mirror), ('mirror (b_C < b_A on the box)', mirror),
                ('Lemma 5.3 (theta_B <= 0.1772)', l53),
                ('Lemma 5.2 (theta_B >= 174.29 deg)', l52), ('C inside the KF sector (certs2/kf_collar, certs2/kf_spike)', kf),
                ('no dart with beta > 0.0339 (Prop 5.1 or not a dart)', nodart)]
    def discard(self, p0, r):
        for why, f in self.preds():
            if _holds(f, p0, r): return why
        return None
    def check(self, reason, p0, r):
        for why, f in self.preds():
            if why == reason: return _holds(f, p0, r)
        return False

# ------------------------------------------------------------------------------------------------ far families
RF = 20.0
class FarMixin:
    def _armA(self, p0, r):
        exact_box(p0, r)
        if self.Afar: return Arm(('far', p0[0] - r[0], p0[0] + r[0], RF), +1)
        return Arm(('polar', p0[0] - r[0], p0[0] + r[0], exp_lo(p0[1] - r[1]), exp_hi(p0[1] + r[1])), +1)
    def _armC(self, p0, r):
        exact_box(p0, r)
        if self.Cfar: return Arm(('far', p0[2] - r[2], p0[2] + r[2], RF), -1)
        return Arm(('polar', p0[2] - r[2], p0[2] + r[2], exp_lo(p0[3] - r[3]), exp_hi(p0[3] + r[3])), -1)

class PF(FarMixin, Family):
    """A polar or far, C polar or far, at least one far.  Coordinates (bA, lA|0, bC, lC|0).  FEM fallback (one far arm only):
    the FARL chart (charts.farl_problem) via fem_setup (mirrored for a far C arm)."""
    @property
    def fem(self): return self.Afar != self.Cfar
    def __init__(self, Afar, Cfar, **kw):
        super().__init__(Afar=Afar, Cfar=Cfar, **kw)
        self.name = ('F' if Afar else 'P') + ('F' if Cfar else 'P')
    def arms(self, p0, r): return self._armA(p0, r), self._armC(p0, r)
    def inners(self, p0, r):
        A_, C_ = self.arms(p0, r)
        Ai = A_.inner_vertex(); Ci = C_.inner_vertex()
        out = []
        if Ai is not None and Ci is not None: out.append((Ai, Ci))
        if Ci is not None: out.append((None, Ci))
        return out
    def dirs(self, p0, r):
        bA, bC = p0[0], p0[2]
        sp = [bA, -bC, 0.5*(bA - bC), 0.5*(bA - bC) + math.pi/2, 0.0, math.pi/2, bA + math.pi/2, -bC + math.pi/2]
        return sp + [i*math.pi/24 - math.pi/2 for i in range(24)]
    def scale(self, p, mode='bo', r=None):
        # a long polar arm (rho >= 4 on the box) is treated as far (BO2 far variant, FARL): split its ln(rho) coordinate late
        rr_ = r if r is not None else np.zeros(4)
        lA = LONG_L_SCALE if (not self.Afar and math.exp(p[1] - rr_[1]) >= 4.0) else 1.0
        lC = LONG_L_SCALE if (not self.Cfar and math.exp(p[3] - rr_[3]) >= 4.0) else 1.0
        if mode == 'fem': return np.array([1.0, lA, 1.0, lC])   # RIGOR: heuristic: split scaling of the cover
        return np.array([ang_scale(p[0]), lA, ang_scale(p[2]), lC])   # RIGOR: heuristic: split scaling of the cover
    def preds(self):
        def bb(p0, r): return ival(p0[0] - r[0], p0[0] + r[0]), ival(p0[2] - r[2], p0[2] + r[2])
        def not_dart_pi(p0, r):
            bA, bC = bb(p0, r); return bool((bA + bC).lower() >= PI.upper())
        def mirror(p0, r):
            # exact (box corners are dyadic): boxes touching the diagonal b_A = b_C must be decidable (the union-based ival
            # rounds outward and would leave them undecided forever)
            return bool(((arb(p0[2]) + arb(r[2])) - (arb(p0[0]) - arb(r[0]))).upper() <= 0)
        def l52(p0, r):
            bA, bC = bb(p0, r); return bool((bA + bC).lower() >= (arb('174.29')*PI/180).upper())
        def l53(p0, r):
            bA, bC = bb(p0, r); return bool((bA + bC).upper() <= arb('0.1772').lower())
        def spike(p0, r):
            bA, bC = bb(p0, r)
            if not bA.upper() < arb('0.3').atan().lower(): return False
            if self.Afar: return bool((arb(RF)*bA.cos()).lower() >= 1)                     # a_x >= RF cos b_A >= 1
            return bool((arb(exp_lo(p0[1] - r[1]))*bA.cos()).lower() >= 1)
        def collar(p0, r):
            bA, bC = bb(p0, r); return bool(bA.upper() < arb('0.02').atan().lower())
        def nodart(p0, r):
            # g = sin bC/rhoA + sin bA/rhoC - 0.9661 sin(bA + bC); far arms: 1/rho in [0, 1/RF]
            def inv(il, far):
                if far: return ival(0.0, 1/arb(RF))                                         # 1/rho in [0, 1/RF] (balls)
                return ival(1/arb(exp_hi(p0[il] + r[il])), 1/arb(exp_lo(p0[il] - r[il])))
            bA, bC = bb(p0, r)
            g = bC.sin()*inv(1, self.Afar) + bA.sin()*inv(3, self.Cfar) - arb('0.9661')*(bA + bC).sin()
            return bool(g.lower() >= 0)
        return [('theta_B >= pi on the box (not a dart)', not_dart_pi), ('mirror (b_C <= b_A on the box)', mirror),
                ('mirror (b_C < b_A on the box)', mirror),
                ('Lemma 5.2 (theta_B >= 174.29 deg)', l52), ('Lemma 5.3 (theta_B <= 0.1772)', l53),
                ('collar family (tan b_A < 0.02 on the box)', collar),
                ('spike family (tan b_A < 0.3 and a_x >= 1 on the box)', spike),
                ('no dart with beta > 0.0339 (Prop 5.1 or not a dart)', nodart)]
    def discard(self, p0, r):
        for why, f in self.preds():
            if _holds(f, p0, r): return why
        return None
    def check(self, reason, p0, r):
        for why, f in self.preds():
            if why == reason: return _holds(f, p0, r)
        return False

class SF(SP):
    """spike A (variant col/spk) x far C.  Coordinates (t, 0, bC, 0)."""
    fem = False
    def __init__(self, variant='spk', **kw):
        super().__init__(variant=variant, **kw); self.name = 'SF_' + variant
    def arms(self, p0, r):
        exact_box(p0, r)
        t1, t2 = pos(p0[0] - r[0]), p0[0] + r[0]
        A_ = Arm(('far', fdown(arb(t1).atan()) if t1 > 0 else 0.0, fup(arb(t2).atan()), 0.0), +1)   # directed atan
        C_ = Arm(('far', p0[2] - r[2], p0[2] + r[2], RF), -1)
        return A_, C_
    def cart_box(self, p0, r): return None
    def dirs(self, p0, r):
        bC = p0[2]
        sp = [0.0, -bC, 0.5*(math.atan(p0[0]) - bC), -bC + math.pi/2, 0.5*(math.atan(p0[0]) - bC) + math.pi/2, math.pi/2]
        return sp + [i*math.pi/24 - math.pi/2 for i in range(24)]
    def preds(self):
        base = [x for x in SP.preds(self) if not x[0].startswith('C inside the KF')]
        return base

FAMILIES = {'PP': PP, 'SP_col': lambda: SP('col'), 'SP_spk': lambda: SP('spk'),
            'PF': lambda: PF(False, True), 'FP': lambda: PF(True, False), 'FF': lambda: PF(True, True),
            'SF_col': lambda: SF('col'), 'SF_spk': lambda: SF('spk')}

# ------------------------------------------------------------------------------------------------ Lemma 12.1 predicate
L121_CMAX = '0.05'              # exact decimal threshold of Lemma 12.1
def stem_chord_bound(armA, armC, K=48):
    """rigorous upper bound (float) of c* = max_{0<=x<=1} (h_A + h_C) over all darts of the arm specs (vertical slices,
    Lemma 11.3(iii) slab extents); None if an arm cannot be bounded (far arm not transversal to the vertical slices)"""
    from bo1d import slab_extent
    u = (arb(1), arb(0))
    pts = []
    for a in (armA, armC):
        if a is None: pts.append(None); continue
        if a.kind == 'polar': pts.append(a.outer_points())
        else:
            s = a.dir_sign(u)
            if s == 0: return None
            pts.append(a.outer_points(u, 1.0 if s > 0 else 0.0))
    nodes = cuts(0.0, 1.0, K)                                # cells cover [0, 1] exactly
    tots = []
    for k in range(K):
        tot = arb(0)
        for P in pts:
            if P is not None: tot = tot + slab_extent(P, u, nodes[k], nodes[k + 1])
        tots.append(afin(tot))
    return amax(*tots)                                       # arb upper bound (exact point)

def lemma121_holds(armA, armC):
    c = stem_chord_bound(armA, armC)
    return c is not None and bool(c <= dec(L121_CMAX).lower())

def _add_l121(cls):
    old = cls.preds
    def preds(self):
        lst = list(old(self))
        def l121(p0, r):
            try: armA, armC = self.arms(np.array(p0, float), np.array(r, float))   # RIGOR: exact: float64 arrays of Python floats (no rounding)
            except Exception: return False
            return lemma121_holds(armA, armC)
        lst.append(('Lemma 12.1 (vertical stem chords <= 0.05 on the box)', l121))
        return lst
    cls.preds = preds
    return cls
for _c in (PP, SP, PF): _add_l121(_c)        # SF inherits SP's list

# ------------------------------------------------------------------------------------------------ earlier certified regions
SHALLOW_BOX = ((0.965, 1.262), (0.085, 0.796), (0.479, 1.213), (-1.091, -0.483))   # certs/shallow_crit.jsonl (Theorem H)
def _cart_inside(cb, box):
    """the (exact-endpoint) Cartesian box cb inside the float box 'box' (compared exactly)"""
    if cb is None: return False
    pc, rc = cb
    exact_box(pc, rc)
    return all(pc[i] - rc[i] >= box[i][0] and pc[i] + rc[i] <= box[i][1] for i in range(4))
def _add_old_regions(cls):
    old = cls.preds
    def preds(self):
        lst = list(old(self))
        def shallow(p0, r):
            return _cart_inside(self.cart_box(np.array(p0, float), np.array(r, float)), SHALLOW_BOX)   # RIGOR: exact: float64 arrays of Python floats (no rounding)
        def maincbox(p0, r):
            from regions import cbox_snapped
            cb = self.cart_box(np.array(p0, float), np.array(r, float))   # RIGOR: exact: float64 arrays of Python floats (no rounding)
            if cb is None: return False
            pc, rc = cb
            exact_box(pc, rc)
            (x0, x1), (y0, y1) = cbox_snapped()
            if not (pc[2] - rc[2] >= x0 and pc[2] + rc[2] <= x1 and pc[3] - rc[3] >= y0 and pc[3] + rc[3] <= y1): return False
            if not pc[1] + rc[1] <= 1.5479: return False
            bA = ival(p0[0] - r[0], p0[0] + r[0])
            return bool(bA.lower() >= arb('0.3').atan().upper())
        # insert before the (expensive) no-dart predicate
        i = [w for w, _ in lst].index('no dart with beta > 0.0339 (Prop 5.1 or not a dart)')
        lst[i:i] = [('inside the shallow box of Theorem H (certs/shallow_crit.jsonl)', shallow),
                    ('C in CBOX, tan b_A >= 0.3, a_y <= 1.5479: Theorem M (certs/main_cbox.jsonl)', maincbox)]
        return lst
    cls.preds = preds
    return cls
_add_old_regions(PP)

# ------------------------------------------------------------------------------------------------ FEM fallback for one far arm
from fractions import Fraction as _Fr
def _dyad_down(x, bits): return _Fr(math.floor(x*2**bits), 2**bits)
def _dyad_up(x, bits): return _Fr(math.ceil(x*2**bits), 2**bits)
def farl_setup(bf1, bf2, Cbox, k, RFv=None):
    """FARL chart (charts.farl_problem) for: far arm X = rho e(b), b in [bf1, bf2] subset (0, pi/2], rho >= RF, in the UPPER half
    plane (the arm 'A' of the chart), and the other triangle's apex in the Cartesian box Cbox = ((x0, x1), (y0, y1)) in the lower
    half plane.  Returns (problem, pc, rc) or None if the validity conditions cannot be verified."""
    from charts import farl_problem
    if not (0 < bf1 < bf2 <= math.pi/2): return None
    u_lo = fdown(1/arb(bf2).tan()); u_hi = fup(1/arb(bf1).tan())            # directed: cot is decreasing
    u0 = 0.5*(u_lo + u_hi)
    k0 = _Fr(round(1.0/(1.0 + u0*u0)*4096), 4096)   # RIGOR: heuristic: free dyadic choice of the cut slope; validity is checked rigorously below
    if not k0 > 0: return None                               # Lemma 9.6 needs k0 > 0 (k0 rounds to 0 for u0 > ~90.5)
    sin_lo = ival(bf1, bf2).sin().lower()                                     # exact arb point
    RA = _dyad_down(fdown(arb(RF if RFv is None else RFv)*sin_lo), 10)       # a_y = rho sin b >= RF sin b >= RA (rho >= RF)
    (x0, x1), (y0, y1) = Cbox
    # cut level: beyond the projections of B, D and the C-box along (u, 1)/(1+u^2), plus a margin
    need = 0.0
    for u in (u_lo, u_hi):
        for (cx, cy) in ((x0, y0), (x0, y1), (x1, y0), (x1, y1), (1.0, 0.0), (0.0, 0.0)):
            need = max(need, (cx*u + cy)/(1 + u*u))   # RIGOR: heuristic: choice of the cut level L; validity is checked rigorously below
    need = max(need, float(k0)*max(u_hi, 0.0))   # RIGOR: heuristic: choice of the cut level L; validity is checked rigorously below
    L = _dyad_up(need + 0.15, 8)
    if not (RA > L): return None
    # rigorous validity: B, D and the C-box corners strictly on B's side of the cut line for every u in [u_lo, u_hi]
    kk = arb(k0.numerator)/arb(k0.denominator); LL = arb(L.numerator)/arb(L.denominator)
    n = 32
    ucut = [u_lo + (u_hi - u_lo)*i/n for i in range(n)] + [u_hi]     # cells share endpoints and cover [u_lo, u_hi] exactly
    for i in range(n):
        uu = ival(ucut[i], ucut[i + 1])
        q1x, q1y = uu*LL, LL
        L2 = LL - kk*uu
        if not (L2.lower() > 0): return None
        q2x, q2y = 1 + uu*L2, L2
        dx, dy = q2x - q1x, q2y - q1y
        def side(px, py): return dx*(arb(py) - q1y) - dy*(arb(px) - q1x)
        sB = side(0.0, 0.0)
        if not (sB.lower() > 0 or sB.upper() < 0): return None
        sg = 1 if sB.lower() > 0 else -1
        for (px, py) in ((1.0, 0.0), (x0, y0), (x0, y1), (x1, y0), (x1, y1)):
            s_ = side(px, py)*sg
            if not (s_.lower() > 0): return None
        # A (on BA beyond Q1) must be on the other side: the point (RA)(u, 1) with RA > L
        sA = side_arb = dx*(arb(RA.numerator)/arb(RA.denominator) - q1y) - dy*(uu*arb(RA.numerator)/arb(RA.denominator) - q1x)
        if not ((sA*sg).upper() < 0): return None
    eu = enclose(u_lo, u_hi); ex = enclose(x0, x1); ey = enclose(y0, y1)       # exact-endpoint boxes containing the ranges
    pc = [eu[0], 0.0, ex[0], ey[0]]; rc = [eu[1], 0.0, ex[1], ey[1]]
    return farl_problem(k, RA, L, k0), pc, rc

def _polar_cart_lower(b1, b2, r1, r2):
    """Cartesian enclosure ((x0,x1),(y0,y1)) of {rho (cos b, -sin b)}, directed floats"""
    b = ival(b1, b2); rho = ival(r1, r2)
    x = afin(rho*b.cos()); y = afin(-rho*b.sin())
    return ((fdown(x), fup(x)), (fdown(y), fup(y)))

def _fp_fem_setup(self, p0, r, k):
    exact_box(p0, r)
    if self.Afar and not self.Cfar:
        bf1, bf2 = p0[0] - r[0], p0[0] + r[0]
        Cbox = _polar_cart_lower(p0[2] - r[2], p0[2] + r[2], exp_lo(p0[3] - r[3]), exp_hi(p0[3] + r[3]))
    elif self.Cfar and not self.Afar:
        # mirror (A, C) -> (C bar, A bar): the far C-arm becomes the far 'A' arm of the chart, A bar the other apex
        bf1, bf2 = p0[2] - r[2], p0[2] + r[2]
        Cbox = _polar_cart_lower(p0[0] - r[0], p0[0] + r[0], exp_lo(p0[1] - r[1]), exp_hi(p0[1] + r[1]))
    else:
        return None
    if Cbox[1][1] >= 0: return None
    return farl_setup(bf1, bf2, Cbox, k)
PF.fem_setup = _fp_fem_setup
