"""Fail-closed, DIRECTED conversion of ball-arithmetic results to floats (adversarial review, 2026-10-02).
 * python-flint can return NaN (e.g. ball**k for a ball whose midpoint is 0), and Python's max()/min() silently drop a NaN that
   comes second: every value entering a rigorous reduction or decision is checked (NaNError on non-finite);
 * float(arb) rounds to NEAREST, which can move a bound inward by up to half an ulp: a lower bound is converted with fdown()
   (largest float <= the exact lower endpoint), an upper bound with fup() (smallest float >= the exact upper endpoint).
Certificate callers treat any exception as failure; family predicates treat NaNError as 'not established'."""
import math
from flint import arb, ctx
if ctx.prec < 106: ctx.prec = 106          # the working precision of every verifier module
class NaNError(ArithmeticError): pass
def afin(a):
    """the arb ball a, raising NaNError unless it is finite"""
    if not a.is_finite(): raise NaNError(f'non-finite ball {a}')
    return a
def fin(x):
    """float(x) for a value that is used only as a heuristic / non-decisive number, raising NaNError unless finite"""
    v = float(x)   # RIGOR: helper: fin() is for non-decisive values only (every call site is annotated)
    if not math.isfinite(v): raise NaNError(f'non-finite value {x!r}')
    return v
def fdown(x):
    """a float <= every point of x (x: arb ball, or a number): the exact lower endpoint rounded DOWN; NaNError if non-finite"""
    a = x if isinstance(x, arb) else arb(x)
    afin(a)
    e = a.lower()                      # exact (radius 0)
    v = float(e)   # RIGOR: helper: nearest float of the exact lower endpoint, corrected DOWN by the exact comparison below
    if not math.isfinite(v): raise NaNError(f'non-finite endpoint {e}')
    if arb(v) > e: v = math.nextafter(v, -math.inf)   # RIGOR: helper: one ulp down when the nearest float lies above the endpoint
    return v
def fup(x):
    """a float >= every point of x: the exact upper endpoint rounded UP; NaNError if non-finite"""
    a = x if isinstance(x, arb) else arb(x)
    afin(a)
    e = a.upper()
    v = float(e)   # RIGOR: helper: nearest float of the exact upper endpoint, corrected UP by the exact comparison below
    if not math.isfinite(v): raise NaNError(f'non-finite endpoint {e}')
    if arb(v) < e: v = math.nextafter(v, math.inf)   # RIGOR: helper: one ulp up when the nearest float lies below the endpoint
    return v
def amax(*xs):
    """max of arb balls as an arb upper bound (ball whose every point >= each x); all must be finite"""
    best = None
    for x in xs:
        u = afin(x if isinstance(x, arb) else arb(x)).upper()
        best = u if best is None or u > best else best
    return arb(best)

# ---------------------------------------------------------------------------------------------- exact boxes and subdivisions
from fractions import Fraction as _Fr
class InexactBox(ArithmeticError): pass
def exact_box(p0, r):
    """fail closed (InexactBox) unless, for every coordinate, the float expressions p0_i - r_i and p0_i + r_i are EXACT (so that
    every predicate / certificate that forms them works with the true endpoints of the represented closed box) and r_i >= 0"""
    for c, h in zip(p0, r):
        c = float(c); h = float(h)   # RIGOR: helper: float() of float values; exactness is checked with Fractions below
        if not (math.isfinite(c) and math.isfinite(h) and h >= 0): raise InexactBox(f'bad box coordinate {c!r} +- {h!r}')
        C = _Fr(c); H = _Fr(h)
        if _Fr(c - h) != C - H or _Fr(c + h) != C + H: raise InexactBox(f'box endpoints of {c!r} +- {h!r} are not exact floats')
    return True
def cuts(lo, hi, n):
    """n consecutive cells [x_i, x_{i+1}] whose union is [lo, hi]: x_0 = lo and x_n = hi EXACTLY, consecutive cells share the
    same float endpoint (interior points may be rounded; any values in between are fine)"""
    return [lo] + [lo + (hi - lo)*i/n for i in range(1, n)] + [hi]
def enclose(lo, hi):
    """floats (c, h), h >= 0, such that c - h and c + h are EXACT floats (no rounding in the float expressions) and the represented
    interval [c - h, c + h] contains [lo, hi] (lo <= hi: floats or arb balls; lo rounded down, hi up): the ends are snapped
    outward to a grid of spacing q = 2^k ulp(max |end|) (k = 0, 1, ... until c, h, c - h, c + h are all exact), all checked in
    rational arithmetic"""
    lo_f = fdown(lo); hi_f = fup(hi)
    if not lo_f <= hi_f: raise NaNError(f'enclose: empty interval {lo!r} {hi!r}')
    L = _Fr(lo_f); H = _Fr(hi_f)
    M = max(abs(lo_f), abs(hi_f), 2.0**-1000)   # RIGOR: helper: grid scale (a free choice; containment is checked exactly)
    q = _Fr(math.ulp(M))   # RIGOR: helper: grid spacing (a power of two)
    for _ in range(12):
        a = (L/q).__floor__()*q; b = (-((-H)/q).__floor__())*q
        c = (a + b)/2; h = (b - a)/2
        cf, hf = float(c), float(h)   # RIGOR: helper: conversions of rationals, accepted only if exact (checked)
        if _Fr(cf) == c and _Fr(hf) == h and _Fr(cf - hf) == a and _Fr(cf + hf) == b and a <= L and b >= H:
            return cf, hf
        q = 2*q
    raise NaNError(f'enclose failed for {lo!r} {hi!r}')
def dec(x):
    """an arb ball containing the EXACT decimal threshold x (a decimal string as written in a reason / lemma, an int, or a float
    LITERAL, read as the decimal of its shortest repr: 1.8 -> 9/5).  Thresholds are compared with directed endpoints of this ball."""
    if isinstance(x, arb): return x
    q = decq(x)
    return arb(q.numerator)/arb(q.denominator)
def decq(x):
    """the exact rational of a decimal threshold (as dec)"""
    if isinstance(x, (str, int)): return _Fr(x)
    return _Fr(repr(float(x)))   # RIGOR: helper: definition: a float LITERAL denotes the decimal of its shortest repr
def amin(*xs):
    """min of arb balls as an arb LOWER bound (exact point <= every point of each x); all must be finite"""
    best = None
    for x in xs:
        l = afin(x if isinstance(x, arb) else arb(x)).lower()
        best = l if best is None or l < best else best
    return arb(best)
def upper_min(*xs):
    """for a quantity bounded above by EACH of the balls xs: the smallest of their upper ends (exact arb point)"""
    best = None
    for x in xs:
        u = afin(x if isinstance(x, arb) else arb(x)).upper()
        best = u if best is None or u < best else best
    return arb(best)
def lower_max(*xs):
    """for a quantity bounded below by EACH of the balls xs: the largest of their lower ends (exact arb point)"""
    best = None
    for x in xs:
        l = afin(x if isinstance(x, arb) else arb(x)).lower()
        best = l if best is None or l > best else best
    return arb(best)
