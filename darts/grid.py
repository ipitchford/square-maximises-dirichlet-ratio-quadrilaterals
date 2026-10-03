"""Exact dyadic snapping of root boxes (used by the replays, the accounting and the production drivers)."""
def dyadic_root(p0, r, bits=16):
    """snap a box OUTWARD to the dyadic grid 2^-bits, in exact rational arithmetic: lo -> floor, hi -> ceil of the EXACT ends
    p0 -/+ r of the represented box.  Then every bisection p +- r/2 is exact in binary floating point (for depth <= ~30), so the
    leaves of a root tile it exactly (checked by the replays).  The returned floats are exact (checked)."""
    import math
    from fractions import Fraction as Fr
    sc = 2**bits
    p0n = []; rn = []
    for c, h in zip(p0, r):
        c = float(c); h = float(h)   # RIGOR: exact: float() of the float box values
        if h == 0: p0n.append(c); rn.append(0.0); continue
        if not h > 0: raise ValueError('negative radius')
        lo = Fr(math.floor((Fr(c) - Fr(h))*sc), sc); hi = Fr(math.ceil((Fr(c) + Fr(h))*sc), sc)
        cm = (lo + hi)/2; hm = (hi - lo)/2
        if Fr(float(cm)) != cm or Fr(float(hm)) != hm: raise ValueError('dyadic_root: snapped box not exactly representable')   # RIGOR: exact: the conversion is checked
        p0n.append(float(cm)); rn.append(float(hm))   # RIGOR: exact: checked exactly representable above
    return p0n, rn
