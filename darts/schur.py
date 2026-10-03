"""Rigorous negativity test of a symmetric 2x2 matrix of Taylor models (the Schur-complement step of gcert / p2up / p2fix).
N(p) < 0 for every p of the box if  N11(p) <= n11 < 0  and  N22(p) + N12(p)^2/|N11(p)| <= N22.upper + n12m^2/|n11| < 0,
where n11 = upper end of the range of N11 and n12m = max(|lo|, |hi|) of the range of N12 -- all in ball arithmetic (no float
conversion; products, not powers)."""
from nanguard import afin, amax

def schur_upper(N):
    """returns (ok, schur, why); schur is an arb upper bound of max over the box of N22 + N12^2/|N11| (when N11 < 0 on the box)"""
    n11 = N[0][0].upper()
    if not (n11 < 0): return False, None, 'N11 >= 0'
    lo, hi = N[0][1].range()
    n12m = amax(abs(afin(lo)), abs(afin(hi)))
    schur = afin(N[1][1].upper() + n12m*n12m/abs(n11))
    if not (schur < 0): return False, schur, 'schur >= 0'
    return True, schur, None
