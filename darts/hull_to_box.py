"""Rigorous check that a set of darts described by hull parameters lies in the BD box of Theorem H (README §2).

Hull placement C=(0,0), A=(1,0), B=(v, 1)/(u+v) with u = cot θ_A, v = cot θ_C (θ_A, θ_C = hull angles at A, C);
D = αA + βB + γC (barycentric; β = notch depth), s := α/(α+γ) = position on AC of the foot of the cevian BD (|CF| : |FA| = s : 1-s).
BD chart (orientation preserving similarity B -> 0, D -> 1):  a = (A-B)/(D-B), c = (C-B)/(D-B).  Since D - B = (1-β)(s - B),
with p := s u - (1-s) v:
   a_x = (u p + 1)/((p²+1)(1-β)),  a_y = (1-s)(u+v)/((p²+1)(1-β)),  c_x = (1 - v p)/((p²+1)(1-β)),  c_y = -s(u+v)/((p²+1)(1-β)).
These are enclosed over parameter sub-boxes by order-3 Taylor models with rigorous remainders (tm.py, arb 106 bits); a sub-box is
accepted only if all four enclosures lie in BOX; otherwise it is bisected; the run aborts if a single POINT maps outside BOX.
The mirror x -> 1-x maps (θ_A, θ_C, s) to (θ_C, θ_A, 1-s) and preserves ξ, so for θ-ranges symmetric in (A, C) one half in s suffices.
Usage: hull_to_box.py TH_LO TH_HI (degrees, for both θ_A and θ_C) BETA_LO BETA_HI S_LO S_HI"""
import sys
from flint import arb
from tm import TM
BOX = ((0.965, 1.262), (0.085, 0.796), (0.479, 1.213), (-1.091, -0.483))
def coords(u, v, s, be):
    p = s*u - (1 - s)*v
    F = (1 - be).recip() if isinstance(be, TM) else 1/(1 - be)
    w = (p*p + 1).recip() if isinstance(p, TM) else 1/(p*p + 1)
    g = w*F
    return [(u*p + 1)*g, (1 - s)*(u + v)*g, (1 - v*p)*g, -(s*(u + v))*g]
class Outside(Exception): pass
def check(b, depth, st):
    from nanguard import enclose
    enc = [enclose(b[i][0], b[i][1]) for i in range(4)]        # exact-endpoint boxes containing [b_lo, b_hi]
    tms = [TM.var(i, enc[i][0], enc[i][1]) for i in range(4)]   # RIGOR: exact: TM.var builds the exact variable c + h t (not a NumPy variance)
    try:
        img = [f.range() for f in coords(*tms)]
        good = all(lo.lower() >= BOX[i][0] and hi.upper() <= BOX[i][1] for i, (lo, hi) in enumerate(img))
    except ValueError:
        good = False
    if good: st[0] += 1; return True
    mid = coords(*[arb((x[0] + x[1])/2) for x in b])
    if any(m.upper() < BOX[i][0] or m.lower() > BOX[i][1] for i, m in enumerate(mid)):
        raise Outside(f'parameter point {[(x[0] + x[1])/2 for x in b]} maps to {[float(m.mid()) for m in mid]}')   # RIGOR: info: failure message
    if depth > 40: raise Outside(f'no decision at depth 40 near {b}')
    w = [(b[i][1] - b[i][0])/SCALE[i] for i in range(4)]; i = w.index(max(w)); m = (b[i][0] + b[i][1])/2   # RIGOR: heuristic: split choice (children share the float midpoint exactly)
    b1 = list(b); b2 = list(b); b1[i] = (b[i][0], m); b2[i] = (m, b[i][1])
    return check(b1, depth + 1, st) and check(b2, depth + 1, st)
if __name__ == '__main__':
    th0, th1, be0, be1, s0, s1 = [float(x) for x in sys.argv[1:7]]   # RIGOR: exact: the parameter box is defined by these floats
    pi = arb.pi()
    def cot(deg): return 1/(arb(deg)*pi/180).tan()
    # u, v range: [cot θ_hi, cot θ_lo] rounded OUTWARD (cot is decreasing on (0, π/2))
    from nanguard import fdown, fup
    ulo = fdown(cot(th1)); uhi = fup(cot(th0))                 # directed
    box = [(ulo, uhi), (ulo, uhi), (s0, s1), (be0, be1)]
    SCALE = [uhi - ulo, uhi - ulo, s1 - s0, be1 - be0]
    st = [0]
    try: ok = check(box, 0, st); why = ''
    except Outside as e: ok = False; why = str(e)
    print(dict(theta_deg=(th0, th1), beta=(be0, be1), s=(s0, s1), u_v_range=(ulo, uhi), contained=ok, subboxes=st[0], why=why))
    sys.exit(0 if ok else 1)
