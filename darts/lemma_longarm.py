"""Lemma 9.5 (very long steep arm) — certified table, ball arithmetic (arb, 128 bits).
BD chart, |BD| = 1, L := |BA|, b := b_A.  For 0 < l < L - 1:
   xi(Q) <= (L+1)^2 / (sin^2 b (L-1-l)^2) + 4/l^2        (LEMMAS.md §9, Lemma 9.5)
The right-hand side is decreasing in L (fixed l) and depends on b only through sin b, which on [b_k, pi - b_k] is >= sin b_k.
For each b_k we choose (float search) l_k and L_k, then certify RHS(L_k, l_k, b_k) < 5/2 in arb; hence every dart with
b_A in [b_k, pi - b_k] and |BA| >= L_k |BD| has xi < 5/2 (and the same with A <-> C)."""
import math
from flint import arb, ctx
ctx.prec = 128
pi = arb.pi()
def rhs(L, l, b_deg):
    s = (arb(b_deg)*pi/180).sin()
    q = arb(L) - 1 - arb(l); p = arb(L) + 1
    return p*p/(s*s*q*q) + 4/(arb(l)*arb(l))
def Lmin(b_deg, l):
    q = math.sin(math.radians(b_deg))*math.sqrt(max(2.5 - 4/l**2, 0))   # RIGOR: heuristic: float search for the table parameters (certified below in arb)
    return (1 + q*(1 + l))/(q - 1) if q > 1 else math.inf
ok = True
for b in [40, 41, 42, 43, 45, 47.5, 50, 55, 60, 65, 70, 80, 90]:
    l = min((x/10 for x in range(13, 2000)), key=lambda x: Lmin(b, x))   # RIGOR: heuristic: float search for the table parameters
    L = math.ceil(Lmin(b, l)*1.001 + 0.01)
    v = rhs(L, l, b); good = v.upper() < 2.5; ok &= bool(good)
    print(f"b_A in [{b:5.1f} deg, {180 - b:5.1f} deg]:  |BA| >= {L:6d} |BD|  (l = {l:5.1f}): xi <= {float(v.upper()):.5f}  {'OK' if good else 'FAIL'}")   # RIGOR: info: printed
th = (arb(2)/5).sqrt().asin()*180/pi
print('necessary for this bound: sin^2 b_A > 2/5, i.e. b_A in (', th, ',', 180 - th, ') deg')
print('ALL OK' if ok else 'SOME FAIL')
import sys; sys.exit(0 if ok else 1)
