"""Rigorous verification for Lemma 5.3 (thin wedge at B, sector version): for every alpha in [0.0997, 0.1772],
   xi(Q) <= (j_{pi/alpha_k, 2} * tan(alpha_{k+1}) / pi)^2 < 5/2   on each grid interval [alpha_k, alpha_{k+1}],
alpha_k = 0.0997 + k * 0.0025 (k = 0..31, exact decimals; alpha_31 = 0.1772 exactly).
Upper bounds X_k > j_{nu,2} (nu = pi/alpha_k) are certified by sign checks in ball arithmetic: J_nu(Y) < 0 and J_nu(X) > 0 with
0 < Y < X checked explicitly (J_nu > 0 on (0, j_{nu,1}), so a sign change in (Y, X) after a negative value at Y puts j_{nu,2}
below X), and X < 2 nu < j_{2 nu, 1} (so lambda_2 of the sector is j_{nu,2}^2/R^2).  The float zeros from mpmath only choose
X and Y; every decision is a ball-arithmetic comparison."""
from fractions import Fraction as Fr
from flint import arb, ctx
import mpmath as mp
ctx.prec = 128
A0, STEP, NSTEP = Fr('0.0997'), Fr('0.0025'), 31
assert A0 + NSTEP*STEP == Fr('0.1772')
pi = arb.pi()
def q2a(q): return arb(q.numerator)/arb(q.denominator)
def J(nu, x): return arb(x).bessel_j(nu)
ok = True; bounds = []
for k in range(NSTEP):
    a = q2a(A0 + k*STEP); a_next = q2a(A0 + (k + 1)*STEP)
    nu = pi/a
    nuf = float(nu.mid())                                  # heuristic only: choose the test points X, Y   # RIGOR: heuristic: chooses the test points X, Y only
    j1 = float(mp.besseljzero(nuf, 1)); j2 = float(mp.besseljzero(nuf, 2))   # RIGOR: heuristic: chooses the test points X, Y only
    Y = arb((j1 + j2)/2); X = arb(j2*(1 + 1e-6))
    cY = J(nu, Y); cX = J(nu, X)
    good = bool(Y > 0) and bool(Y < X) and bool(cY < 0) and bool(cX > 0) and bool(X < 2*nu)
    t = X*a_next.tan()/pi
    bound = t*t
    good = good and bool(bound < arb(5)/2)
    ok = ok and good; bounds.append(bound)
    print(f"alpha in [{float(A0 + k*STEP):.4f}, {float(A0 + (k + 1)*STEP):.4f}]: nu={nuf:.3f} X={float(X.mid()):.5f} "   # RIGOR: info: printed
          f"0<Y<X:{bool(Y > 0) and bool(Y < X)} J(Y)<0:{bool(cY < 0)} J(X)>0:{bool(cX > 0)} bound<={float(bound.upper()):.9f} {'OK' if good else 'FAIL'}")   # RIGOR: info: printed
from nanguard import amax
worst = amax(*bounds)
print('ALL OK' if ok else 'FAILED', 'worst bound <=', worst)
import sys; sys.exit(0 if ok else 1)
