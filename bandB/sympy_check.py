"""sympy_check.py — symbolic re-derivation of every formula the band certificate relies on (BAND_SPEC §1, §3).
Author: Anonymous.  Run: ~/.venvs/analysis-py313/bin/python sympy_check.py   (exits non-zero on any failure)

1. the weighted 1-D form: for u = f(x) sin(pi (y - a(x))/h(x)), h = b - a,
     int_a^b |grad u|^2 dy = f'^2 h/2 + f f' h'/2 + f^2 (pi^2/h) (1/2 + a'b'/2 + h'^2 (1/6 + 1/(4 pi^2)))
     int_a^b u^2 dy       = f^2 h/2
   (general smooth a, b — pointwise in x);
2. the per-element closed forms used in cert.jl / cert_float.jl (P1 f, g; linear a, b on an element);
3. the f g / h integrals: log closed forms, series coefficients, end elements;
4. multilinearity of all vertex coordinates and all c_i in the chart coordinates; the c_i formulas of BAND_SPEC §1;
5. the inner tent minimum kink / apex, and positivity of the coefficient (1/2 + a'b'/2 + h'^2 C6) >= 1/2.
"""
import sys
import sympy as sp

ok = True
def check(name, expr):
    global ok
    e = sp.simplify(expr)
    good = (e == 0)
    print(("PASS " if good else "FAIL ") + name + ("" if good else f"   residual = {e}"))
    ok = ok and good

x, y, eta = sp.symbols('x y eta', real=True)
pi = sp.pi
f = sp.Function('f')(x); a = sp.Function('a')(x); b = sp.Function('b')(x)
h = b - a
u = f * sp.sin(pi * (y - a) / h)
grad2 = sp.diff(u, x)**2 + sp.diff(u, y)**2
# integrate over y in (a, b): substitute y = a + eta h, dy = h d eta
G = sp.simplify(grad2.subs(y, a + eta * h) * h)
lhs = sp.integrate(sp.expand(G), (eta, 0, 1))
fp, ap, bp = sp.diff(f, x), sp.diff(a, x), sp.diff(b, x)
hp = bp - ap
C6 = sp.Rational(1, 6) + 1 / (4 * pi**2)
rhs = fp**2 * h / 2 + f * fp * hp / 2 + f**2 * (pi**2 / h) * (sp.Rational(1, 2) + ap * bp / 2 + hp**2 * C6)
check("1a  int |grad u|^2 dy  = weighted form N", sp.simplify(lhs - rhs))
lhsD = sp.integrate(sp.expand((u**2).subs(y, a + eta * h) * h), (eta, 0, 1))
check("1b  int u^2 dy = f^2 h / 2", sp.simplify(lhsD - f**2 * h / 2))
# the four elementary integrals quoted in the spec
check("1c  int sin^2 = 1/2", sp.integrate(sp.sin(pi * eta)**2, (eta, 0, 1)) - sp.Rational(1, 2))
check("1d  int eta sin(pi eta) pi cos(pi eta) = -1/4", sp.integrate(eta * sp.sin(pi * eta) * pi * sp.cos(pi * eta), (eta, 0, 1)) + sp.Rational(1, 4))
check("1e  int eta cos^2 = 1/4", sp.integrate(eta * sp.cos(pi * eta)**2, (eta, 0, 1)) - sp.Rational(1, 4))
check("1f  int eta^2 cos^2 = 1/6 + 1/(4 pi^2)", sp.integrate(eta**2 * sp.cos(pi * eta)**2, (eta, 0, 1)) - C6)

# 2. element formulas: s in [0,1], x = xa + L s; f, g, a, b linear
s, L = sp.symbols('s L', positive=True)
fa, fb, ga, gb, aa, ab, ba, bb = sp.symbols('f_a f_b g_a g_b a_a a_b b_a b_b', real=True)
F = fa + (fb - fa) * s; Gf = ga + (gb - ga) * s
A = aa + (ab - aa) * s; Bf = ba + (bb - ba) * s; Hs = Bf - A
ha, hb = ba - aa, bb - ab
dF, dG = (fb - fa) / L, (gb - ga) / L
dA, dB = (ab - aa) / L, (bb - ba) / L
dH = dB - dA
ce = sp.Rational(1, 2) + dA * dB / 2 + dH**2 * C6
# bilinear form N(f,g) = int f'g' h/2 + (f'g + f g') h'/4 + f g (pi^2/h) ce ; D(f,g) = int f g h / 2
K1 = sp.integrate(dF * dG * Hs / 2 * L, (s, 0, 1))
check("2a  int f'g'h/2 = (fb-fa)(gb-ga)(ha+hb)/(4L)", K1 - (fb - fa) * (gb - ga) * (ha + hb) / (4 * L))
K2 = sp.integrate((dF * Gf + F * dG) * dH / 4 * L, (s, 0, 1))
check("2b  int (f'g+fg')h'/4 = (h'/4)(fb gb - fa ga)", K2 - dH / 4 * (fb * gb - fa * ga))
M = sp.integrate(F * Gf * Hs / 2 * L, (s, 0, 1))
Mform = L / 24 * (fa * ga * (3 * ha + hb) + (fa * gb + fb * ga) * (ha + hb) + fb * gb * (ha + 3 * hb))
check("2c  int f g h/2 = L/24 [..]", M - Mform)
# f g / h: int_0^1 F G / H ds = fa ga I00 + (fa gb + fb ga) I01 + fb gb I11
hA, hB = sp.symbols('h_A h_B', positive=True)
Hlin = hA * (1 - s) + hB * s
W00, W01, W11 = (1 - s)**2, s * (1 - s), s**2
D_ = hB - hA; lg = sp.log(hB / hA)
I11c = (D_ * (hB - 3 * hA) / 2 + hA**2 * lg) / D_**3
I00c = (D_ * (hA - 3 * hB) / 2 + hB**2 * lg) / D_**3
I01c = (D_ - hA * lg) / D_**2 - I11c
for nm, W, cf in (("I00", W00, I00c), ("I01", W01, I01c), ("I11", W11, I11c)):
    # exact integral at rational (hA, hB) vs the closed form, compared to 40 digits
    bad = 0
    for (p, q) in ((1, 3), (2, 7), (5, 2), (9, 4), (1, 10)):
        v1 = sp.N(sp.integrate(W / (p * (1 - s) + q * s), (s, 0, 1)), 40); v2 = sp.N(cf.subs({hA: p, hB: q}), 40)
        bad += int(bool(abs(v1 - v2) > 1e-30))
    check(f"3a  closed form {nm} (log) at 5 rational points", sp.Integer(bad))
check("3b  F G/H integrand expansion", sp.simplify(sp.expand(F * Gf) - (fa * ga * (1 - s)**2 + (fa * gb + fb * ga) * s * (1 - s) + fb * gb * s**2)))
# series: 1/h = (1/hm) sum (-eps sigma)^n,  sigma = s - 1/2 ; coefficients C_n = int sigma^n (1/2+sigma)^2, B_n = int sigma^n (1/4 - sigma^2)
sig = sp.symbols('sigma', real=True)
for n in range(0, 9):
    Cn = sp.integrate(sig**n * (sp.Rational(1, 2) + sig)**2, (sig, -sp.Rational(1, 2), sp.Rational(1, 2)))
    Bn = sp.integrate(sig**n * (sp.Rational(1, 4) - sig**2), (sig, -sp.Rational(1, 2), sp.Rational(1, 2)))
    An = sp.integrate(sig**n * (sp.Rational(1, 2) - sig)**2, (sig, -sp.Rational(1, 2), sp.Rational(1, 2)))
    mom = lambda k: sp.Rational(1, 2**k * (k + 1)) if k % 2 == 0 else 0
    check(f"3c  series coefficients n={n}: C_n, B_n, A_n = (-1)^n C_n",
          sp.Abs(Cn - (mom(n) / 4 + mom(n + 1) + mom(n + 2))) + sp.Abs(Bn - (mom(n) / 4 - mom(n + 2))) + sp.Abs(An - (-1)**n * Cn))
# the series reproduces I00, I01, I11 (truncated at order 12 vs exact, eps = 1/5)
epsv = sp.Rational(1, 5); hm = 1
ser = lambda w: sum(sp.integrate(w * (-epsv * sig)**n, (sig, -sp.Rational(1, 2), sp.Rational(1, 2))) for n in range(0, 40)) / hm
ex = lambda w: sp.integrate(w / (hm * (1 + epsv * sig)), (sig, -sp.Rational(1, 2), sp.Rational(1, 2)))
for nm, w in (("I00", (sp.Rational(1, 2) - sig)**2), ("I01", sp.Rational(1, 4) - sig**2), ("I11", (sp.Rational(1, 2) + sig)**2)):
    check(f"3d  series of {nm} converges to the exact value (40 terms, |err| < 1e-25)", sp.Integer(int(bool(abs(sp.N(ser(w) - ex(w), 50)) > 1e-25))))
# end elements: h(0) = 0 with f(0) = 0 => int_0^L f g / h dx = L fb gb /(2 hb)
check("3e  end element [0,L]: int f g/h = L fb gb/(2 hb)", sp.integrate((fb * s) * (gb * s) / (hB * s) * L, (s, 0, 1)) - L * fb * gb / (2 * hB))
check("3f  end element [.,1]: int f g/h = L fa ga/(2 ha)", sp.integrate((fa * (1 - s)) * (ga * (1 - s)) / (hA * (1 - s)) * L, (s, 0, 1)) - L * fa * ga / (2 * hA))

# 4. charts: vertices and c_i, multilinearity
y3, t, x1, x3, y2, x2 = sp.symbols('y3 t x1 x3 y2 x2', real=True)
def cvals(P):
    return [sp.expand((P[(i + 1) % 4][0] - P[i][0]) * (P[(i + 2) % 4][1] - P[i][1]) - (P[(i + 1) % 4][1] - P[i][1]) * (P[(i + 2) % 4][0] - P[i][0])) for i in range(4)]
def multilinear(e, vs):
    pe = sp.Poly(sp.expand(e), *vs)
    return all(max(m) <= 1 for m in pe.monoms()) if pe.monoms() else True
P1 = [(0, 0), (x1, -t * y3), (1, 0), (x3, y3)]
c = cvals(P1)
check("4a  chart I: c0 = y1, c2 = y3", sp.Abs(c[0] - t * y3) + sp.Abs(c[2] - y3))
check("4b  chart I: c1 = (P2-P1)x(P3-P1) = y3(1 + t - x1 - t x3)", c[1] - y3 * (1 + t - x1 - t * x3))
check("4c  chart I: c3 = (P0-P3)x(P1-P3) = y3(x1 + t x3)", c[3] - y3 * (x1 + t * x3))
v1 = [y3, t, x1, x3]
allml = all(multilinear(e, v1) for e in c) and all(multilinear(P1[j][k] - P1[i][k], v1) for i in range(4) for j in range(4) for k in range(2))
allml = allml and all(multilinear(c[i] - sp.Rational(7, 200) * c[(i + 2) % 4], v1) for i in range(4))
check("4d  chart I: every c_i, c_i - (7/200)c_{i+2} and vertex-difference component is multilinear", sp.Integer(0 if allml else 1))
P2 = [(0, 0), (1, 0), (x2, y2), (x3, t * y2)]
c = cvals(P2)
check("4e  chart II: c0 = y2, c3 = y3", sp.Abs(c[0] - y2) + sp.Abs(c[3] - t * y2))
check("4f  chart II: c1 = y2(1-x3) - y3(1-x2)", c[1] - (y2 * (1 - x3) - t * y2 * (1 - x2)))
check("4g  chart II: c2 = x2 y3 - x3 y2", c[2] - (x2 * t * y2 - x3 * y2))
v2 = [y2, t, x2, x3]
allml = all(multilinear(e, v2) for e in c) and all(multilinear(P2[j][k] - P2[i][k], v2) for i in range(4) for j in range(4) for k in range(2))
allml = allml and all(multilinear(c[i] - sp.Rational(7, 200) * c[(i + 2) % 4], v2) for i in range(4))
check("4h  chart II: every c_i, c_i - (7/200)c_{i+2} and vertex-difference component is multilinear", sp.Integer(0 if allml else 1))
# width normal to edge i: distances of P_{i+2}, P_{i+3} to line(P_i P_{i+1}) are c_i/|e_i| and c_{i+3}/|e_i|
for P, nm in ((P1, "I"), (P2, "II")):
    c = cvals(P); res = 0
    for i in range(4):
        e = (P[(i + 1) % 4][0] - P[i][0], P[(i + 1) % 4][1] - P[i][1])
        for k, ci in ((2, c[i]), (3, c[(i + 3) % 4])):
            w = (P[(i + k) % 4][0] - P[i][0], P[(i + k) % 4][1] - P[i][1])
            res += sp.Abs(sp.expand(e[0] * w[1] - e[1] * w[0] - ci))
    check(f"4i  chart {nm}: e_i x (P_(i+2)-P_i) = c_i and e_i x (P_(i+3)-P_i) = c_(i+3)", res)

# 5. tents
X, l, uu = sp.symbols('X l u', positive=True)
k = uu / (1 + uu - l)
check("5a  T-(x;l,u) = min(x/u, (1-x)/(1-l)) has its kink at u/(1+u-l)", sp.simplify(k / uu - (1 - k) / (1 - l)))
ylo = sp.symbols('ylo', positive=True)
check("5b  apex height ylo*k/u = ylo/(1+u-l)", sp.simplify(ylo * k / uu - ylo / (1 + uu - l)))
p_, q_ = sp.symbols('p q', real=True)
Qf = sp.expand(p_ * q_ / 2 + (q_ - p_)**2 * C6)    # a'b'/2 + h'^2 C6 with a' = p, b' = q
Mq = sp.Matrix([[C6, sp.Rational(1, 4) - C6], [sp.Rational(1, 4) - C6, C6]])
check("5c  a'b'/2 + h'^2 C6 is the quadratic form of [[C6, 1/4-C6],[1/4-C6, C6]]", sp.expand(Qf - (sp.Matrix([p_, q_]).T * Mq * sp.Matrix([p_, q_]))[0]))
check("5d  that form is positive semidefinite (C6 > 0, det = C6^2-(1/4-C6)^2 = (2C6-1/4)/4 > 0) => ce >= 1/2",
      sp.Integer(0 if (bool(sp.N(C6) > 0) and bool(sp.N(Mq.det()) > 0)) else 1))
print("ALL PASS" if ok else "SOME CHECKS FAILED")
sys.exit(0 if ok else 1)
