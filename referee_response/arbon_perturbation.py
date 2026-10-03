"""Exact first-order perturbation of the square's eigenvalues under single-vertex moves (Hadamard formula), and the
failure of the rule 'first-order shift of a double eigenvalue under a sum of vertex moves = sum of the separate minimal
shifts' (Arbon 2022, eq. (3.16)) along the isosceles-trapezoid direction.  Exact arithmetic (sympy).
Unit square, vertices v1=(0,0), v2=(1,0), v3=(1,1), v4=(0,1) (counter-clockwise). Dirichlet eigenpairs:
lambda_1 = 2 pi^2, u0 = 2 sin(pi x) sin(pi y);  lambda_2 = lambda_3 = 5 pi^2, u1 = 2 sin(pi x) sin(2 pi y), u2 = 2 sin(2 pi x) sin(pi y).
Hadamard: under a boundary motion with outward normal speed V, dlambda = -int (d_n u)^2 V ds for a simple eigenvalue, and the
first-order splitting of the double eigenvalue is given by the eigenvalues of M_ij = -int d_n u_i d_n u_j V ds.
Moving one vertex by delta (others fixed) gives V = (hat function, 1 at that vertex, 0 at the neighbours) * (delta . n) on its two
edges."""
import sympy as sp
from sympy.simplify.fu import TR8
sp.TR8 = TR8
x, y, s = sp.symbols('x y s', real=True)
pi = sp.pi
u = [2*sp.sin(pi*x)*sp.sin(pi*y), 2*sp.sin(pi*x)*sp.sin(2*pi*y), 2*sp.sin(2*pi*x)*sp.sin(pi*y)]
V = [sp.Matrix([0, 0]), sp.Matrix([1, 0]), sp.Matrix([1, 1]), sp.Matrix([0, 1])]
grad = [sp.Matrix([sp.diff(f, x), sp.diff(f, y)]) for f in u]

def integ(expr):
    """exact int_0^1 expr ds for (polynomial in s) x (trigonometric polynomial in pi s): product-to-sum, then termwise"""
    e = sp.expand(sp.fu(sp.expand(expr), measure=lambda z: 0) if False else sp.expand(sp.TR8(sp.expand(expr))))
    return sp.simplify(sum(sp.integrate(t, (s, 0, 1)) for t in sp.Add.make_args(e)))

def vertex_matrices(k, delta):
    """(dlambda_1, 2x2 M) for moving vertex k by delta (sympy 2-vector), exact"""
    d1 = 0; M = sp.zeros(2, 2)
    for (P, Q, w) in ((k, (k + 1) % 4, 1 - s), ((k - 1) % 4, k, s)):     # edge leaving v_k (weight 1-s), edge entering v_k (weight s)
        A, B = V[P], V[Q]; t = B - A
        n = sp.Matrix([t[1], -t[0]])                                      # outward normal of a CCW unit edge
        pt = {x: A[0] + s*t[0], y: A[1] + s*t[1]}
        dn = [sp.simplify((g.T*n)[0].subs(pt)) for g in grad]
        Vn = w*(delta.T*n)[0]
        d1 += -integ(dn[0]**2*Vn)
        for i in range(2):
            for j in range(2):
                M[i, j] += -integ(dn[1+i]*dn[1+j]*Vn)
    return sp.simplify(d1), sp.simplify(M)

def eo_derivative(mode):
    """vertex velocities d v_k / d(mode) at p = 0 for the Endo-Osting parametrisation (paper section 2.1)"""
    a, b, c, d = sp.symbols('a b c d')
    v = [sp.Matrix([0, 0]) + sp.Rational(1, 2)*sp.Matrix([a + d - b, -a + d + c]),
         sp.Matrix([1, 0]) + sp.Rational(1, 2)*sp.Matrix([-a + d + b, -a - d - c]),
         sp.Matrix([1, 1]) + sp.Rational(1, 2)*sp.Matrix([-a - d - b, a - d + c]),
         sp.Matrix([0, 1]) + sp.Rational(1, 2)*sp.Matrix([a - d + b, a + d - c])]
    sym = dict(a=a, b=b, c=c, d=d)[mode]
    return [vk.diff(sym) for vk in v]

def mineig(M):
    ev = list(M.eigenvals().keys())
    return sp.simplify(sp.Min(*ev)) if len(ev) > 1 else sp.simplify(ev[0])

out = {}
for mode in "abcd":
    vel = eo_derivative(mode)
    parts = [vertex_matrices(k, vel[k]) for k in range(4)]
    d1 = sp.simplify(sum(p[0] for p in parts)); H = sp.simplify(sum((p[1] for p in parts), sp.zeros(2, 2)))
    out[mode] = (vel, parts, d1, H)
    print(f"mode {mode}: dlambda_1 = {d1};  H = {H.tolist()}")
# the trapezoid direction b: true first-order slope of xi vs the sum-of-minima rule
vel, parts, d1, H = out['b']
lam1, xi = 2*pi**2, sp.Rational(5, 2)
print("\nb-direction (isosceles trapezoid with parallel sides 1+b and 1-b):")
for k, (d1k, Mk) in enumerate(parts):
    print(f"  vertex v{k+1} velocity {list(vel[k])}: dlambda_1 = {sp.nsimplify(d1k)}, M = {Mk.tolist()}, min eig = {mineig(Mk)}")
true_slope = sp.simplify((mineig(H) - xi*d1)/lam1)
rule4 = sp.simplify(sum(mineig(Mk) - xi*d1k for d1k, Mk in parts)/lam1)
top = parts[2][1] + parts[3][1]; bot = parts[0][1] + parts[1][1]
d1top = parts[2][0] + parts[3][0]; d1bot = parts[0][0] + parts[1][0]
rule2 = sp.simplify((mineig(top) - xi*d1top + mineig(bot) - xi*d1bot)/lam1)
print(f"  true first-order slope d xi/db at 0 = {true_slope}")
print(f"  sum over the 4 vertices of separate first-order slopes = {rule4} = {sp.N(rule4, 12)}")
print(f"  sum over the 2 vertex pairs (top v3+v4, bottom v1+v2) = {rule2} = {sp.N(rule2, 12)}")
# two-vertex analogue of (3.16): move only the top pair; compare min eig of the sum with the sum of the minima
two_true = sp.simplify((mineig(top) - xi*d1top)/lam1)
two_rule = sp.simplify((mineig(parts[2][1]) - xi*parts[2][0] + mineig(parts[3][1]) - xi*parts[3][0])/lam1)
print(f"\nmoving only v3 and v4 (top pair): true slope = {two_true} = {sp.N(two_true, 12)}; sum of the separate slopes = {two_rule} = {sp.N(two_rule, 12)}")
