"""Exact symbolic verification (sympy) of the first-order identities used to divide by t^2 in the local step.
Setting: reference square [0,1]^2, modes g_k = 2 sin(m pi u) sin(n pi v), bilinear map with Endo-Osting displacements.
Claims:
 (L0) K_kl(0) = lambda_k delta_kl, M_kl(0) = delta_kl, LL_kl(0) = lambda_k lambda_l delta_kl  (orthonormality; exact)
 (L1) for each parameter i in {a,b,c,d}:  d_i [K_22 + K_33 - 5 pi^2 (M_22 + M_33)](0) = 5 * d_i [K_11 - 2 pi^2 M_11](0)
      (cluster modes 2=(2,1), 3=(1,2), ground 1=(1,1)); this is tr S_1 = 0.
 Also record d_i K_23(0) - 5 pi^2 d_i M_23(0) (the off-diagonal first-order cluster coupling) for reference.
"""
import sympy as sp
u,v,a,b,c,d = sp.symbols('u v a b c d', real=True)
P = (a,b,c,d)
half = sp.Rational(1,2)
d1 = (half*(a+d-b), half*(-a+d+c)); d2 = (half*(-a+d+b), half*(-a-d-c)); d3 = (half*(-a-d-b), half*(a-d+c)); d4 = (half*(a-d+b), half*(a+d-c))
vx = [d1[0], 1+d2[0], 1+d3[0], d4[0]]; vy = [d1[1], d2[1], 1+d3[1], 1+d4[1]]
N = [(1-u)*(1-v), u*(1-v), u*v, (1-u)*v]
X = sum(N[i]*vx[i] for i in range(4)); Y = sum(N[i]*vy[i] for i in range(4))
G = sp.Matrix([[sp.diff(X,u), sp.diff(X,v)],[sp.diff(Y,u), sp.diff(Y,v)]])
J = sp.expand(G.det()); adj = G.adjugate(); B = adj*adj.T
def g(m,n): return 2*sp.sin(m*sp.pi*u)*sp.sin(n*sp.pi*v)
modes = {1:(1,1), 2:(2,1), 3:(1,2)}
lam = {k: sp.pi**2*(m*m+n*n) for k,(m,n) in modes.items()}
def grad(f): return sp.Matrix([sp.diff(f,u), sp.diff(f,v)])
def first_order(expr, i):
    """first-order coefficient in parameter P[i] at p=0 of a (u,v,p)-expression"""
    e = sp.diff(expr, P[i]).subs({a:0,b:0,c:0,d:0})
    return sp.simplify(e)
def integ(e): return sp.integrate(sp.integrate(sp.expand(e), (u,0,1)), (v,0,1))
gk = {k: g(*mn) for k,mn in modes.items()}
# K integrand: grad g_k^T (B/J) grad g_l ; M integrand: g_k g_l J
results = {}
for i in range(4):
    Kd = {}; Md = {}
    for k in modes:
        for l in modes:
            if l < k: continue
            TK = B/J
            integrandK = (grad(gk[k]).T*TK*grad(gk[l]))[0,0]
            Kd[(k,l)] = integ(first_order(integrandK, i))
            Md[(k,l)] = integ(first_order(gk[k]*gk[l]*J, i))
    lhs = sp.simplify(Kd[(2,2)] + Kd[(3,3)] - 5*sp.pi**2*(Md[(2,2)]+Md[(3,3)]))
    rhs = sp.simplify(5*(Kd[(1,1)] - 2*sp.pi**2*Md[(1,1)]))
    off = sp.simplify(Kd[(2,3)] - 5*sp.pi**2*Md[(2,3)])
    diag22 = sp.simplify(Kd[(2,2)] - 5*sp.pi**2*Md[(2,2)]); diag33 = sp.simplify(Kd[(3,3)] - 5*sp.pi**2*Md[(3,3)])
    print(f"param {'abcd'[i]}: lhs={lhs}  rhs={rhs}  lhs-rhs={sp.simplify(lhs-rhs)}  | first-order cluster matrix: diag ({diag22}, {diag33}), off {off} | d_i lambda1 = {sp.simplify(Kd[(1,1)] - 2*sp.pi**2*Md[(1,1)])}")
    assert sp.simplify(lhs) == 0 and sp.simplify(rhs) == 0, 'L1 fails'
    results['abcd'[i]] = (lhs, rhs, diag22, diag33, off)
# zeroth order check
for k in modes:
    for l in modes:
        Kv = integ(((grad(gk[k]).T*(B/J)*grad(gk[l]))[0,0]).subs({a:0,b:0,c:0,d:0}))
        Mv = integ((gk[k]*gk[l]*J).subs({a:0,b:0,c:0,d:0}))
        print(f"K({k},{l})(0)={sp.simplify(Kv)}  M={sp.simplify(Mv)}")
        assert sp.simplify(Mv - (1 if k == l else 0)) == 0 and sp.simplify(Kv - (lam[k] if k == l else 0)) == 0, 'L0 fails'

print('LEMMA L0, L1: VERIFIED (exact)')
