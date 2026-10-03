"""Exact check of (L2): for each parameter i,  d_i LL_11(0) - lambda_1^2 d_i M_11(0) - 2 lambda_1 (d_i K_11(0) - lambda_1 d_i M_11(0)) = 0,
where LL_11 = int (Lg_1)^2 J, Lg = (1/J) d_al ( (B^{al be}/J) d_be g ).  (Reference square [0,1]^2.)"""
import sympy as sp
u,v,a,b,c,d = sp.symbols('u v a b c d', real=True); P=(a,b,c,d); half=sp.Rational(1,2)
d1 = (half*(a+d-b), half*(-a+d+c)); d2 = (half*(-a+d+b), half*(-a-d-c)); d3 = (half*(-a-d-b), half*(a-d+c)); d4 = (half*(a-d+b), half*(a+d-c))
vx = [d1[0], 1+d2[0], 1+d3[0], d4[0]]; vy = [d1[1], d2[1], 1+d3[1], 1+d4[1]]
N = [(1-u)*(1-v), u*(1-v), u*v, (1-u)*v]
X = sum(N[i]*vx[i] for i in range(4)); Y = sum(N[i]*vy[i] for i in range(4))
G = sp.Matrix([[sp.diff(X,u), sp.diff(X,v)],[sp.diff(Y,u), sp.diff(Y,v)]]); J = sp.expand(G.det()); adj = G.adjugate(); B = adj*adj.T
g1 = 2*sp.sin(sp.pi*u)*sp.sin(sp.pi*v); lam1 = 2*sp.pi**2
T = B/J
grad = sp.Matrix([sp.diff(g1,u), sp.diff(g1,v)])
flux = T*grad
Lg = (sp.diff(flux[0],u) + sp.diff(flux[1],v))/J
def fo(expr, i): return sp.diff(expr, P[i]).subs({a:0,b:0,c:0,d:0})
def integ(e): return sp.integrate(sp.integrate(sp.expand(e), (u,0,1)), (v,0,1))
for i in range(4):
    dLL = integ(sp.simplify(fo(Lg*Lg*J, i)))
    dM = integ(fo(g1*g1*J, i))
    dK = integ(sp.simplify(fo((grad.T*T*grad)[0,0], i)))
    val = sp.simplify(dLL - lam1**2*dM - 2*lam1*(dK - lam1*dM))
    print(f"param {'abcd'[i]}: dLL={sp.simplify(dLL)} dM={dM} dK={sp.simplify(dK)}  ->  (L2) residual identity value = {val}")
    assert val == 0, 'L2 fails'

# L0 for LL (LL(0) check): LL_kl(0) = int (Lg_k)(Lg_l) J at p=0 = lambda_k lambda_l delta_kl for the base modes
zero = {a:0, b:0, c:0, d:0}
modes = {1:(1,1), 2:(2,1), 3:(1,2)}
def Lg_of(m, n):
    g = 2*sp.sin(m*sp.pi*u)*sp.sin(n*sp.pi*v)
    gr = sp.Matrix([sp.diff(g,u), sp.diff(g,v)]); fl = (B/J)*gr
    return ((sp.diff(fl[0],u) + sp.diff(fl[1],v))/J).subs(zero)
Ls = {k: sp.simplify(Lg_of(*mn)) for k, mn in modes.items()}
for k in modes:
    for l in modes:
        if l < k: continue
        val = sp.simplify(sp.integrate(sp.integrate(sp.expand(Ls[k]*Ls[l]), (u,0,1)), (v,0,1)))
        lk = sp.pi**2*(modes[k][0]**2+modes[k][1]**2); ll = sp.pi**2*(modes[l][0]**2+modes[l][1]**2)
        assert sp.simplify(val - (lk*ll if k == l else 0)) == 0, ('LL(0) fails', k, l, val)
print('LL(0) check: L0 for LL VERIFIED (exact)')
print('LEMMA L2: VERIFIED (exact)')
