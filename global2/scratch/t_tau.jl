include("../geom.jl"); include("../rules.jl"); include("../fsolve.jl")
using LinearAlgebra
# r = Δu1 + ρ0 u1 at an arbitrary point, from the RR solution s
function u1_at(s, V, x, y)
    sc = s.sc; n = s.n; K = n
    Pu = zeros(K+1); dPu = zeros(K+1); d2Pu = zeros(K+1); Pv = zeros(K+1); dPv = zeros(K+1); d2Pv = zeros(K+1)
    legendre012!(Pu, dPu, d2Pu, (x-sc.xc)/sc.sx, K); legendre012!(Pv, dPv, d2Pv, (y-sc.yc)/sc.sy, K)
    b, bx, by, lapb = bubble_f(V, x, y)
    c1 = s.C[:,1]; u=0.0; lap=0.0; ux=0.0; uy=0.0
    for (k,(i,l)) in enumerate(basis_index(n))
        P = Pu[i+1]*Pv[l+1]; Px = dPu[i+1]*Pv[l+1]/sc.sx; Py = Pu[i+1]*dPv[l+1]/sc.sy
        LP = d2Pu[i+1]*Pv[l+1]/sc.sx^2 + Pu[i+1]*d2Pv[l+1]/sc.sy^2
        u += c1[k]*b*P; lap += c1[k]*(lapb*P + 2*(bx*Px+by*Py) + b*LP); ux += c1[k]*(bx*P+b*Px); uy += c1[k]*(by*P+b*Py)
    end
    return u, lap, ux, uy
end
# generic: tau = phi*e, phi(x) = -(1/rho0) ∫_{t0(x)}^{0} r(x + t e) dt  (t0 <= 0 is where the line from x hits the start set)
function lg_with_tau(s, V, n, kpsi, e, startfun; beta=70.0)
    m = npoints(n); R = rule_f(m); nd = nodes_f(V, R); Q = length(nd.x)
    g, gw = rule_f(24).X, rule_f(24).W
    sc = s.sc; rho0 = s.rho0
    gx = zeros(Q); gy = zeros(Q)
    for q in 1:Q
        x = nd.x[q]; y = nd.y[q]
        u, lap, ux, uy = u1_at(s, V, x, y)
        t0 = startfun(x, y)   # signed parameter of the start point along e from x
        acc = 0.0
        for k in eachindex(g)
            tt = t0*(1-g[k]); uu, ll, _, _ = u1_at(s, V, x+tt*e[1], y+tt*e[2]); acc += gw[k]*(ll + rho0*uu)
        end
        phi = -(1/rho0) * (0 - t0) * acc
        gx[q] = ux/rho0 + phi*e[1]; gy[q] = uy/rho0 + phi*e[2]
    end
    pidx = psi_index(kpsi); Np = length(pidx); K = kpsi
    Cx = zeros(Np,Q); Cy = zeros(Np,Q)
    Pu = zeros(K+1); dPu = zeros(K+1); d2Pu = zeros(K+1); Pv = zeros(K+1); dPv = zeros(K+1); d2Pv = zeros(K+1)
    for q in 1:Q
        legendre012!(Pu, dPu, d2Pu, (nd.x[q]-sc.xc)/sc.sx, K); legendre012!(Pv, dPv, d2Pv, (nd.y[q]-sc.yc)/sc.sy, K)
        for (k,(i,l)) in enumerate(pidx); Cx[k,q] = Pu[i+1]*dPv[l+1]/sc.sy; Cy[k,q] = -dPu[i+1]*Pv[l+1]/sc.sx; end
    end
    W = nd.w; sw = sqrt.(abs.(W))
    A = vcat((Cx .* sw')', (Cy .* sw')'); rhs = -vcat(sw .* gx, sw .* gy)
    d = qr(A, ColumnNorm()) \ rhs
    wx = gx .+ Cx'*d; wy = gy .+ Cy'*d
    A2 = dot(W, wx.^2 .+ wy.^2)
    return lg_bound(s.A0, s.A1, A2, beta)
end
cr(a, b) = a[1]*b[2]-a[2]*b[1]
mm = [0.25, 0.4330127]; nrm = [-0.8660254, 0.5]
for dd in (0.03, 0.1, 0.2), (n,k) in ((10,10),(12,14),(10,14))
    V = Float64.([0 0; 1 0; 0.5 0.8660254; (mm+dd*nrm)'])
    s = fsolve_V(V, n; kpsi=k)
    # (A) current piecewise
    LA = lg_bound(s.A0, s.A1, s.A2, 70.0)
    v1 = V[1,:]; v2 = V[2,:]; v3 = V[3,:]; v4 = V[4,:]
    e = v3 - v1
    # (B) start on line v2v4
    sB(x,y) = (dv = v4 - v2; -cr(dv, [x,y]-v2)/cr(dv, e))
    # (C) start on line through v1 perpendicular to e
    sC(x,y) = -dot([x,y]-v1, e)/dot(e,e)
    # (D) x-direction from x_ref = min x
    xref = minimum(V[:,1])
    sD(x,y) = (xref - x)
    LB = lg_with_tau(s, V, n, k, e, sB)
    LC = lg_with_tau(s, V, n, k, e, sC)
    LD = lg_with_tau(s, V, n, k, [1.0,0.0], sD)
    println("d=$dd n=$n k=$k  A(piecewise)=", round(LA,digits=4), "  B(v2v4)=", round(LB,digits=4), "  C(perp@v1)=", round(LC,digits=4), "  D(x-dir)=", round(LD,digits=4))
end
