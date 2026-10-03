include("../geom.jl"); include("../rules.jl"); include("../fsolve.jl")
using Random, Printf
import TOML
cfg = TOML.parsefile(joinpath(@__DIR__, "..", "config.toml"))
cube = cfg["cube"]; eta = cfg["eta_neartri"]; s0 = cfg["thin_s0"]
# inradius of convex quadrilateral by LP vertices (equidistant to 3 edge lines)
function inradius(V)
    nrm = []; off = []
    for k in 1:4
        j = mod1(k+1,4); ex = V[j,1]-V[k,1]; ey = V[j,2]-V[k,2]; L = hypot(ex,ey)
        n = [-ey, ex] ./ L   # inward normal (ccw)
        push!(nrm, n); push!(off, n[1]*V[k,1] + n[2]*V[k,2])   # n·x - off >= 0 inside
    end
    best = 0.0
    for drop in 1:4
        idx = [k for k in 1:4 if k != drop]
        A = [nrm[idx[i]][j] for i in 1:3, j in 1:2]; A = hcat(A, -ones(3)); b = [off[i] for i in idx]
        abs(det(A)) < 1e-14 && continue
        sol = A \ b; x = sol[1:2]; r = sol[3]
        r > 0 || continue
        dd = nrm[drop][1]*x[1] + nrm[drop][2]*x[2] - off[drop]
        dd >= r - 1e-12 && (best = max(best, r))
    end
    return best
end
rng = MersenneTwister(7)
lo = [0.0,0,0,-0.905]; hi = [0.905,1.4,1.4,0.905]
res = []
while length(res) < 1500
    p = lo .+ rand(rng, 4) .* (hi .- lo)
    classify_point(p; cube, eta, s0) == :R || continue
    s = fsolve(p, 10); V = verts(p)
    beta = beta_point_f(p); Lm = lmono_point_f(p)
    r = inradius(V); LH = π^2/(4r^2)
    w, D = width_diam(V)
    L1 = max(Lm, LH)
    # AC bootstrap: beta_AC = L + 3pi^2/D^2 (iterate with LG)
    b = max(beta, L1 + 3π^2/D^2); LG = -Inf
    for it in 1:20
        b * s.A1 > s.A0 || break
        LG = lg_bound(s.A0, s.A1, s.A2, b)
        bn = max(beta, max(LG, L1) + 3π^2/D^2)
        bn <= b + 1e-12 && break
        b = bn
    end
    xb = s.theta[2] / max(L1, LG)
    push!(res, (p=p, x=xb, xm=s.theta[2]/Lm, xh=s.theta[2]/LH, xlg=LG > 0 ? s.theta[2]/LG : Inf, wd=w/D, ang=maximum(angles_deg(V)), b_rho=b/s.rho0))
end
x = sort([r.x for r in res]); println("best-bound quantiles 50/90/99/100%: ", [round(x[clamp(round(Int,q*length(x)),1,length(x))],digits=4) for q in (0.5,0.9,0.99,1.0)])
println("count >= 2.4: ", count(r->r.x>=2.4, res), "   >= 2.3: ", count(r->r.x>=2.3, res))
for r in sort(res, by=r->-r.x)[1:12]
    @printf("  p=%s best=%.3f mono=%.3f hersch=%.3f lg(AC)=%.3f w/D=%.3f ang=%.1f beta/rho=%.3f\n", string(round.(r.p,digits=3)), r.x, r.xm, r.xh, r.xlg, r.wd, r.ang, r.b_rho)
end
