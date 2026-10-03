include("../geom.jl"); include("../rules.jl"); include("../fsolve.jl")
using Random, Printf
import TOML
cfg = TOML.parsefile(joinpath(@__DIR__, "..", "config.toml"))
cube = cfg["cube"]; eta = cfg["eta_neartri"]; s0 = cfg["thin_s0"]
rng = MersenneTwister(7)
lo = [0.0,0,0,-0.905]; hi = [0.905,1.4,1.4,0.905]
res = []
t0 = time()
while length(res) < 1500
    p = lo .+ rand(rng, 4) .* (hi .- lo)
    classify_point(p; cube, eta, s0) == :R || continue
    s = fsolve(p, 10)
    beta = beta_point_f(p); Lm = lmono_point_f(p)
    LG = beta * s.A1 > s.A0 ? lg_bound(s.A0, s.A1, s.A2, beta) : -Inf
    xb_lg = LG > 0 ? s.theta[2]/LG : Inf; xb_m = s.theta[2]/Lm
    w, D = width_diam(verts(p))
    push!(res, (p=p, xlg=xb_lg, xm=xb_m, x=min(xb_lg, xb_m), beta_over_rho=beta/s.rho0, wd=w/D, maxang=maximum(angles_deg(verts(p))), th2=s.theta[2], rho=s.rho0))
end
println("time ", time()-t0)
x = [r.x for r in res]
println("best xi-bound quantiles: ", [round(sort(x)[clamp(round(Int, q*length(x)),1,length(x))], digits=4) for q in (0.5, 0.9, 0.99, 1.0)]) 
nlg = count(r -> isfinite(r.xlg), res); println("LG available at ", nlg, " of ", length(res))
bad = filter(r -> r.x >= 2.4, res)
println("points with best bound >= 2.4: ", length(bad))
for r in sort(bad, by = r -> -r.x)[1:min(end,25)]
    @printf("  p=%s xi_LG=%.3f xi_mono=%.3f beta/rho=%.3f w/D=%.3f maxang=%.1f\n", string(round.(r.p,digits=3)), r.xlg, r.xm, r.beta_over_rho, r.wd, r.maxang)
end
# classify which bound wins
println("LG wins: ", count(r -> r.xlg < r.xm, res), "  mono wins: ", count(r -> r.xm <= r.xlg, res))
