include("../geom.jl"); include("../rules.jl"); include("../fsolve.jl"); include("../tm.jl"); include("../assemble.jl"); include("../criterion.jl")
using Printf, Random, LinearAlgebra
import TOML
set_order!(3)
cfg = TOML.parsefile(joinpath(@__DIR__, "..", "config.toml")); cube = cfg["cube"]; eta = cfg["eta_neartri"]; s0 = cfg["thin_s0"]
function point_kappa(Mv, Kv, A2v, beta, L0)
    M = [Mv[1] Mv[2] Mv[3]; Mv[2] Mv[4] Mv[5]; Mv[3] Mv[5] Mv[6]]
    K = [Kv[1] Kv[2] Kv[3]; Kv[2] Kv[4] Kv[5]; Kv[3] Kv[5] Kv[6]]
    th = eigvals(Symmetric(K), Symmetric(M))
    L = beta*M[1,1] > K[1,1] ? (beta*M[1,1] - K[1,1])/(beta*A2v - M[1,1]) : 0.0
    return th[2]/max(L, L0)
end
rng = MersenneTwister(11)
lo = [0.0,0,0,-0.905]; hi = [0.905,1.4,1.4,0.905]
aw2 = newasm(Float64)
radii = (0.0025, 0.005, 0.01, 0.02, 0.04)
res = []
while length(res) < 120
    p0 = lo .+ rand(rng, 4) .* (hi .- lo)
    classify_point(p0; cube, eta, s0) == :R || continue
    best = 0.0
    sol = fsolve(p0, 10)
    rmax = 0.0
    for r in radii
        rad = fill(r, 4)
        sbb = scalar_bounds(p0, rad); L0b = max(sbb.L_rect, sbb.L_hersch); betab = max(sbb.beta_rect_ks, L0b + sbb.gapAC)
        kt = 0.0
        for s in 1:40
            t = s <= 16 ? [((s-1)>>k)&1 == 1 ? 1.0 : -1.0 for k in 0:3] : 2 .* rand(rng, 4) .- 1
            b = assemble_box!(aw2, p0 .+ rad .* t, zeros(4), sol)
            kt = max(kt, point_kappa([x.c[1] for x in b.M], [x.c[1] for x in b.K], b.A2.c[1], betab, L0b))
        end
        kt < 2.45 ? (rmax = r) : break
    end
    V = verts(p0); w, D = width_diam(V)
    push!(res, (p0=p0, rmax=rmax, ang=maximum(angles_deg(V)), wd=w/D))
end
using Statistics
for r in (0.0, radii...)
    println("rmax == $r: ", count(x -> x.rmax == r, res))
end
# volume-weighted box count estimate: each sample represents vol(R)/N; boxes needed ~ (vol/N)/(2 rmax)^4
volR = 0.29
est = sum(x -> x.rmax > 0 ? (volR/length(res))/(2x.rmax)^4 : 0.0, res)
println("box count estimate (true loss only, cube boxes): ", round(est, sigdigits=3), "   points with rmax=0: ", count(x->x.rmax==0, res))
println("max angle of rmax<=0.0025 points: ", [round(x.ang) for x in res if x.rmax <= 0.0025])
println("max angle of rmax>=0.02 points: ", [round(x.ang) for x in res if x.rmax >= 0.02])
