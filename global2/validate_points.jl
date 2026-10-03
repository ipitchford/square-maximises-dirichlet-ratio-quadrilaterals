# Float64 pilot validation: pick >= 24 points of the uncovered region R (near square, near-triangle,
# width/diam ~ 0.1-0.2, random), evaluate the Float64 certificate pipeline (RR theta, LG bound with the
# point beta), write them for the FEM comparison (fem_compare.py), and check D4 invariance.
include("geom.jl"); include("rules.jl"); include("fsolve.jl")
using Random, Printf
import TOML
cfg = TOML.parsefile(joinpath(@__DIR__, "config.toml"))
cube = cfg["cube"]; eta = cfg["eta_neartri"]; s0 = cfg["thin_s0"]
rng = MersenneTwister(20261001)
pts = Vector{Vector{Float64}}(); tags = String[]
for p in ([0,0.05,0,0], [0,0.042,0.042,0], [0.06,0.03,0.02,0.01], [0,0.1,0.1,0], [0.05,0,0,0.05], [0.02,0.07,0.01,-0.03], [0.0,0.0,0.06,0.0])
    push!(pts, Float64.(p)); push!(tags, "nearsq")
end
function ratio_min(p); c = cvals(p); minimum(c[i]/c[mod1(i+2,4)] for i=1:4); end
lo = [0.0,0,0,-0.9]; hi = [0.9,1.3,1.3,0.9]
cnt = Dict("neartri"=>0, "thin"=>0, "random"=>0)
while cnt["neartri"] < 6 || cnt["thin"] < 6 || cnt["random"] < 6
    p = lo .+ rand(rng, 4) .* (hi .- lo)
    classify_point(p; cube, eta, s0) == :R || continue
    rm = ratio_min(p); w, D = width_diam(verts(p))
    if rm < 0.06 && cnt["neartri"] < 6
        push!(pts, p); push!(tags, "neartri"); cnt["neartri"] += 1
    elseif w/D < 0.2 && cnt["thin"] < 6
        push!(pts, p); push!(tags, "thin"); cnt["thin"] += 1
    elseif cnt["random"] < 6 && rm >= 0.06 && w/D >= 0.2
        push!(pts, p); push!(tags, "random"); cnt["random"] += 1
    end
end
open(joinpath(@__DIR__, "validation_points.csv"), "w") do io
    println(io, "tag,a,b,c,d,v1x,v1y,v2x,v2y,v3x,v3y,v4x,v4y,n,theta1,theta2,theta3,A2,beta_pt,LG_pt,Lmono_pt,maxangle,wd,ratio_min")
    for (p, tag) in zip(pts, tags)
        V = verts(p); w, D = width_diam(V)
        for n in (8, 12)
            s = fsolve(p, n)
            beta = beta_point_f(p)
            L = beta * s.A1 > s.A0 ? lg_bound(s.A0, s.A1, s.A2, beta) : NaN
            @printf(io, "%s,%.17g,%.17g,%.17g,%.17g,%s,%d,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.6g,%.6g,%.6g\n", tag, p..., join([@sprintf("%.17g", V[j,k]) for j=1:4 for k=1:2], ","), n,
                s.theta[1], s.theta[2], s.theta[3], s.A2, beta, L, lmono_point_f(p), maximum(angles_deg(V)), w/D, ratio_min(p))
        end
    end
end
println("wrote ", length(pts), " points")
# D4 invariance of the Float64 pipeline
println("D4 invariance (n=10): max relative spread of theta1..3 and A2 over the 8 images")
for p in pts[[1, 3, 9, 15, 21]]
    imgs = [p]
    for _ in 1:3; push!(imgs, rho_rot(imgs[end])); end
    imgs = vcat(imgs, sigma_h.(imgs))
    vals = [(s = fsolve(q, 10); [s.theta[1:3]; s.A2]) for q in imgs]
    M = hcat(vals...)
    spread = [maximum(M[k,:]) - minimum(M[k,:]) for k in 1:4] ./ abs.(M[:,1])
    @printf("  p=%s  spread = %s\n", string(round.(p, digits=3)), string(round.(spread, sigdigits=2)))
end
