include("../geom.jl"); include("../rules.jl"); include("../fsolve.jl"); include("../tm.jl"); include("../assemble.jl"); include("../criterion.jl")
using Random, Printf
rng = MersenneTwister(9)
pts = []
lo = [0.0,0,0,-0.905]; hi = [0.905,1.4,1.4,0.905]
while length(pts) < 10
    p = lo .+ rand(rng, 4) .* (hi .- lo)
    classify_point(p; cube = 0.049, eta = 0.035, s0 = 0.1) == :R || continue
    c = cvals(p); r = minimum(c[i]/c[mod1(i+2,4)] for i in 1:4)
    0.035 < r < 0.045 || continue
    push!(pts, p)
end
for p in pts
    sb = scalar_bounds(p, fill(1e-12,4)); L0 = max(sb.L_rect, sb.L_hersch); beta = max(sb.beta_rect_ks, L0 + sb.gapAC)
    c = cvals(p); r = minimum(c[i]/c[mod1(i+2,4)] for i in 1:4)
    line = @sprintf("p=%s ratio=%.4f maxang=%.1f w/D=%.2f:", string(round.(p,digits=3)), r, maximum(angles_deg(verts(p))), (x->x[1]/x[2])(width_diam(verts(p))))
    for n in (8, 12, 16)
        s = fsolve(p, n)
        LG = beta*s.A1 > s.A0 ? lg_bound(s.A0, s.A1, s.A2, beta) : 0.0
        # AC bootstrap
        b2 = max(beta, LG + sb.gapAC); LG2 = b2*s.A1 > s.A0 ? lg_bound(s.A0, s.A1, s.A2, b2) : 0.0
        line *= @sprintf("  n=%d k=%.4f", n, s.theta[2]/max(LG2, L0))
    end
    println(line)
end
