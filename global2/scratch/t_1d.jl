include("../cover_t.jl"); include("../bound1d.jl")
using Printf
for p0 in ([0.5965, 0.15, 0.15, -0.15], [0.321, 0.606, 1.013, 0.668], [0.5, 0.3, 0.9, -0.2], [0.3, 0.5, 0.8, 0.5], [0.3,0.2,0.4,0.1], [0.74, 0.025, 0.15, 0.222])
    s = fsolve(p0, 10); w, D = width_diam(verts(p0))
    line = @sprintf("p0=%-26s w/D=%.3f theta1=%7.3f theta2=%7.3f |", string(p0), w/D, s.theta[1], s.theta[2])
    for r in (1e-9, 0.005, 0.01, 0.02)
        t = @elapsed L1 = lambda1_1d(p0, fill(r, 4))
        sb = scalar_bounds(p0, fill(r, 4))
        line *= @sprintf(" r=%.3f: 1D %.2f (%.0f%%) old %.2f  [%.2fs]", r, L1, 100L1/s.theta[1], max(sb.L_rect, sb.L_hersch), t)
    end
    println(line)
end
