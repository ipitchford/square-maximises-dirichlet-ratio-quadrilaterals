include("../cover_t.jl")
using Printf
for p0 in ([0.321, 0.606, 1.013, 0.668], [0.343, 0.671, 1.061, 0.712], [0.3, 0.5, 0.8, 0.5], [0.5, 0.3, 0.9, -0.2])
    s = fsolve(p0, 10); V = verts(p0); w, D = width_diam(V)
    for r in (1e-12, 0.01, 0.02)
        sb = scalar_bounds(p0, fill(r, 4)); L0 = max(sb.L_rect, sb.L_hersch)
        @printf("p0=%s r=%.2f: theta1=%.2f theta2=%.2f | beta_rect/KS=%.2f L_rect=%.2f L_hersch=%.2f gapAC=%.2f | w/D=%.3f maxang=%.0f  theta2/L0=%.3f\n",
            string(p0), r, s.theta[1], s.theta[2], sb.beta_rect_ks, sb.L_rect, sb.L_hersch, sb.gapAC, w/D, maximum(angles_deg(V)), s.theta[2]/L0)
    end
end
