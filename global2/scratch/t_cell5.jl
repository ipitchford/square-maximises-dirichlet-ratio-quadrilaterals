include("../cover_t.jl")
using Printf, LinearAlgebra
lines = filter(l -> occursin("pilot_pass", l), readlines(joinpath(@__DIR__, "..", "runs", "pilot_t1", "cell_5.jsonl")))
W = TWorker()
for ln in lines[1:4]
    lo = parse.(Float64, split(match(r"\"lo\":\[([^\]]+)\]", ln).captures[1], ",")); hi = parse.(Float64, split(match(r"\"hi\":\[([^\]]+)\]", ln).captures[1], ","))
    p0, rad = box_center_radius(lo, hi); V = verts(p0); w, D = width_diam(V)
    sb = scalar_bounds(p0, rad); L0 = max(sb.L_rect, sb.L_hersch); beta0 = max(sb.beta_rect_ks, addlo(L0, sb.gapAC))
    bt = bound_tms(Float64, p0, rad, W.wsF)
    line = @sprintf("p0=%s rad=%s w/D=%.3f maxang=%.0f |", string(round.(p0,digits=3)), string(round.(rad,digits=4)), w/D, maximum(angles_deg(V)))
    for sh in (0, 1)
        s1 = tsolve1(p0, 8; shift = sh)
        s0 = merge(s1, (; dC = zero(s1.dC), ddpsi = zero(s1.ddpsi)))
        for (lab, s) in (("plain", s0), ("1st", s1))
            R = transport1_refs(Float64, s); transport1_base!(W.tpF, R, p0, rad, W.wsF, sh)
            best = Inf
            for (mode, B, L) in vcat([(:LG, beta0, L0), (:MONO, beta0, L0)], [(:LG, b[2], L0) for b in bt.betas], [(:MONO, beta0, l[2]) for l in bt.L0s])
                rr = crit_eval!(W.cpF, W.tpF.M, W.tpF.K, W.tpF.A2, W.wsF; mode, beta = B, L0 = L); rr.ok || continue
                best = min(best, kappa_star_t!(W.cpF, (; M = W.tpF.M, K = W.tpF.K, A2 = W.tpF.A2), W.wsF, mode, B, L; steps = 8))
            end
            line *= @sprintf(" sh%d %s %s (|dC| %.1e)", sh, lab, isinf(best) ? "FAIL " : @sprintf("%.3f", best), maximum(abs, s.dC))
        end
    end
    println(line)
end
