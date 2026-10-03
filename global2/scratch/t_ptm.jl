include("../cover_t.jl"); include("../bounds_tm.jl")
using Printf
W = TWorker()
function kstar(cp, tp, ws, mode, beta, L0; steps = 10)
    r = crit_eval!(cp, tp.M, tp.K, tp.A2, ws; mode, beta, L0); r.ok || return NaN
    lo, hi = 1.0, 2.5
    for _ in 1:steps
        mid = (lo+hi)/2
        rr = crit_eval!(cp, tp.M, tp.K, tp.A2, ws; mode, beta, L0, kappa = rationalize(mid; tol=1e-9)); rr.ok ? (hi = mid) : (lo = mid)
    end
    hi
end
for p0 in ([0.321, 0.606, 1.013, 0.668], [0.3, 0.5, 0.8, 0.5], [0.5, 0.3, 0.9, -0.2], [0.3,0.2,0.4,0.1], [0.171, 0.374, 0.092, -0.438], [0.343, 0.671, 1.061, 0.712])
    line = @sprintf("p0=%-26s", string(p0))
    for r in (0.005, 0.01, 0.02, 0.04)
        rad = fill(r, 4)
        sb = scalar_bounds(p0, rad); L0c = max(sb.L_rect, sb.L_hersch); betac = max(sb.beta_rect_ks, addlo(L0c, sb.gapAC))
        bt = bound_tms(Float64, p0, rad, W.wsF)
        best_c = Inf; best_t = Inf
        for sh in (0, 1)
            s = tsolve(p0, 8; shift = sh); R = transport_refs(Float64, s)
            transport_base!(W.tpF, R, p0, rad, W.wsF, sh)
            for mode in (:LG, :MONO)
                (k = kstar(W.cpF, W.tpF, W.wsF, mode, betac, L0c); isnan(k) || (best_c = min(best_c, k)))
            end
            for (nm, L) in bt.L0s; k = kstar(W.cpF, W.tpF, W.wsF, :MONO, 0.0, L); isnan(k) || (best_t = min(best_t, k)); end
            for (nm, B) in bt.betas; k = kstar(W.cpF, W.tpF, W.wsF, :LG, B, 0.0); isnan(k) || (best_t = min(best_t, k)); end
        end
        line *= @sprintf("  r=%.3f const %s | TM %s", r, isinf(best_c) || isnan(best_c) ? "FAIL " : @sprintf("%.3f", best_c), isinf(best_t) || isnan(best_t) ? "FAIL " : @sprintf("%.3f", best_t))
    end
    println(line)
end
