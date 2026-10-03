include("../geom.jl"); include("../rules.jl"); include("../fsolve.jl"); include("../tm.jl"); include("../assemble.jl"); include("../criterion.jl")
using Printf
function kappa_star(cpool, base, ws, mode, beta, L0, z)
    r = crit_eval!(cpool, base.M, base.K, base.A2, ws; mode, beta, L0, kappa = 5//2, z)
    r.ok || return (r, NaN)
    lo, hi = 1.0, 2.5
    for it in 1:12
        mid = (lo + hi)/2
        rr = crit_eval!(cpool, base.M, base.K, base.A2, ws; mode, beta, L0, kappa = rationalize(mid; tol=1e-9), z)
        rr.ok ? (hi = mid) : (lo = mid)
    end
    return (r, hi)
end
for K in (3, 4)
    set_order!(K)
    awF = newasm(Float64); cpF = critpool(Float64); wsF = newwork(Float64)
    println("=== ORDER $K")
    for (p0, n) in (([0.3, 0.3, 0.3, 0.1], 8), ([0.1,0.6,0.7,0.5], 10), ([0.2, 0.5, 0.1, -0.3], 12), ([0.0, 0.1, 0.05, 0.0], 8), ([0.0, 0.2, 0.1, 0.0], 8), ([0.5,0.3,0.2,0.3],10))
        sol = fsolve(p0, n)
        line = @sprintf("p0=%-22s n=%d :", string(p0), n)
        for r in (0.0, 0.005, 0.01, 0.02, 0.03, 0.05)
            rad = fill(max(r, 1e-14), 4)
            base = assemble_box!(awF, p0, rad, sol)
            sbb = scalar_bounds(p0, rad)
            L0b = max(sbb.L_rect, sbb.L_hersch); betab = max(sbb.beta_rect_ks, L0b + sbb.gapAC)
            res, ks = kappa_star(cpF, base, wsF, :LG, betab, L0b, (1.0, 0.0))
            line *= @sprintf("  r=%.3f:%s", r, isnan(ks) ? "  FAIL" : @sprintf("%.4f", ks))
        end
        println(line)
    end
end
