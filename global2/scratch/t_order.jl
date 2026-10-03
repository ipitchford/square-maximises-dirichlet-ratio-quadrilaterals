include("../geom.jl"); include("../rules.jl"); include("../fsolve.jl"); include("../tm.jl"); include("../assemble.jl"); include("../criterion.jl")
using Printf, Random
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
for K in (2, 3, 4)
    set_order!(K)
    awF = newasm(Float64); cpF = critpool(Float64); wsF = newwork(Float64)
    println("=== ORDER $K, nc=$(TAB[].nc)")
    for (p0, n) in (([0.0, 0.06, 0.0, 0.0], 8), ([0.3, 0.3, 0.3, 0.1], 8), ([0.2, 0.5, 0.1, -0.3], 10))
        sol = fsolve(p0, n)
        for r in (0.005, 0.01, 0.02, 0.04)
            rad = fill(r, 4)
            t = @elapsed base = assemble_box!(awF, p0, rad, sol)
            # correctness check: TM vs exact at random points
            rng = MersenneTwister(3); worst = 0.0
            if r == 0.02
                aw2 = newasm(Float64)
                for s in 1:10
                    tt = 2 .* rand(rng, 4) .- 1; p = p0 .+ rad .* tt
                    rp = assemble_box!(aw2, p, zeros(4), sol)
                    worst = max(worst, abs(tm_eval(base.K[4], tt) - rp.K[4].c[1]) / remval(base.K[4]))
                end
            end
            sbb = scalar_bounds(p0, rad)
            L0b = max(sbb.L_rect, sbb.L_hersch); betab = max(sbb.beta_rect_ks, L0b + sbb.gapAC)
            res, ks = kappa_star(cpF, base, wsF, :LG, betab, L0b, (1.0, 0.0))
            @printf("  p0=%s n=%d r=%.3f: kappa*=%.4f ok=%s | rem K22=%.3g M11=%.3g A2=%.3g | t=%.3fs %s\n", string(p0), n, r, ks, res.ok, remval(base.K[4]), remval(base.M[1]), remval(base.A2), t,
                r == 0.02 ? @sprintf("max|TM-exact|/rem=%.3f", worst) : "")
        end
    end
end
