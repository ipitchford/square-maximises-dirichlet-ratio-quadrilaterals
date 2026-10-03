include("../geom.jl"); include("../rules.jl"); include("../fsolve.jl"); include("../tm.jl"); include("../assemble.jl"); include("../criterion.jl")
using Printf
println(config_rat(joinpath(@__DIR__, "..", "config.toml"), "cube"), "  ", config_rat(joinpath(@__DIR__, "..", "config.toml"), "eta_neartri"))
x = Arb((Arblib.lbound(Arb(0.3)), Arblib.ubound(Arb(0.5)))); println("Arb from Arf tuple: ", x)
function kappa_star(cpool, base, ws, mode, beta, L0, z)
    # largest certified gap: bisection on kappa in [1, 2.5]
    r = crit_eval!(cpool, base.M, base.K, base.A2, ws; mode, beta, L0, kappa = 5//2, z)
    r.ok || return (r, NaN)
    lo, hi = 1.0, 2.5
    for it in 1:14
        mid = (lo + hi)/2
        rr = crit_eval!(cpool, base.M, base.K, base.A2, ws; mode, beta, L0, kappa = rationalize(mid; tol=1e-9), z)
        rr.ok ? (hi = mid) : (lo = mid)
    end
    return (r, hi)
end
awF = newasm(Float64); cpF = critpool(Float64); wsF = newwork(Float64)
for (p0, n) in (([0.0, 0.06, 0.0, 0.0], 8), ([0.0, 0.1, 0.05, 0.0], 8), ([0.3, 0.3, 0.3, 0.1], 8), ([0.2, 0.5, 0.1, -0.3], 10), ([0.1,0.6,0.7,0.5], 10))
    sol = fsolve(p0, n)
    sb = scalar_bounds(p0, fill(1e-12, 4))
    beta = max(sb.beta_rect_ks, max(sb.L_rect, sb.L_hersch) + sb.gapAC)
    L0 = max(sb.L_rect, sb.L_hersch)
    @printf("p0=%s n=%d  rho=%.4f th2=%.4f beta_pt=%.3f (KS/rect %.3f, L0 %.3f + gap %.3f) L_LG=%.4f  xi_pt(LG)=%.4f xi_pt(mono)=%.4f\n", string(p0), n, sol.rho0, sol.theta[2], beta, sb.beta_rect_ks, L0, sb.gapAC,
        lg_bound(sol.A0, sol.A1, sol.A2, beta), sol.theta[2]/lg_bound(sol.A0, sol.A1, sol.A2, beta), sol.theta[2]/L0)
    for r in (0.0, 0.002, 0.005, 0.01, 0.02, 0.04)
        rad = fill(r, 4); rad[1] = p0[1] == 0 ? r : r   # same radius in all coords
        base = assemble_box!(awF, p0, rad, sol)
        sbb = scalar_bounds(p0, rad .+ 1e-15)
        L0b = max(sbb.L_rect, sbb.L_hersch); betab = max(sbb.beta_rect_ks, L0b + sbb.gapAC)
        z = (1.0, 0.0)
        res, ks = kappa_star(cpF, base, wsF, :LG, betab, L0b, z)
        resm, ksm = kappa_star(cpF, base, wsF, :MONO, betab, L0b, z)
        @printf("   r=%.3f beta=%.3f: LG ok=%s kappa*=%.4f which=%s | MONO ok=%s kappa*=%.4f | maxrem K22=%.3g M11=%.3g A2=%.3g\n", r, betab, res.ok, ks, res.which, resm.ok, ksm, remval(base.K[4]), remval(base.M[1]), remval(base.A2))
    end
end
