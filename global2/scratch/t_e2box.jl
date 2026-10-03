include("../geom.jl"); include("../rules.jl"); include("../fsolve.jl"); include("../tm.jl"); include("../assemble.jl"); include("../criterion.jl")
using Printf, Random, LinearAlgebra
set_order!(3)
function point_kappa(Mv, Kv, A2v, beta, L0)
    M = [Mv[1] Mv[2] Mv[3]; Mv[2] Mv[4] Mv[5]; Mv[3] Mv[5] Mv[6]]
    K = [Kv[1] Kv[2] Kv[3]; Kv[2] Kv[4] Kv[5]; Kv[3] Kv[5] Kv[6]]
    th = eigvals(Symmetric(K), Symmetric(M))
    L = beta*M[1,1] > K[1,1] ? (beta*M[1,1] - K[1,1])/(beta*A2v - M[1,1]) : 0.0
    return th[2]/max(L, L0)
end
function kappa_star(cpool, base, ws, mode, beta, L0)
    r = crit_eval!(cpool, base.M, base.K, base.A2, ws; mode, beta, L0, kappa = 5//2)
    r.ok || return NaN
    lo, hi = 1.0, 2.5
    for it in 1:12
        mid = (lo + hi)/2
        rr = crit_eval!(cpool, base.M, base.K, base.A2, ws; mode, beta, L0, kappa = rationalize(mid; tol=1e-9))
        rr.ok ? (hi = mid) : (lo = mid)
    end
    hi
end
awF = newasm(Float64); aw2 = newasm(Float64); cpF = critpool(Float64); wsF = newwork(Float64)
rng = MersenneTwister(5)
degn(f) = [sum(abs, f.c[TAB[].degrange[d+1]]) for d in 0:TAB[].order]
for (p0, n) in (([0.171, 0.374, 0.092, -0.438], 8), ([0.262, 0.31, 0.297, -0.291], 8), ([0.211, 0.446, 0.262, -0.162], 8))
    sol = fsolve_d(p0, n)
    for r in (0.0025, 0.005, 0.01, 0.02)
        rad = fill(r, 4)
        sbb = scalar_bounds(p0, rad); L0b = max(sbb.L_rect, sbb.L_hersch); betab = max(sbb.beta_rect_ks, L0b + sbb.gapAC)
        base = assemble_box!(awF, p0, rad, sol)
        kc = kappa_star(cpF, base, wsF, :LG, betab, L0b)
        kt = 0.0
        for s in 1:40
            t = s <= 16 ? [((s-1)>>k)&1 == 1 ? 1.0 : -1.0 for k in 0:3] : 2 .* rand(rng, 4) .- 1
            b = assemble_box!(aw2, p0 .+ rad .* t, zeros(4), sol)
            kt = max(kt, point_kappa([x.c[1] for x in b.M], [x.c[1] for x in b.K], b.A2.c[1], betab, L0b))
        end
        @printf("p0=%s r=%.4f beta=%.2f L0=%.2f: true(1st order) kappa=%.4f certified=%s | rem/val K11 %.2e M11 %.2e A2 %.2e; deg1/val K11 %.3f\n", string(p0), r, betab, L0b, kt,
            isnan(kc) ? "FAIL" : @sprintf("%.4f", kc), remval(base.K[1])/base.K[1].c[1], remval(base.M[1])/base.M[1].c[1], remval(base.A2)/base.A2.c[1], degn(base.K[1])[2]/base.K[1].c[1])
    end
end
