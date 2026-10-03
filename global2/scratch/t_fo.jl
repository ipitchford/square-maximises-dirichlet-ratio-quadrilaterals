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
for (p0, n) in (([0.5,0.3,0.2,0.3], 10), ([0.3, 0.3, 0.3, 0.1], 8), ([0.0, 0.1, 0.05, 0.0], 8), ([0.2, 0.5, 0.1, -0.3], 12))
    t0 = time(); sol = fsolve_d(p0, n); tsol = time() - t0
    sol0 = fsolve(p0, n)
    for r in (0.005, 0.01, 0.02, 0.04)
        rad = fill(r, 4)
        sbb = scalar_bounds(p0, rad); L0b = max(sbb.L_rect, sbb.L_hersch); betab = max(sbb.beta_rect_ks, L0b + sbb.gapAC)
        t1 = time(); base = assemble_box!(awF, p0, rad, sol); tas = time() - t1
        kc = kappa_star(cpF, base, wsF, :LG, betab, L0b)
        kt = 0.0; kt0 = 0.0
        for s in 1:60
            t = s <= 16 ? [((s-1)>>k)&1 == 1 ? 1.0 : -1.0 for k in 0:3] : 2 .* rand(rng, 4) .- 1
            b = assemble_box!(aw2, p0 .+ rad .* t, zeros(4), sol)
            kt = max(kt, point_kappa([x.c[1] for x in b.M], [x.c[1] for x in b.K], b.A2.c[1], betab, L0b))
            b0 = assemble_box!(aw2, p0 .+ rad .* t, zeros(4), sol0)
            kt0 = max(kt0, point_kappa([x.c[1] for x in b0.M], [x.c[1] for x in b0.K], b0.A2.c[1], betab, L0b))
        end
        @printf("p0=%-22s n=%d r=%.3f: true kappa fixed=%.4f  first-order=%.4f  certified(1st-order)=%s  [solve %.2fs, asm %.3fs]\n", string(p0), n, r, kt0, kt, isnan(kc) ? "FAIL" : @sprintf("%.4f", kc), tsol, tas)
    end
end
