include("../cover_t.jl"); include("../transport1.jl")
using Printf, Random, LinearAlgebra
W = TWorker(); rng = MersenneTwister(9)
for K in (3, 4)
    set_order!(K); tpF = tpool(Float64); cpF = critpool(Float64); wsF = newwork(Float64)
    for p0 in ([0.3,0.2,0.4,0.1], [0.171, 0.374, 0.092, -0.438])
        sh = choose_shift(p0); s = tsolve1(p0, 8; shift = sh)
        R = transport1_refs(Float64, s)
        # point check vs pencil_at
        Kp, Mp, TT, TC, CC = s.pencil_at(p0)
        Kd = s.C' * Kp * s.C
        b0 = transport1_base!(tpF, R, p0, fill(1e-14, 4), wsF, sh)
        D0 = tpF.D.c[1]
        @printf("ORDER %d p0=%s: K11 refs %.10f vs full %.10f | A2 %.10f vs %.10f\n", K, string(p0), tpF.K[1].c[1]/D0, Kd[1,1],
            tpF.A2.c[1]/D0, s.C[:,1]'*TT*s.C[:,1] + 2*s.C[:,1]'*TC*s.dpsi + s.dpsi'*CC*s.dpsi)
        for r in (0.02, 0.04, 0.06)
            rad = fill(r, 4)
            transport1_base!(tpF, R, p0, rad, wsF, sh)
            # TM vs direct at random points
            worst = 0.0
            for smp in 1:10
                t = 2 .* rand(rng, 4) .- 1; p = p0 .+ rad .* t; δ = rad .* t
                Kp, Mp, TT, TC, CC = s.pencil_at(p)
                C1 = s.C .+ sum(δ[i] .* s.dC[:, :, i] for i in 1:4); psi1 = s.dpsi .+ sum(δ[i] .* s.ddpsi[:, i] for i in 1:4)
                Kx = C1' * Kp * C1; d = tm_eval(tpF.D, t)
                worst = max(worst, abs(tm_eval(tpF.K[4], t) - d*Kx[2,2]) / remval(tpF.K[4]))
            end
            bt = bound_tms(Float64, p0, rad, wsF)
            best = Inf
            for (nm, B) in bt.betas
                rr = crit_eval!(cpF, tpF.M, tpF.K, tpF.A2, wsF; mode = :LG, beta = B, L0 = 0.0); rr.ok || continue
                k = kappa_star_t!(cpF, (; M = tpF.M, K = tpF.K, A2 = tpF.A2), wsF, :LG, B, 0.0; steps = 10); best = min(best, k)
            end
            for (nm, L) in bt.L0s
                rr = crit_eval!(cpF, tpF.M, tpF.K, tpF.A2, wsF; mode = :MONO, beta = 0.0, L0 = L); rr.ok || continue
                k = kappa_star_t!(cpF, (; M = tpF.M, K = tpF.K, A2 = tpF.A2), wsF, :MONO, 0.0, L; steps = 10); best = min(best, k)
            end
            @printf("    r=%.2f: certified kappa* %s   (TM err/rem %.2f)\n", r, isinf(best) ? "FAIL" : @sprintf("%.4f", best), worst)
        end
    end
end
