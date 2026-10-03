include("../cover_t.jl")
using Printf
W = TWorker()
p0 = [0.5965, 0.15, 0.15, -0.15]
for r in (1e-9, 0.003, 0.006, 0.0125)
    rad = fill(r, 4); rad[1] = r/2
    sh = choose_shift(p0); s = tsolve1(p0, 8; shift = sh); R = transport1_refs(Float64, s); transport1_base!(W.tpF, R, p0, rad, W.wsF, sh)
    base = (; M = W.tpF.M, K = W.tpF.K, A2 = W.tpF.A2)
    bt = bound_tms(Float64, p0, rad, W.wsF)
    sb = scalar_bounds(p0, rad); L0c = max(sb.L_rect, sb.L_hersch)
    iKS = findfirst(x -> x[1] == "KS", bt.betas); BKS = bt.betas[iKS][2]
    r1 = crit_eval!(W.cpF, base.M, base.K, base.A2, W.wsF; mode = :LG, beta = BKS, L0 = 0.0)
    @printf("r=%.4f theta1=%.3f theta2=%.3f KS(p0)=%.3f [lo %.3f] gap(p0)=%.3f L0const=%.3f | LG(KS) ok=%s g1ok=%s ups.tr=%.3g\n", r, s.theta[1], s.theta[2], BKS.c[1], tm_lower(BKS), bt.gap.c[1], L0c, r1.ok, r1.g1ok, r1.ups.tr)
    if r1.g1ok
        L1lo = r1.g2lo / tm_upper(W.cpF.g1)
        b2 = boot_tm(W.cpF, base, BKS, bt.gap, W.wsF)
        if b2 !== nothing
            @printf("     L1 lower(const ratio)=%.3f  beta2(p0)=%.3f [lo %.3f, rem %.3g]\n", L1lo, b2.c[1], tm_lower(b2), remval(b2))
            r2 = crit_eval!(W.cpF, base.M, base.K, base.A2, W.wsF; mode = :LG, beta = b2, L0 = 0.0)
            k2 = r2.ok ? kappa_star_t!(W.cpF, base, W.wsF, :LG, b2, 0.0; steps = 8) : NaN
            @printf("     LG(beta2) ok=%s kappa*=%.4f g1lo=%.3g g2lo=%.3g a11up=%.3g tr=%.3g det=%.3g\n", r2.ok, k2, r2.g1lo, r2.g2lo, r2.a11up, r2.ups.tr, r2.ups.det)
            b3 = boot_tm(W.cpF, base, b2, bt.gap, W.wsF)
            b3 === nothing || @printf("     beta3(p0)=%.3f [lo %.3f]\n", b3.c[1], tm_lower(b3))
        else
            println("     boot_tm returned nothing")
        end
    end
end
