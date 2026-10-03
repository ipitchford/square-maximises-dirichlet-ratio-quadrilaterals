include("../cover_t.jl"); include("../transport1.jl")
using Printf, Random, LinearAlgebra
W = TWorker(); rng = MersenneTwister(9)
function kq(Kp, Mp, TT, TC, CC, C, psi, beta, L0)
    K = C' * Kp * C; M = C' * Mp * C; c1 = C[:, 1]
    A2 = c1' * TT * c1 + 2 * c1' * TC * psi + psi' * CC * psi
    th = eigvals(Symmetric(K), Symmetric(M)); L = beta*M[1,1] > K[1,1] ? (beta*M[1,1]-K[1,1])/(beta*A2-M[1,1]) : 0.0
    th[2]/max(L, L0)
end
for p0 in ([0.3,0.2,0.4,0.1], [0.171, 0.374, 0.092, -0.438], [0.5, 0.3, 0.2, 0.3], [0.321, 0.606, 1.013, 0.668], [0.262, 0.31, 0.297, -0.291])
    sh = choose_shift(p0); t0 = time(); s = tsolve1(p0, 8; shift = sh); ts = time() - t0
    line = @sprintf("p0=%-28s [%.2fs]", string(p0), ts)
    for r in (0.02, 0.04, 0.06)
        rad = fill(r, 4); bt = bound_tms(Float64, p0, rad, W.wsF)
        kfix = 0.0; kfo = 0.0
        for smp in 1:24
            t = smp <= 16 ? [((smp-1)>>k)&1 == 1 ? 1.0 : -1.0 for k in 0:3] : 2 .* rand(rng, 4) .- 1
            p = p0 .+ rad .* t; minimum(cvals(p)) > 0 || continue
            beta = maximum(tm_eval(B, t) - remval(B) for (nm, B) in bt.betas; init = 0.0)
            L0 = maximum(tm_eval(L, t) - remval(L) for (nm, L) in bt.L0s; init = 0.0)
            Kp, Mp, TT, TC, CC = s.pencil_at(p)
            kfix = max(kfix, kq(Kp, Mp, TT, TC, CC, s.C, s.dpsi, beta, L0))
            δ = rad .* t
            C1 = s.C .+ sum(δ[i] .* s.dC[:, :, i] for i in 1:4); psi1 = s.dpsi .+ sum(δ[i] .* s.ddpsi[:, i] for i in 1:4)
            kfo = max(kfo, kq(Kp, Mp, TT, TC, CC, C1, psi1, beta, L0))
        end
        line *= @sprintf("  r=%.2f: fixed %.3f  first-order %.3f", r, kfix, kfo)
    end
    println(line)
end
