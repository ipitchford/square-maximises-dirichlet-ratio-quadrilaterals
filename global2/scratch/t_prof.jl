include("../cover_t.jl")
using Printf
W = TWorker()
p0 = [0.3, 0.2, 0.4, 0.1]; rad = fill(0.03, 4); sh = 1
s = tsolve1(p0, 8; shift = sh)
for rep in 1:2
    t1 = @elapsed RA = transport1_refs(Arb, s)
    t2 = @elapsed bA = transport1_base!(W.tpA, RA, p0, rad, W.wsA, sh)
    t3 = @elapsed btA = bound_tms(Arb, p0, rad, W.wsA)
    t4 = @elapsed sb = scalar_bounds(p0, rad; ws = W.wsA)
    t5 = @elapsed r = crit_eval!(W.cpA, W.tpA.M, W.tpA.K, W.tpA.A2, W.wsA; mode = :LG, beta = btA.betas[1][2], L0 = 0.0)
    t6 = @elapsed (s2 = tsolve1(p0, 8; shift = sh))
    t7 = @elapsed (RF = transport1_refs(Float64, s))
    t8 = @elapsed bt = bound_tms(Float64, p0, rad, W.wsF)
    t9 = @elapsed transport1_base!(W.tpF, RF, p0, rad, W.wsF, sh)
    t10 = @elapsed crit_eval!(W.cpF, W.tpF.M, W.tpF.K, W.tpF.A2, W.wsF; mode = :LG, beta = bt.betas[1][2], L0 = 0.0)
    @printf("Arb: refs %.3f base %.3f bounds %.3f scalar %.3f crit %.4f | F64: tsolve1 %.3f refs %.3f bounds %.3f base %.4f crit %.4f\n", t1, t2, t3, t4, t5, t6, t7, t8, t9, t10)
end
