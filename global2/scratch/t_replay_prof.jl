include("../replay_t.jl")
using Printf
recs = [parse_rec(ln) for ln in Iterators.take(eachline("runs/mr_test/sample.jsonl"), 20)]
W = TWorker()
function prof(r, W)
    lo = r["lo"]; hi = r["hi"]; n = Int(r["n"]); sh = Int(r["shift"])
    p0, rad = box_center_radius(lo, hi)
    t = zeros(6)
    t[1] = @elapsed (sb = scalar_bounds(p0, rad; ws = W.wsA))
    t[2] = @elapsed (sol = tsolve1(p0, n; shift = sh))
    t[3] = @elapsed (RA = transport1_refs(Arb, sol))
    t[4] = @elapsed transport1_base!(W.tpA, RA, p0, rad, W.wsA, sh)
    t[5] = @elapsed (btA = bound_tms(Arb, p0, rad, W.wsA))
    t[6] = @elapsed crit_eval!(W.cpA, W.tpA.M, W.tpA.K, W.tpA.A2, W.wsA; mode = :LG, beta = 30.0, L0 = 10.0)
    t
end
prof(recs[1], W)
T = zeros(6)
for r in recs; T .+= prof(r, W); end
@printf("per certificate (s): scalar_bounds %.3f  tsolve1 %.3f  refs(Arb) %.3f  base(Arb) %.3f  bound_tms(Arb) %.3f  crit %.3f  total %.3f\n", (T ./ length(recs))..., sum(T) / length(recs))
# exclusion check timing: E2 with k up to 8
