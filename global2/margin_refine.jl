# margin_refine.jl — for the certified boxes with the largest logged κ*, find the certified κ* to high resolution
# (bisection on κ in (2.5 − 2^-3, 2.5] with 40 steps; each evaluation is a full Arb criterion on the same base TMs).
include("cover.jl")
logf = ARGS[1]; K = parse(Int, ARGS[2])
recs = [ln for ln in eachline(logf) if occursin("\"certified\"", ln)]
kap(ln) = parse(Float64, match(r"\"kappa_cert\":([0-9.eE+-]+)", ln).captures[1])
sort!(recs, by = ln -> -kap(ln))
W = Worker()
open(joinpath(dirname(logf), "margin_refine.txt"), "w") do io
    println(io, "path  logged_kappa  refined_kappa*  margin(5/2-kappa*)  which")
    for ln in recs[1:min(K, end)]
        path = match(r"\"path\":\"([^\"]+)\"", ln).captures[1]
        lo = parse.(Float64, split(match(r"\"lo\":\[([^\]]+)\]", ln).captures[1], ","))
        hi = parse.(Float64, split(match(r"\"hi\":\[([^\]]+)\]", ln).captures[1], ","))
        n = parse(Int, match(r"\"n\":(\d+)", ln).captures[1]); beta = parse(Float64, match(r"\"beta\":([0-9.eE+-]+)", ln).captures[1])
        mode = Symbol(match(r"\"mode\":\"([A-Z]+)\"", ln).captures[1])
        p0, rad = box_center_radius(lo, hi)
        sb = scalar_bounds(p0, rad; ws = W.wsA); L0 = max(sb.L_rect, sb.L_hersch)
        sol = fsolve_d(p0, n); base = assemble_box!(W.awA, p0, rad, sol)
        r0 = crit_eval!(W.cpA, base.M, base.K, base.A2, W.wsA; mode, beta, L0, kappa = 5 // 2)
        lo_k = big(5) / 2 - big(1) / 8; hi_k = big(5) / 2
        for _ = 1:40
            mid = (lo_k + hi_k) / 2
            rr = crit_eval!(W.cpA, base.M, base.K, base.A2, W.wsA; mode, beta, L0, kappa = rationalize(BigInt, mid; tol = big(1e-30)))
            rr.ok ? (hi_k = mid) : (lo_k = mid)
        end
        @printf(io, "%s  %.6f  %.10f  %.3e  %s ok@5/2=%s\n", path, kap(ln), Float64(hi_k), Float64(big(5) / 2 - hi_k), r0.which, r0.ok)
        @printf(stderr, "%s  logged %.6f  refined %.10f  margin %.3e\n", path, kap(ln), Float64(hi_k), Float64(big(5) / 2 - hi_k))
    end
end
