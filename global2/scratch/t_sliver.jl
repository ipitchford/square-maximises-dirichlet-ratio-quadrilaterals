include("../cover.jl")
ln = filter(l -> occursin("pilot_pass", l), readlines(joinpath(@__DIR__, "..", "runs", "shell025_pilot.jsonl")))[end]
lo = parse.(Float64, split(match(r"\"lo\":\[([^\]]+)\]", ln).captures[1], ","))
hi = parse.(Float64, split(match(r"\"hi\":\[([^\]]+)\]", ln).captures[1], ","))
p0, rad = box_center_radius(lo, hi)
println("p0=", p0, " rad=", rad, " angles=", round.(angles_deg(verts(p0)), digits=1))
for scale in (1.0, 4.0, 16.0, 64.0)
    r = copy(rad); r[1] *= scale; r[4] *= scale
    sb = scalar_bounds(p0, r); L0 = max(sb.L_rect, sb.L_hersch); beta = max(sb.beta_rect_ks, addlo(L0, sb.gapAC))
    s = fsolve_d(p0, 8)
    W = Worker()
    b = assemble_box!(W.awF, p0, r, s)
    res = crit_eval!(W.cpF, b.M, b.K, b.A2, W.wsF; mode = :LG, beta, L0)
    println("a,d radius x$scale: beta=", round(beta, digits=3), " (KS/rect ", round(sb.beta_rect_ks,digits=3), " L0+gap ", round(L0+sb.gapAC,digits=3), ") rho=", round(s.rho0,digits=3), " ok=", res.ok, " which=", res.which, " g2lo=", round(res.g2lo, digits=4), " g1lo=", round(res.g1lo,digits=4), " infl s11=", round.(influence(W.cpF.s11), sigdigits=3), " infl g2=", round.(influence(W.cpF.g2), sigdigits=3))
end
