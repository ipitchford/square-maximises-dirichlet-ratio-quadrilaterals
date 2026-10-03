include("../geom.jl"); include("../rules.jl"); include("../fsolve.jl"); include("../tm.jl"); include("../assemble.jl"); include("../criterion.jl")
using Printf
set_order!(4)
awF = newasm(Float64); cpF = critpool(Float64); wsF = newwork(Float64)
degnorms(f) = [sum(abs, f.c[TAB[].degrange[d+1]]) for d in 0:TAB[].order]
p0 = [0.5,0.3,0.2,0.3]; n = 10
sol = fsolve(p0, n)
println("theta=", round.(sol.theta[1:4],digits=3), " maxangle=", maximum(angles_deg(verts(p0))), " cvals=", cvals(p0))
for r in (1e-14, 0.005)
    rad = fill(r, 4)
    base = assemble_box!(awF, p0, rad, sol)
    sbb = scalar_bounds(p0, rad)
    L0b = max(sbb.L_rect, sbb.L_hersch); betab = max(sbb.beta_rect_ks, L0b + sbb.gapAC)
    println("r=$r: scalar bounds ", sbb, " beta=", betab)
    res = crit_eval!(cpF, base.M, base.K, base.A2, wsF; mode = :LG, beta = betab, L0 = L0b, kappa = 2//1)
    println("  kappa=2 res: ", res)
    for (nm, f) in (("M11", base.M[1]), ("K11", base.K[1]), ("K22", base.K[4]), ("A2", base.A2), ("Mn11", cpF.Mn[1]), ("Kn11", cpF.Kn[1]), ("Kn22", cpF.Kn[4]), ("A2n", cpF.A2n), ("g1", cpF.g1), ("g2", cpF.g2), ("At11", cpF.At[1]), ("At22", cpF.At[4]), ("s11", cpF.s11))
        lo, hi = tm_bounds(f)
        @printf("   %-5s c0=% .4e  degnorms=%s rem=%.3e  range=[% .4e, % .4e]\n", nm, f.c[1], string(round.(degnorms(f)[2:end], sigdigits=3)), remval(f), lo, hi)
    end
end
