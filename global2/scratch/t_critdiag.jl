include("../geom.jl"); include("../rules.jl"); include("../fsolve.jl"); include("../tm.jl"); include("../assemble.jl"); include("../criterion.jl")
using Printf
set_order!(parse(Int, get(ENV, "ORD", "3")))
awF = newasm(Float64); cpF = critpool(Float64); wsF = newwork(Float64)
function degnorms(f)
    [sum(abs, f.c[TAB[].degrange[d+1]]) for d in 0:TAB[].order]
end
for (p0, n, r) in (([0.3, 0.3, 0.3, 0.1], 8, 0.02), ([0.0, 0.06, 0.0, 0.0], 8, 0.01))
    sol = fsolve(p0, n)
    println("p0=$p0 theta=", round.(sol.theta[1:4], digits=4))
    rad = fill(r, 4)
    base = assemble_box!(awF, p0, rad, sol)
    sbb = scalar_bounds(p0, rad)
    L0b = max(sbb.L_rect, sbb.L_hersch); betab = max(sbb.beta_rect_ks, L0b + sbb.gapAC)
    res = crit_eval!(cpF, base.M, base.K, base.A2, wsF; mode = :LG, beta = betab, L0 = L0b)
    println("  beta=$betab res: ", res)
    for (nm, f) in (("M11", base.M[1]), ("K11", base.K[1]), ("K22", base.K[4]), ("K33", base.K[6]), ("M22", base.M[4]), ("A2", base.A2), ("g1", cpF.g1), ("g2", cpF.g2), ("At11", cpF.At[1]), ("At22", cpF.At[4]), ("At33", cpF.At[6]), ("At12", cpF.At[2]), ("At23", cpF.At[5]), ("s11", cpF.s11), ("s22", cpF.s22), ("s12", cpF.s12), ("tr", cpF.tr), ("det", cpF.det))
        lo, hi = tm_bounds(f)
        @printf("   %-5s c0=% .4e  degnorms=%s rem=%.3e  range=[% .4e, % .4e]\n", nm, f.c[1], string(round.(degnorms(f)[2:end], sigdigits=3)), remval(f), lo, hi)
    end
end
