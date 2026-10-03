include("../geom.jl"); include("../rules.jl"); include("../fsolve.jl"); include("../tm.jl"); include("../assemble.jl"); include("../criterion.jl")
using Printf
# instrumented: global factors on the Taylor remainder (P) and on product truncation
const FP = Ref(1.0); const FH = Ref(1.0); const FR = Ref(1.0)
function taylor_tm!(h::TM{Float64}, D::Matrix{Float64}, AD::Matrix{Float64}, su::Int, sv::Int, d::Int,
                    hu::Vector{Float64}, hv::Vector{Float64}, QA, QB, QC, HU, HV)
    c = h.c
    c[1] = D[su+1, sv+1]
    D10 = D[su+2, sv+1]; D01 = D[su+1, sv+2]
    for i = 1:4; c[1+i] = D10 * hu[i] + D01 * hv[i]; end
    D20 = D[su+3, sv+1]; D11 = D[su+2, sv+2]; D02 = D[su+1, sv+3]
    for k = 6:NC; c[k] = D20 * QA[k] + D11 * QB[k] + D02 * QC[k]; end
    rem = 0.0; dd = d - su - sv
    for a = 0:dd, b = max(0, 3 - a):dd-a
        rem += AD[a+su+1, b+sv+1] * HU[a+1] * HV[b+1]
    end
    h.r[1] = FP[] * rem
    h
end
function _mulrem(f::TM{Float64}, g::TM{Float64})
    f0, fl, fq = _norms(f); g0, gl, gq = _norms(g)
    Pf = f0 + fl + fq; Pg = g0 + gl + gq
    return FR[] * (Pf * g.r[1] + Pg * f.r[1] + f.r[1] * g.r[1]) + FH[] * (fl * gq + fq * gl + fq * gq)
end
awF = newasm(Float64)
for (p0, n, r) in (([0.3,0.3,0.3,0.1], 8, 0.01), ([0.2,0.5,0.1,-0.3], 10, 0.002), ([0.0,0.1,0.05,0.0], 8, 0.01))
    sol = fsolve(p0, n)
    println("p0=$p0 n=$n r=$r  max|C|=", maximum(abs, sol.C), " max|dpsi|=", maximum(abs, sol.dpsi))
    for (fp, fh, label) in ((1.0, 1.0, "full"), (0.0, 1.0, "no P-taylor rem"), (1.0, 0.0, "no product trunc"), (0.0, 0.0, "neither (only rem propagation)"))
        FP[] = fp; FH[] = fh
        base = assemble_box!(awF, p0, fill(r, 4), sol)
        @printf("   %-32s rem: M11=%.3g K11=%.3g K22=%.3g K33=%.3g A2=%.3g | lin K22=%.3g quad K22=%.3g\n", label, remval(base.M[1]), remval(base.K[1]), remval(base.K[4]), remval(base.K[6]), remval(base.A2),
            sum(abs, base.K[4].c[2:5]), sum(abs, base.K[4].c[6:15]))
    end
    FP[] = 1.0; FH[] = 1.0
end
