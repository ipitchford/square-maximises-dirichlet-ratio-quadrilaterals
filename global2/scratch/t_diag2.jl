include("../geom.jl"); include("../rules.jl"); include("../fsolve.jl"); include("../tm.jl"); include("../assemble.jl"); include("../criterion.jl")
using Printf
const FP = Ref(1.0); const FH = Ref(1.0); const FPS = Ref(1.0)
# instrumented Float64 taylor (P vs psi distinguished by degree d > n)
const NDEG = Ref(8)
function taylor_tm!(h::TM{Float64}, D::Matrix{Float64}, AD::Matrix{Float64}, su::Int, sv::Int, d::Int,
                    HUV::Matrix{TM{Float64,Float64}}, HU, HV)
    tab = TAB[]; K = tab.order
    fill!(h.c, 0.0)
    for a = 0:K, b = 0:K-a
        Dab = D[a+su+1, b+sv+1]; Dab == 0.0 && continue
        hc = HUV[a+1, b+1].c
        for k in tab.degrange[a+b+1]; h.c[k] += Dab * hc[k]; end
    end
    rem = 0.0; dd = d - su - sv
    for a = 0:dd, b = max(0, K + 1 - a):dd-a
        rem += AD[a+su+1, b+sv+1] * HU[a+1] * HV[b+1]
    end
    h.r[1] = (d > NDEG[] ? FPS[] : FP[]) * rem
    h
end
function _mulrem(f::TM{Float64}, g::TM{Float64}, ws::TMWork{Float64})
    K = TAB[].order
    nf = _degnorms!(ws.nf, f); ng = _degnorms!(ws.ng, g)
    Pf = sum(nf[1:K+1]); Pg = sum(ng[1:K+1])
    drop = 0.0
    for d1 = 1:K, d2 = K+1-d1:K; drop += nf[d1+1] * ng[d2+1]; end
    return Pf * g.r[1] + Pg * f.r[1] + f.r[1] * g.r[1] + FH[] * drop
end
for K in (4,)
    set_order!(K)
    awF = newasm(Float64)
    for (p0, n, r) in (([0.5,0.3,0.2,0.3], 10, 0.005), ([0.2, 0.5, 0.1, -0.3], 12, 0.005), ([0.3,0.3,0.3,0.1], 8, 0.02))
        sol = fsolve(p0, n); NDEG[] = n
        println("p0=$p0 n=$n r=$r: A2=", sol.A2, " max|C1|=", round(maximum(abs, sol.C[:,1]),sigdigits=3), " max|dpsi|=", round(maximum(abs, sol.dpsi), sigdigits=3), " A2_nopsi=", sol.A2_nopsi)
        for (fp, fps, fh, lab) in ((1,1,1,"full"), (0,1,1,"no P tails"), (1,0,1,"no psi tails"), (1,1,0,"no product drop"), (0,0,1,"no tails at all"))
            FP[] = fp; FPS[] = fps; FH[] = fh
            base = assemble_box!(awF, p0, fill(r,4), sol)
            @printf("   %-16s rem A2/A2=%.3g  K11/K11=%.3g  K22/K22=%.3g M11=%.3g\n", lab, remval(base.A2)/base.A2.c[1], remval(base.K[1])/base.K[1].c[1], remval(base.K[4])/base.K[4].c[1], remval(base.M[1]))
        end
        FP[] = 1; FPS[] = 1; FH[] = 1
    end
end
