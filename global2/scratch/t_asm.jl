include("../geom.jl"); include("../rules.jl"); include("../fsolve.jl"); include("../tm.jl"); include("../assemble.jl")
using Random, Printf, LinearAlgebra
function tm_eval(f::TM{Float64}, t)
    s = f.c[1]
    for i in 1:4; s += f.c[1+i]*t[i]; end
    for i in 1:4, j in i:4; s += f.c[QIDX[i,j]]*t[i]*t[j]; end
    s
end
p0 = [0.05, 0.2, 0.1, -0.05]
sol = fsolve(p0, 8)
awF = newasm(Float64)
r0 = assemble_box!(awF, p0, zeros(4), sol)
Mv = [x.c[1] for x in r0.M]; Kv = [x.c[1] for x in r0.K]
println("point box: M = ", round.(Mv, sigdigits=6), "\n   K = ", round.(Kv, sigdigits=8), " theta=", sol.theta[1:3])
println("   A2 = ", r0.A2.c[1], " vs fsolve ", sol.A2, "  rel diff ", (r0.A2.c[1]-sol.A2)/sol.A2)
println("   max remainder at point box: ", maximum(remval.(vcat(r0.M, r0.K, [r0.A2]))))
# box test
rad = [0.01, 0.02, 0.015, 0.01]
t0 = @elapsed rB = assemble_box!(awF, p0, rad, sol)
t0 = @elapsed rB = assemble_box!(awF, p0, rad, sol)
println("Float64 box assembly time: ", t0)
names = ["M11","M12","M13","M22","M23","M33","K11","K12","K13","K22","K23","K33","A2"]
TMs = vcat(deepcopy(rB.M), deepcopy(rB.K), [deepcopy(rB.A2)])
rng = MersenneTwister(1)
maxratio = zeros(13)
aw2 = newasm(Float64)
for s in 1:30
    t = 2 .* rand(rng, 4) .- 1
    if s <= 2; t = sign.(t); end
    p = p0 .+ rad .* t
    rp = assemble_box!(aw2, p, zeros(4), sol)
    vals = vcat([x.c[1] for x in rp.M], [x.c[1] for x in rp.K], [rp.A2.c[1]])
    for q in 1:13
        d = abs(tm_eval(TMs[q], t) - vals[q])
        maxratio[q] = max(maxratio[q], d / max(remval(TMs[q]), 1e-300))
    end
end
for q in 1:13
    @printf("%s: c0=%.6g lin=%.3g quad=%.3g rem=%.3g  max|TM(t)-exact|/rem = %.3g\n", names[q], TMs[q].c[1], sum(abs, TMs[q].c[2:5]), sum(abs, TMs[q].c[6:15]), remval(TMs[q]), maxratio[q])
end
# Arb enclosure test
awA = newasm(Arb)
tA = @elapsed rA = assemble_box!(awA, p0, rad, sol)
tA = @elapsed rA = assemble_box!(awA, p0, rad, sol)
println("Arb box assembly time: ", tA)
TA = vcat(rA.M, rA.K, [rA.A2])
for q in 1:13
    # check Float64 coefficients inside Arb balls (up to rounding) and remainders comparable
    dc = maximum(abs(Float64(Arblib.midref(TA[q].c[k])) - TMs[q].c[k]) for k in 1:15)
    rads = maximum(Float64(Arblib.radref(TA[q].c[k])) for k in 1:15)
    @printf("%s: Arb rem=%.3g (F64 rem %.3g), max coeff diff %.2g, max ball radius %.2g\n", names[q], remval(TA[q]), remval(TMs[q]), dc, rads)
end
