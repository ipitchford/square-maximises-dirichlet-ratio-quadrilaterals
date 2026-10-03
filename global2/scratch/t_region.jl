include("../geom.jl"); include("../rules.jl"); include("../fsolve.jl")
using Random
rng = MersenneTwister(3)
for R in (0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
    cnt = Dict{Symbol,Int}(); mxa = 0.0; minwd = 1.0; minratio = Inf
    for _ in 1:200000
        p = [R*rand(rng), R*rand(rng), R*rand(rng), R*(2rand(rng)-1)]
        c = classify_point(p; cube = 0.049, eta = 0.035, s0 = 0.1)
        cnt[c] = get(cnt, c, 0) + 1
        if c == :R
            V = verts(p); mxa = max(mxa, maximum(angles_deg(V))); w, D = width_diam(V); minwd = min(minwd, w/D)
            cv = cvals(p); minratio = min(minratio, minimum(cv[i]/cv[mod1(i+2,4)] for i in 1:4))
        end
    end
    println("|p|inf <= $R: ", cnt, "  max angle in R: ", round(mxa, digits=1), "  min w/D: ", round(minwd, digits=3), " min c-ratio: ", round(minratio, digits=3))
end
