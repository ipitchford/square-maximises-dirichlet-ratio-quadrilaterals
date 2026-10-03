include("../geom.jl")
using Random
rng = MersenneTwister(17)
lo = [0.0,0,0,-0.905]; hi = [0.905,1.4,1.4,0.905]; volbox = prod(hi .- lo)
N = 2_000_000; cR = 0; layer = Dict(0.05=>0, 0.08=>0, 0.12=>0, 0.2=>0); ang = Dict(160=>0, 165=>0, 170=>0)
for _ in 1:N
    p = lo .+ rand(rng, 4) .* (hi .- lo)
    classify_point(p; cube = 0.049, eta = 0.035, s0 = 0.1) == :R || continue
    global cR += 1
    c = cvals(p); r = minimum(c[i]/c[mod1(i+2,4)] for i in 1:4)
    for k in keys(layer); r < k && (layer[k] += 1); end
    a = maximum(angles_deg(verts(p)))
    for k in keys(ang); a > k && (ang[k] += 1); end
end
volR = volbox * cR / N
println("vol(R) = ", round(volR, digits=4))
for k in sort(collect(keys(layer))); println("  c-ratio < $k: volume ", round(volbox*layer[k]/N, sigdigits=3), "  (", round(100layer[k]/cR, digits=1), "% of R)"); end
for k in sort(collect(keys(ang))); println("  max angle > $k deg: volume ", round(volbox*ang[k]/N, sigdigits=3), "  (", round(100ang[k]/cR, digits=1), "% of R)"); end
for (k, w) in ((0.08, 0.005), (0.12, 0.01))
    v = volbox*layer[k]/N
    println("  boxes if c-ratio<$k needs width $w: ", round(v / w^4, sigdigits=3))
end
