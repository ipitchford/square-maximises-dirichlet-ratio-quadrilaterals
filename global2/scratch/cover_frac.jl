using Random
f = ARGS[1]
boxes = Tuple{Vector{Float64},Vector{Float64}}[]
for ln in eachline(f)
    (occursin("pilot_pass", ln) || occursin("\"certified\"", ln)) || continue
    lo = parse.(Float64, split(match(r"\"lo\":\[([^\]]+)\]", ln).captures[1], ","))
    hi = parse.(Float64, split(match(r"\"hi\":\[([^\]]+)\]", ln).captures[1], ","))
    push!(boxes, (lo, hi))
end
println("boxes: ", length(boxes))
rng = MersenneTwister(1)
for rho in (0.1, 0.15, 0.175, 0.2, 0.225, 0.25)
    N = 4000; cov = 0; tot = 0; nb = count(b -> all(b[1][i] < (i == 4 ? rho : rho) for i in 1:3) && b[1][4] < rho && b[2][4] > -rho, boxes)
    while tot < N
        p = [rho*rand(rng), rho*rand(rng), rho*rand(rng), rho*(2rand(rng)-1)]
        maximum(abs, p) <= 0.049 && continue
        tot += 1
        any(b -> all(b[1][i] <= p[i] <= b[2][i] for i in 1:4), boxes) && (cov += 1)
    end
    println("rho=$rho: covered fraction ", round(cov/tot, digits=4), "  boxes intersecting: ", nb)
end
