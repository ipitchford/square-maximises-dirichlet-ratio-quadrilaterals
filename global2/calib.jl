# calib.jl — calibrate the rigorous (Arb) cost and the Float64→Arb agreement on boxes that passed the pilot.
# usage: julia -t 1 calib.jl <pilotdir> <count> <outfile>
include("cover.jl")
using Random

pdir = ARGS[1]; cnt = parse(Int, ARGS[2]); outf = ARGS[3]
cands = Tuple{String,Vector{Float64},Vector{Float64},Int}[]
for f in readdir(pdir; join = true)
    endswith(f, ".jsonl") || continue
    for ln in eachline(f)
        occursin("pilot_pass", ln) || continue
        path = match(r"\"path\":\"([^\"]+)\"", ln).captures[1]
        lo = parse.(Float64, split(match(r"\"lo\":\[([^\]]+)\]", ln).captures[1], ","))
        hi = parse.(Float64, split(match(r"\"hi\":\[([^\]]+)\]", ln).captures[1], ","))
        n = parse(Int, match(r"\"n\":(\d+)", ln).captures[1])
        push!(cands, (path, lo, hi, n))
    end
end
rng = MersenneTwister(7)
# stratify by n so every degree is timed
picked = []
for nn in (8, 10, 12)
    c = filter(x -> x[4] == nn, cands)
    isempty(c) && continue
    append!(picked, c[randperm(rng, length(c))[1:min(length(c), cld(cnt, 3))]])
end
W = Worker()
# warm-up (compilation)
process(picked[1][1], picked[1][2], picked[1][3], W, false)
open(outf, "w") do io
    println(io, "path n status t_total t_f64 kappa_cert")
    for (path, lo, hi, n) in picked
        t = @elapsed rec, dim = process(path, lo, hi, W, false)
        st = rec === nothing ? "split(arb fail)" : rec["status"]
        tf = rec === nothing ? NaN : get(rec, "tF", NaN)
        kc = rec === nothing ? NaN : get(rec, "kappa_cert", NaN)
        @printf(io, "%s %d %s %.3f %.3f %.5f\n", path, n, st, t, tf, kc)
        @printf(stderr, "%s n=%d %s t=%.2fs (f64 %.2fs) kappa=%.4f\n", path, n, st, t, tf, kc)
    end
end
