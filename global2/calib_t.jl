# calib_t.jl — rigorous (Arb) cost and Float64→Arb agreement of the transport certificate on pilot-passing boxes.
# usage: julia -t 1 calib_t.jl <pilotdir> <count> <outfile>
include("cover_t.jl")
using Random
pdir = ARGS[1]; cnt = parse(Int, ARGS[2]); outf = ARGS[3]
cands = Tuple{String,Vector{Float64},Vector{Float64}}[]
for f in readdir(pdir; join = true)
    endswith(f, ".jsonl") || continue
    for ln in eachline(f)
        occursin("pilot_pass", ln) || continue
        push!(cands, (match(r"\"path\":\"([^\"]+)\"", ln).captures[1],
                      parse.(Float64, split(match(r"\"lo\":\[([^\]]+)\]", ln).captures[1], ",")),
                      parse.(Float64, split(match(r"\"hi\":\[([^\]]+)\]", ln).captures[1], ","))))
    end
end
rng = MersenneTwister(7)
picked = cands[randperm(rng, length(cands))[1:min(cnt, length(cands))]]
W = TWorker()
process_transport(picked[1][1], picked[1][2], picked[1][3], W, false)      # warm-up
open(outf, "w") do io
    println(io, "path status t_total t_f64 kappa_cert bound")
    for (path, lo, hi) in picked
        t = @elapsed rec, dim = process_transport(path, lo, hi, W, false)
        st = rec === nothing ? "split" : rec["status"]
        @printf(io, "%s %s %.3f %.3f %.4f %s\n", path, st, t, rec === nothing ? NaN : get(rec, "tF", NaN), rec === nothing ? NaN : get(rec, "kappa_cert", NaN), rec === nothing ? "-" : get(rec, "bound", "-"))
    end
end
