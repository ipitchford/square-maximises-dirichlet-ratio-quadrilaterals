# pilot_t.jl — Float64 pilot of the transport cover (cover_t.jl --pilot logic on a random sample of coarse cells).
# usage: julia -t 8 pilot_t.jl <K> <cap> <outdir>
include("cover_t.jl")
using Random, Statistics
function coarse_cells_t()
    cells = Tuple{String,Vector{Float64},Vector{Float64}}[]; excl = Dict{String,Int}()
    stack = root_boxes_t(); W = newwork(Arb)
    while !isempty(stack)
        path, lo, hi = pop!(stack)
        if is_shell_box(lo, hi); excl["SHELL"] = get(excl, "SHELL", 0) + 1; continue; end
        ex = exclusion_t(lo, hi; ws = W)
        if ex !== :none; excl[String(ex)] = get(excl, String(ex), 0) + 1; continue; end
        w = hi .- lo
        if maximum(w) > WTRY_T
            c1, c2 = children(path, lo, hi, argmax(w)); push!(stack, c1, c2)
        else
            push!(cells, (path, lo, hi))
        end
    end
    cells, excl
end
K = parse(Int, ARGS[1]); cap = parse(Int, ARGS[2]); outdir = ARGS[3]; mkpath(outdir)
t0 = time(); cells, excl = coarse_cells_t()
vol(c) = prod(c[3] .- c[2])
msg = @sprintf("coarse cells: %d (volume %.4f), excluded at coarse level: %s, %.1fs", length(cells), sum(vol, cells), string(excl), time() - t0)
println(msg); flush(stdout)
open(joinpath(outdir, "coarse.txt"), "w") do io; println(io, msg); end
rng = MersenneTwister(20261001)
sample = cells[randperm(rng, length(cells))[1:min(K, length(cells))]]
for (j, c) in enumerate(sample)
    logf = joinpath(outdir, "cell_$(j).jsonl"); rm(logf; force = true)
    tc = time()
    stats = run_cover(logf; pilot = true, maxboxes = cap, roots = [c], procfun = process_transport, mkworker = TWorker)
    @printf("cell %d/%d %s: %s %.1fs\n", j, length(sample), c[1], string(stats), time() - tc); flush(stdout)
end
println("done ", time() - t0)
