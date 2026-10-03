# pilot.jl — Float64 pilot of the whole cover (no rigour): predicts the number of boxes and the rigorous time.
#
# 1. Coarse level: starting from the a-priori box (split at the cube faces), apply the rigorous exclusion tests
#    (E1–E4, Arb) and bisect the widest side until every surviving cell has all widths ≤ WTRY (where the cover
#    starts attempting certificates).  This is exactly what cover.jl does above the attempt level.
# 2. Sample K surviving cells uniformly at random and run the full adaptive Float64 pilot (cover.jl --pilot
#    logic) inside each, with a cap on processed boxes per cell.
# 3. Extrapolate: total leaves ≈ N_cells · mean(leaves per cell); rigorous time ≈ Σ (Arb time per certified box,
#    calibrated by calib.jl) + Float64 time, divided by the thread count.
#
# usage: julia -t 8 pilot.jl <K> <cap> <outdir>

include("cover.jl")
using Random, Statistics

function coarse_cells()
    cells = Tuple{String,Vector{Float64},Vector{Float64}}[]
    excl = Dict{String,Int}()
    stack = root_boxes()
    W = newwork(Arb)
    while !isempty(stack)
        path, lo, hi = pop!(stack)
        ex = exclusion(lo, hi; cube = CUBE, ws = W)
        if ex !== :none
            excl[String(ex)] = get(excl, String(ex), 0) + 1
            continue
        end
        w = hi .- lo
        if maximum(w) > WTRY
            c1, c2 = children(path, lo, hi, argmax(w)); push!(stack, c1, c2)
        else
            push!(cells, (path, lo, hi))
        end
    end
    return cells, excl
end

if abspath(PROGRAM_FILE) == @__FILE__
    K = parse(Int, ARGS[1]); cap = parse(Int, ARGS[2]); outdir = ARGS[3]
    mkpath(outdir)
    t0 = time()
    cells, excl = coarse_cells()
    vol(c) = prod(c[3] .- c[2])
    @printf(stderr, "coarse cells: %d (volume %.4f), excluded at coarse level: %s, %.1fs\n", length(cells), sum(vol, cells), string(excl), time() - t0)
    rng = MersenneTwister(20261001)
    sample = cells[randperm(rng, length(cells))[1:min(K, length(cells))]]
    rows = []
    for (j, c) in enumerate(sample)
        logf = joinpath(outdir, "cell_$(j).jsonl")
        rm(logf; force = true)
        tc = time()
        stats = run_cover(logf; pilot = true, maxboxes = cap, roots = [c])
        el = time() - tc
        # per-n counts and Float64 times of passing leaves
        npass = Dict(8 => 0, 10 => 0, 12 => 0); tF = 0.0; kmin = Inf
        for ln in eachline(logf)
            occursin("pilot_pass", ln) || continue
            n = parse(Int, match(r"\"n\":(\d+)", ln).captures[1]); npass[n] += 1
            tF += parse(Float64, match(r"\"t\":([0-9.eE+-]+)", ln).captures[1])
            kmin = min(kmin, 2.5 - parse(Float64, match(r"\"kappa_f64\":([0-9.eE+-]+)", ln).captures[1]))
        end
        hit = get(stats, "split", 0) + sum(v for (k, v) in stats if k != "split"; init = 0) >= cap
        push!(rows, (cell = c[1], lo = c[2], hi = c[3], stats, npass, elapsed = el, hitcap = hit, minmargin = kmin))
        @printf(stderr, "cell %d/%d %s: %s  n-split %s  %.1fs %s\n", j, length(sample), c[1], string(stats), string(npass), el, hit ? "CAP HIT" : "")
    end
    open(joinpath(outdir, "pilot_cells.txt"), "w") do io
        println(io, "# coarse cells: ", length(cells), " total volume ", sum(vol, cells), "; excluded at coarse level: ", excl)
        for r in rows
            println(io, r.cell, " lo=", r.lo, " hi=", r.hi, " stats=", r.stats, " npass=", r.npass, " elapsed=", round(r.elapsed, digits = 1), " hitcap=", r.hitcap, " minmargin=", r.minmargin)
        end
    end
    println("done ", time() - t0)
end
