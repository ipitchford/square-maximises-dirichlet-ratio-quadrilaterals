# pilot_summary.jl — statistics and cost prediction from the per-cell pilot logs (runs/pilot1/cell_*.jsonl)
using Printf, Statistics
dir = ARGS[1]; ncells = parse(Int, ARGS[2]); cap = parse(Int, ARGS[3])
tarb = Dict(8 => parse(Float64, ARGS[4]), 10 => parse(Float64, ARGS[5]), 12 => parse(Float64, ARGS[6]))   # Arb s/box (1 thread)
rows = []
for f in readdir(dir; join = true)
    m = match(r"cell_(\d+)\.jsonl$", f); m === nothing && continue
    lines = readlines(f)
    length(lines) == 0 && continue
    st = Dict{String,Int}(); npass = Dict(8 => 0, 10 => 0, 12 => 0); tF = 0.0; kmax = 0.0
    for ln in lines
        s = match(r"\"status\":\"([A-Za-z0-9_]+)\"", ln).captures[1]
        st[s] = get(st, s, 0) + 1
        if s == "pilot_pass"
            n = parse(Int, match(r"\"n\":(\d+)", ln).captures[1]); npass[n] += 1
            tF += parse(Float64, match(r"\"t\":([0-9.eE+-]+)", ln).captures[1])
            kmax = max(kmax, parse(Float64, match(r"\"kappa_f64\":([0-9.eE+-]+)", ln).captures[1]))
        end
    end
    capped = length(lines) >= cap
    push!(rows, (cell = parse(Int, m.captures[1]), n = length(lines), st, npass, capped, tF, kmax))
end
sort!(rows, by = r -> r.cell)
K = length(rows); ncap = count(r -> r.capped, rows)
leaves(r) = sum(v for (k, v) in r.st if k != "split"; init = 0)
passes(r) = get(r.st, "pilot_pass", 0)
@printf("sampled cells: %d of %d coarse cells; capped at %d processed boxes: %d (%.0f%%)\n", K, ncells, cap, ncap, 100ncap / K)
for r in rows
    @printf("  cell %2d: processed %5d  leaves %5d  pass %4d (n8 %d, n10 %d, n12 %d)  excl %4d  %s  max kappa_f64 %.4f\n", r.cell, r.n, leaves(r), passes(r),
            r.npass[8], r.npass[10], r.npass[12], leaves(r) - passes(r), r.capped ? "CAPPED" : "done  ", r.kmax)
end
done = filter(r -> !r.capped, rows)
mean_done_pass = mean(passes.(done)); mean_done_proc = mean(getfield.(done, :n))
# lower bound: capped cells counted at the cap with their observed pass fraction
lb_proc = ncells * mean(min.(getfield.(rows, :n), cap))
lb_pass = ncells * mean(passes.(rows))
@printf("finished cells: mean processed %.0f, mean certified-candidate leaves %.0f\n", mean_done_proc, mean_done_pass)
@printf("LOWER BOUND (capped cells counted at the cap): processed boxes >= %.3g, certificate leaves >= %.3g\n", lb_proc, lb_pass)
avg_tarb = sum(tarb[n] * sum(r.npass[n] for r in rows) for n in (8, 10, 12)) / max(1, sum(passes.(rows)))
@printf("mean Arb time per certified box (1 thread): %.2f s\n", avg_tarb)
