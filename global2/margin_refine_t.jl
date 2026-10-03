# margin_refine_t.jl — certified κ* to high resolution for the tightest transport certificates.
#
# usage: julia -t 8 margin_refine_t.jl <log> <K> [outfile]
#        julia -t 4 margin_refine_t.jl <log> <K> --from-replay <replay_results.jsonl>   (fine phase only, from the replay's κ*)
#
# The logged `kappa_cert` is the first rung of a short ladder above the Float64 estimate that the Arb certificate
# accepted (2.5 if no rung was accepted; κ = 5/2 itself is always certified with strict inequalities).  Here every
# selected certificate is rebuilt from its log record alone (`rebuild_cert_t` of replay_t.jl: Float64 choice, Arb
# reference integrals and Taylor models, logged lower-bound candidate in Arb) and κ is bisected in EXACT rationals.
# The invariant "the Arb criterion holds at hi" makes every reported κ* a certified upper bound of ξ = λ2/λ1 on the
# whole box, whether or not the test is monotone in κ.
# Selected: every certificate with kappa_cert = 2.5 and the K certificates with the largest kappa_cert < 2.5.
# Phase A: 14 bisection steps in [1, 5/2] for all selected boxes; phase B: 40 further steps for the 40 boxes with the
# largest phase-A κ*.  The smallest certified margin of the log is min(refined margins, 5/2 − largest unrefined kappa_cert).
include("replay_t.jl")
logf = ARGS[1]; K = parse(Int, ARGS[2])
outf = length(ARGS) >= 3 ? ARGS[3] : joinpath(dirname(logf), "margin_refine_t.txt")
certs = [r for r in (parse_rec(ln) for ln in eachline(logf)) if get(r, "status", "") == "certified"]
top = filter(r -> r["kappa_cert"] >= 2.5, certs)
rest = sort(filter(r -> r["kappa_cert"] < 2.5, certs), by = r -> -r["kappa_cert"])
sel = vcat(top, rest[1:min(K, end)])
unref = length(rest) > K ? rest[K+1]["kappa_cert"] : -Inf
@printf(stderr, "%d certified, %d with the ladder fallback 2.5, refining %d\n", length(certs), length(top), length(sel))

"rebuild the certificate of record r in W, then bisect κ in [lo, hi] (exact rationals; hi must be certified)"
function refine(r, W, lo::Rational{BigInt}, hi::Rational{BigInt}, steps)
    ok, msg, mode, betaA, L0A = try
        rebuild_cert_t(r, W)
    catch err
        (false, "error: " * first(sprint(showerror, err), 100), :LG, nothing, nothing)
    end
    ok || return (nothing, nothing, "REBUILD FAILED: " * msg)
    crit(k) = crit_eval!(W.cpA, W.tpA.M, W.tpA.K, W.tpA.A2, W.wsA; mode, beta = betaA, L0 = L0A, kappa = k).ok
    crit(hi) || return (nothing, nothing, "FAILS AT UPPER END " * string(Float64(hi)))
    for _ = 1:steps
        mid = (lo + hi) / 2
        crit(mid) ? (hi = mid) : (lo = mid)
    end
    return (lo, hi, r["bound"] * "/" * r["mode"])
end
function prun(items, f)
    out = Vector{Any}(undef, length(items))
    chunks = collect(Iterators.partition(eachindex(items), max(1, cld(length(items), 4 * Threads.nthreads()))))
    Threads.@threads for ch in chunks
        W = TWorker()
        for i in ch; out[i] = f(items[i], W); end
    end
    out
end
if (ir = findfirst(==("--from-replay"), ARGS)) !== nothing
    # fine refinement of the K certificates with the largest κ* recorded by the replay (14-step grid on [2, 5/2]):
    # 40 bisection steps in [κ* − 2^-15, κ*] (the lower end failed or is below the grid; hi = κ* is certified)
    resf = ARGS[ir+1]; outf2 = joinpath(dirname(logf), "margin_refine_t.txt")
    recs = Dict{String,Dict{String,Any}}()
    for ln in eachline(logf); r = parse_rec(ln); haskey(r, "path") && (recs[r["path"]] = r); end
    ks = Tuple{Float64,String}[]
    for ln in eachline(resf)
        endswith(ln, "\"end\":1}") || continue
        q = parse_rec(ln); k = get(q, "kstar", NaN); k >= 0 && push!(ks, (k, q["path"]))
    end
    ks = sort(unique(ks), rev = true)[1:min(K, end)]
    step = Rational{BigInt}(1, 2^15)
    out = Vector{Any}(undef, length(ks))
    Threads.@threads for i in eachindex(ks)
        W = TWorker(); k, path = ks[i]
        hi = Rational{BigInt}(BigFloat(k; precision = 256)); lo = hi - step
        a = refine(recs[path], W, lo, hi, 40)
        out[i] = (path, k, a[2] === nothing ? NaN : Float64(BigFloat(a[2]; precision = 256), RoundUp), a[3])
    end
    sort!(out, by = x -> -x[3])
    open(outf2, "w") do io
        @printf(io, "# fine refinement (40 bisection steps) of the %d certificates with the largest replay κ* (%s)\n", length(out), resf)
        @printf(io, "# SMALLEST CERTIFIED MARGIN 5/2 - kappa* = %.4e at %s\n", 2.5 - out[1][3], out[1][1])
        println(io, "path  replay_kappa*(14 steps)  refined_kappa*  margin  bound/mode")
        for x in out; @printf(io, "%s  %.10f  %.14f  %.4e  %s\n", x[1], x[2], x[3], 2.5 - x[3], x[4]); end
    end
    @printf(stderr, "smallest certified margin %.4e at %s\n", 2.5 - out[1][3], out[1][1])
    exit(0)
end
t0 = time()
A = prun(sel, (r, W) -> refine(r, W, Rational{BigInt}(1), Rational{BigInt}(5, 2), 14))
tA = time() - t0
okA = [i for i in eachindex(sel) if A[i][2] !== nothing]
order = sort(okA, by = i -> -A[i][2])
fine = order[1:min(40, end)]
B = prun(fine, (i, W) -> refine(sel[i], W, A[i][1], A[i][2], 40))
res = Dict{Int,Any}()
for i in eachindex(sel)
    a = A[i]
    res[i] = (sel[i]["path"], sel[i]["kappa_cert"], a[2] === nothing ? NaN : Float64(BigFloat(a[2]; precision = 256), RoundUp), a[3], "A")
end
for (j, i) in enumerate(fine)
    b = B[j]
    b[2] === nothing && (res[i] = (sel[i]["path"], sel[i]["kappa_cert"], NaN, b[3], "B"); continue)
    res[i] = (sel[i]["path"], sel[i]["kappa_cert"], Float64(BigFloat(b[2]; precision = 256), RoundUp), b[3], "B")
end
rows = sort(collect(values(res)), by = x -> isnan(x[3]) ? -Inf : -x[3])
nbad = count(x -> isnan(x[3]), rows)
kmax = maximum((x[3] for x in rows if !isnan(x[3])); init = -Inf)
open(outf, "w") do io
    @printf(io, "# %s: %d certified; refined %d (all %d with kappa_cert = 2.5 + the %d largest below); phase A %.0f s, total %.0f s\n",
            logf, length(certs), length(sel), length(top), length(sel) - length(top), tA, time() - t0)
    @printf(io, "# largest refined kappa* = %.10f  -> smallest certified margin among refined = %.4e\n", kmax, 2.5 - kmax)
    @printf(io, "# largest unrefined logged kappa_cert = %.6f (margin %.4e)\n", unref, 2.5 - unref)
    @printf(io, "# SMALLEST CERTIFIED MARGIN (5/2 - kappa*) OF THE LOG = %.4e; rebuild/upper-end failures: %d\n",
            2.5 - max(kmax, unref), nbad)
    sA = sort([x[3] for x in rows if !isnan(x[3])], rev = true)
    @printf(io, "# refined kappa* quantiles (largest first): %s\n", join([@sprintf("%.4f", sA[clamp(round(Int, q * length(sA)), 1, length(sA))]) for q in (0.0, 0.01, 0.1, 0.5, 0.9, 1.0)], " "))
    println(io, "path  logged_kappa_cert  refined_kappa*  margin  bound/mode  phase")
    for x in rows
        @printf(io, "%s  %.6f  %.10f  %.4e  %s  %s\n", x[1], x[2], x[3], 2.5 - x[3], x[4], x[5])
    end
end
@printf(stderr, "smallest certified margin %.4e (refined max kappa* %.10f, unrefined max %.6f), failures %d, %.0f s\n",
        2.5 - max(kmax, unref), kmax, unref, nbad, time() - t0)
