# replay.jl — independent re-verification of a cover log (JSONL written by cover.jl), from the log alone.
#
# usage: julia -t 8 replay.jl <logfile> [--certs-sample K] [--region a0,a1,b0,b1,c0,c1,d0,d1]
#
# Fail closed (logcheck.jl): strict JSON, status whitelist {split, certified, E1–E4} with per-status schemas, duplicate
# paths, pending/uncovered nodes and unreached records give a nonzero exit; "REPLAY OK" only for exhaustive runs
# (a --certs-sample run ends with "PARTIAL REPLAY PASSED").
#
# 1. Coverage: the root boxes are recomputed from config.toml and the a-priori box; the tree is rebuilt from the
#    logged split records (each split halves one coordinate at the Float64 midpoint, so children are recomputed
#    exactly); every node must be either split or a logged leaf, and every leaf's logged bounds must equal the
#    recomputed ones.  Any node without a record, and any FAIL/ERROR leaf, is reported as UNCOVERED.
# 2. Exclusions: each E1–E4 leaf is re-checked with an implementation independent of cover.jl's Taylor-model
#    tests: exact rational comparison (E1) and natural ball-arithmetic extension of the polynomial conditions on a
#    subdivision of the box (E2–E4; one fixed condition index for the whole box, as required).
# 3. Certificates: each certified leaf is recomputed from scratch (Float64 choice of the test functions at the box
#    centre with the logged degree n, rigorous Arb Taylor models, rigorous β/L0, Arb criterion with the logged
#    mode).  The logged β is accepted only if it is ≤ the rigorously recomputed admissible shift.
include("cover.jl")
using Random
include("logcheck.jl")       # fail-closed log loading (strict JSON, status whitelist, schemas, duplicates)

# ---------------------------------------------------------------- log parsing (flat JSON subset)
function parse_rec(ln)
    r = Dict{String,Any}()
    for m in eachmatch(r"\"([A-Za-z0-9_]+)\":(\"[^\"]*\"|\[[^\]]*\]|[-+0-9.eE]+|true|false|null)", ln)
        k = m.captures[1]; v = m.captures[2]
        if startswith(v, "\"")
            r[k] = v[2:end-1]
        elseif startswith(v, "[")
            r[k] = parse.(Float64, split(v[2:end-1], ","))
        elseif v == "true" || v == "false"
            r[k] = v == "true"
        elseif v == "null"
            r[k] = NaN
        else
            r[k] = parse(Float64, v)
        end
    end
    r
end

# ---------------------------------------------------------------- independent exclusion checks
# vertex coordinate enclosures over a box (affine in p: exact ranges)
function vert_balls(lo, hi)
    [let c = Arb(V0[j, k]; prec = PREC)
         for i = 1:4
             c += Arb(DV[j, k, i]; prec = PREC) * Arb((lo[i], hi[i]); prec = PREC)
         end
         c
     end for j = 1:4, k = 1:2]
end
crossA(ax, ay, bx, by) = ax * by - ay * bx
cA(V, i) = (j = mod1(i + 1, 4); k = mod1(i + 2, 4);
            crossA(V[j, 1] - V[i, 1], V[j, 2] - V[i, 2], V[k, 1] - V[i, 1], V[k, 2] - V[i, 2]))
function subboxes(lo, hi, k)
    out = Tuple{Vector{Float64},Vector{Float64}}[]
    for idx in Iterators.product(ntuple(_ -> 0:k-1, 4)...)
        l = [lo[i] + (hi[i] - lo[i]) * idx[i] / k for i = 1:4]
        h = [idx[i] == k - 1 ? hi[i] : lo[i] + (hi[i] - lo[i]) * (idx[i] + 1) / k for i = 1:4]
        # make the sub-boxes cover [lo,hi]: round outward
        l = [idx[i] == 0 ? lo[i] : prevfloat(l[i]) for i = 1:4]
        h = [idx[i] == k - 1 ? hi[i] : nextfloat(h[i]) for i = 1:4]
        push!(out, (l, h))
    end
    out
end
function check_E(kind, lo, hi; maxk = 8)
    if kind == "E1"
        return all(abs(Rational{BigInt}(lo[i])) <= CUBE && abs(Rational{BigInt}(hi[i])) <= CUBE for i = 1:4)
    end
    eta = Arb(numerator(RAT_ETA); prec = PREC) / Arb(denominator(RAT_ETA); prec = PREC)
    s02 = Arb(numerator(RAT_S0^2); prec = PREC) / Arb(denominator(RAT_S0^2); prec = PREC)
    k = 1
    while k <= maxk
        sbs = subboxes(lo, hi, k)
        Vs = [vert_balls(l, h) for (l, h) in sbs]
        if kind == "E4"
            for i = 1:4
                all(V -> Arblib.ubound(cA(V, i)) <= 0, Vs) && return true
            end
        elseif kind == "E2"
            for i = 1:4
                all(V -> Arblib.ubound(cA(V, i) - eta * cA(V, mod1(i + 2, 4))) <= 0, Vs) && return true
            end
        elseif kind == "E3"
            for i = 1:4, kk = 1:4, ll = kk+1:4
                ok = all(Vs) do V
                    ip = mod1(i + 1, 4)
                    ex = V[ip, 1] - V[i, 1]; ey = V[ip, 2] - V[i, 2]
                    e2 = ex^2 + ey^2; d2 = (V[kk, 1] - V[ll, 1])^2 + (V[kk, 2] - V[ll, 2])^2
                    all(j -> Arblib.ubound(crossA(ex, ey, V[j, 1] - V[i, 1], V[j, 2] - V[i, 2])^2 - s02 * e2 * d2) <= 0, (mod1(i + 2, 4), mod1(i + 3, 4)))
                end
                ok && return true
            end
        end
        k *= 2
    end
    return false
end

# ---------------------------------------------------------------- certificate recomputation
function check_cert(rec, W::Worker)
    lo = rec["lo"]; hi = rec["hi"]; n = Int(rec["n"]); mode = Symbol(rec["mode"]); beta_log = rec["beta"]
    p0, rad = box_center_radius(lo, hi)
    sb = scalar_bounds(p0, rad; ws = W.wsA)
    L0 = max(sb.L_rect, sb.L_hersch)
    beta0 = max(sb.beta_rect_ks, addlo(L0, sb.gapAC))
    sol = fsolve_d(p0, n)
    baseA = assemble_box!(W.awA, p0, rad, sol)
    baseA.ok || return (false, "assembly", NaN)
    beta = beta0
    if mode === :LG && beta_log > beta0
        r0 = crit_eval!(W.cpA, baseA.M, baseA.K, baseA.A2, W.wsA; mode = :LG, beta = beta0, L0)
        r0.g1ok || return (false, "bootstrap g", NaN)
        g1up = tm_upper(W.cpA.g1)
        bA = Arb(r0.g2lo; prec = PREC) / Arb(g1up; prec = PREC) + Arb(sb.gapAC; prec = PREC)
        bmax = Float64(Arblib.lbound(bA), RoundDown)
        beta_log <= bmax || return (false, "beta not justified", NaN)
        beta = beta_log
    end
    res = crit_eval!(W.cpA, baseA.M, baseA.K, baseA.A2, W.wsA; mode, beta, L0)
    res.ok || return (false, "criterion", NaN)
    return (true, String(res.which), get(rec, "kappa_cert", NaN))
end

function replay_main_s()
    isempty(ARGS) && (println("usage: julia -t N replay.jl <log> [--certs-sample K] [--region a0,a1,b0,b1,c0,c1,d0,d1]"); exit(2))
    logf = ARGS[1]
    roots = root_boxes(); region = ""
    if (i = findfirst(==("--region"), ARGS)) !== nothing
        v = parse.(Float64, split(ARGS[i+1], ","))
        roots = root_boxes([v[1], v[3], v[5], v[7]], [v[2], v[4], v[6], v[8]]); region = ARGS[i+1]
    end
    certsample = (i = findfirst(==("--certs-sample"), ARGS)) === nothing ? typemax(Int) : parse(Int, ARGS[i+1])
    @printf("log %s sha256 %s\n", logf, file_sha(logf))
    # fail closed: strict JSON, status whitelist (split, certified, E1–E4), per-status schema, duplicate paths
    recs, rawsha, problems = load_log_strict(logf, schema_shell)
    if !isempty(problems)
        for p in problems[1:min(end, 50)]; println("  LOG PROBLEM ", p); end
        @printf("%d log problems\nREPLAY FOUND PROBLEMS (log rejected by the strict loader)\n", length(problems)); exit(2)
    end
    # 1. coverage (every record must be reached exactly once; nodes without a record are pending/uncovered)
    leaves = Dict{String,Any}[]; uncovered = String[]; mismatch = String[]; reached = 0
    stack = copy(roots)
    while !isempty(stack)
        path, lo, hi = pop!(stack)
        r = get(recs, path, nothing)
        if r === nothing
            push!(uncovered, path * " (no record)"); continue
        end
        reached += 1
        if r["status"] == "split"
            c1, c2 = children(path, lo, hi, Int(r["dim"])); push!(stack, c1, c2)
        else
            (r["lo"] == lo && r["hi"] == hi) || push!(mismatch, path)
            push!(leaves, r)
        end
    end
    unreached = length(recs) - reached
    @printf("coverage: %d leaves, %d uncovered, %d bound mismatches, %d records unreached\n", length(leaves), length(uncovered),
            length(mismatch), unreached)
    if !isempty(uncovered) || !isempty(mismatch) || unreached != 0
        for u in uncovered[1:min(end, 50)]; println("  UNCOVERED ", u); end
        for m in mismatch[1:min(end, 50)]; println("  MISMATCH ", m); end
        println("REPLAY FOUND PROBLEMS (coverage)"); exit(3)
    end
    # 2./3. verification (threaded); every leaf status is one of certified, E1–E4 (schema)
    certs = filter(r -> r["status"] == "certified", leaves)
    excls = filter(r -> r["status"] in ("E1", "E2", "E3", "E4"), leaves)
    length(certs) + length(excls) == length(leaves) || (println("REPLAY FOUND PROBLEMS (unclassified leaves)"); exit(1))
    ncert_all = length(certs)
    if length(certs) > certsample
        certs = certs[sort(randperm(MersenneTwister(1), length(certs))[1:certsample])]
    end
    bad = Threads.Atomic{Int}(0); badlist = String[]; lk = ReentrantLock()
    Threads.@threads for r in excls
        ok = check_E(r["status"], r["lo"], r["hi"])
        ok || (Threads.atomic_add!(bad, 1); lock(lk) do; push!(badlist, r["path"] * " " * r["status"]); end)
    end
    @printf("exclusions re-checked: %d, failed: %d\n", length(excls), bad[])
    badc = Threads.Atomic{Int}(0); minmargin = Ref(Inf); t0 = time()
    chunks = collect(Iterators.partition(certs, max(1, cld(length(certs), 8 * Threads.nthreads()))))
    Threads.@threads for ch in chunks
        W = Worker()
        for r in ch
            ok, why, kc = try
                check_cert(r, W)
            catch err
                (false, "error: " * first(sprint(showerror, err), 120), NaN)
            end
            lock(lk) do
                if ok
                    minmargin[] = min(minmargin[], 2.5 - kc)
                else
                    Threads.atomic_add!(badc, 1); push!(badlist, r["path"] * " certified: " * why)
                end
            end
        end
    end
    @printf("certificates re-verified: %d of %d, failed: %d, min logged margin 2.5-kappa = %.5f, time %.0fs\n",
            length(certs), ncert_all, badc[], minmargin[], time() - t0)
    for b in badlist[1:min(end, 50)]; println("  BAD ", b); end
    if bad[] != 0 || badc[] != 0
        println("REPLAY FOUND PROBLEMS"); exit(1)
    elseif length(certs) < ncert_all
        @printf("PARTIAL REPLAY PASSED (not exhaustive: %d of %d certificates sampled)\n", length(certs), ncert_all); exit(0)
    else
        println(isempty(region) ? "REPLAY OK" : "REPLAY OK (exhaustive for the region " * region * ")"); exit(0)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    replay_main_s()
end
