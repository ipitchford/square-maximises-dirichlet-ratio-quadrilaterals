# replay_t.jl — independent re-verification of a transport cover log (cover_t.jl), from the log alone.
#
# usage: julia -t N replay_t.jl <logfile> [--results file] [--certs-sample K] [--max-leaves K] [--region a0,a1,...,d1]
#
# Fail closed: the log is loaded by the strict parser with per-status schemas (logcheck.jl); unknown/FAIL/ERROR/pending
# statuses, malformed lines, duplicate paths, uncovered nodes, unreached records -> nonzero exit.  The unqualified
# "REPLAY OK" (and the receipt <results>_RECEIPT.json) only for an exhaustive run with every leaf verified.
#
# Resumable and record-bound: every result line carries the SHA-256 of the exact log line of its leaf, of the log, and
# of the verifier (all sources + config + Julia/Arblib versions); a cached result is reused only if all three hashes and
# the kind match.  Coverage is always redone.
#
# 1. Coverage: roots recomputed (a-priori box split at S = Float64(1/5)); tree rebuilt from the logged splits
#    (children recomputed exactly); every node must be split or a logged leaf; leaf bounds must match.
# 2. SHELL leaf: must lie inside [0, 0.225]³ × [−0.225, 0.225] (exact rationals) — that region is certified by
#    runs/shell0225_cert.jsonl (replayed separately by replay.jl) together with the local theorem (|p|∞ ≤ 1/20).
# 3. E1–E4 / E3T leaves: independent checks for one fixed condition index: centred forms of the exact polynomials
#    with adaptive bisection (xcheck_cf.jl), then ball arithmetic on a uniform subdivision; only if both are
#    inconclusive, the cover's own Taylor-model test (counted as TM-only).
# 4. COMPOUND leaves: the logged tree is re-expanded (widest-side bisection at the Float64 midpoint) and every leaf
#    is re-checked with its condition by the independent check.
# 5. Certified leaves (cert = transport1): recomputed from scratch: Float64 choice (tsolve1 at the box centre with
#    the logged degree and diagonal), Arb reference integrals, Arb Taylor models, the logged lower-bound candidate
#    rebuilt in Arb (box constants / bootstrap / named Taylor-model bound), Arb criterion in the logged mode.
include("cover_t.jl")
using Random, SHA, Dates
include("xcheck_cf.jl")      # independent centred-form exclusion checks
include("logcheck.jl")       # fail-closed log loading, record and verifier hashes

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

# ---- independent exclusion check (ball arithmetic on a subdivision; no Taylor models)
function vert_balls(lo, hi)
    [let c = Arb(V0[j, k]; prec = PREC)
         for i = 1:4; c += Arb(DV[j, k, i]; prec = PREC) * Arb((lo[i], hi[i]); prec = PREC); end
         c
     end for j = 1:4, k = 1:2]
end
crossA(ax, ay, bx, by) = ax * by - ay * bx
cA(V, i) = (j = mod1(i + 1, 4); k = mod1(i + 2, 4); crossA(V[j, 1] - V[i, 1], V[j, 2] - V[i, 2], V[k, 1] - V[i, 1], V[k, 2] - V[i, 2]))
function subboxes(lo, hi, k)
    out = Tuple{Vector{Float64},Vector{Float64}}[]
    for idx in Iterators.product(ntuple(_ -> 0:k-1, 4)...)
        l = [lo[i] + (hi[i] - lo[i]) * idx[i] / k for i = 1:4]
        h = [idx[i] == k - 1 ? hi[i] : lo[i] + (hi[i] - lo[i]) * (idx[i] + 1) / k for i = 1:4]
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
        Vs = [vert_balls(l, h) for (l, h) in subboxes(lo, hi, k)]
        if kind == "E4"
            for i = 1:4; all(V -> Arblib.ubound(cA(V, i)) <= 0, Vs) && return true; end
        elseif kind == "E2"
            for i = 1:4; all(V -> Arblib.ubound(cA(V, i) - eta * cA(V, mod1(i + 2, 4))) <= 0, Vs) && return true; end
        elseif kind == "E3" || kind == "E3T"
            s02 = kind == "E3" ? s02 : Arb(numerator(RAT_S1^2); prec = PREC) / Arb(denominator(RAT_S1^2); prec = PREC)
            for i = 1:4, kk = 1:4, ll = kk+1:4
                ok = all(Vs) do V
                    ip = mod1(i + 1, 4); ex = V[ip, 1] - V[i, 1]; ey = V[ip, 2] - V[i, 2]
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
# compound tree "dim(left,right)" | "E2" ...
const TMONLY = Threads.Atomic{Int}(0)
"""independent checks (centred forms, then a uniform ball subdivision); if both are inconclusive, the rigorous
Taylor-model exclusion test of the cover itself (counted as TM-only).  Returns (ok, number of TM-only fallbacks)."""
function check_E_or_tm2(kind, lo, hi)
    check_E_cf(kind, lo, hi) && return (true, 0)            # centred forms of the exact polynomials (xcheck_cf.jl)
    check_E(kind, lo, hi; maxk = 4) && return (true, 0)     # ball arithmetic on a uniform subdivision
    if String(exclusion_t(lo, hi)) == kind
        Threads.atomic_add!(TMONLY, 1); return (true, 1)
    end
    return (false, 0)
end
check_E_or_tm(kind, lo, hi) = check_E_or_tm2(kind, lo, hi)[1]
"compound tree check: (ok, TM-only fallbacks, number of sub-leaves)"
function check_compound2(tree::AbstractString, lo, hi)
    tree = strip(tree)
    if startswith(tree, "E")
        ok, ntm = check_E_or_tm2(String(tree), lo, hi)
        return (ok, ntm, 1)
    end
    dimc = tree[1]; dim = findfirst(==(dimc), DIMNAME)
    inner = tree[3:end-1]
    depth = 0; cut = 0
    for (i, ch) in enumerate(inner)
        ch == '(' && (depth += 1); ch == ')' && (depth -= 1)
        if ch == ',' && depth == 0; cut = i; break; end
    end
    dim == argmax(hi .- lo) || return (false, 0, 0)       # the logged split must be the deterministic one
    mid = (lo[dim] + hi[dim]) / 2
    (lo[dim] < mid < hi[dim]) || return (false, 0, 0)
    hi1 = copy(hi); hi1[dim] = mid; lo2 = copy(lo); lo2[dim] = mid
    a = check_compound2(inner[1:cut-1], lo, hi1); a[1] || return (false, a[2], a[3])
    b = check_compound2(inner[cut+1:end], lo2, hi)
    return (b[1], a[2] + b[2], a[3] + b[3])
end
check_compound(tree::AbstractString, lo, hi) = check_compound2(tree, lo, hi)[1]

# ---- certificate recomputation
_rfail(msg) = (false, msg, :LG, nothing, nothing)
"""
    rebuild_cert_t(rec, W) -> (ok, msg, mode, betaA, L0A)

Rebuild a logged transport certificate from the log record alone: Float64 choice at the box centre (logged degree,
diagonal, first-order flag), Arb reference integrals and base Taylor models (left in `W.tpA`), and the logged
lower-bound candidate in Arb.  The criterion itself is evaluated by the caller (any κ).
"""
function rebuild_cert_t(rec, W::TWorker)
    lo = rec["lo"]; hi = rec["hi"]; n = Int(rec["n"]); sh = Int(rec["shift"]); mode = Symbol(rec["mode"]); nm = rec["bound"]
    p0, rad = box_center_radius(lo, hi)
    sb = scalar_bounds(p0, rad; ws = W.wsA)
    L0 = max(sb.L_rect, sb.L_hersch)
    beta0 = max(sb.beta_rect_ks, addlo(L0, sb.gapAC))
    sol = tsolve1(p0, n; shift = sh)
    get(rec, "fo", true) || (sol = merge(sol, (; dC = zero(sol.dC), ddpsi = zero(sol.ddpsi))))
    RA = transport1_refs(Arb, sol)
    b = transport1_base!(W.tpA, RA, p0, rad, W.wsA, sh)
    b.ok || return _rfail("fold")
    betaA = beta0; L0A = L0
    if nm == "const"
    elseif nm == "boot"
        r0 = crit_eval!(W.cpA, W.tpA.M, W.tpA.K, W.tpA.A2, W.wsA; mode = :LG, beta = beta0, L0)
        r0.g1ok || return _rfail("bootstrap g")
        g1up = tm_upper(W.cpA.g1)
        bb = Arb(r0.g2lo; prec = PREC) / Arb(g1up; prec = PREC) + Arb(sb.gapAC; prec = PREC)
        bmax = Float64(Arblib.lbound(bb), RoundDown)
        rec["beta"] <= bmax || return _rfail("beta not justified")
        betaA = max(beta0, rec["beta"])
    elseif startswith(nm, "bootTM")
        btA = bound_tms(Arb, p0, rad, W.wsA)
        lev = startswith(nm, "bootTM2:") ? 2 : 1
        nm1 = split(nm, ":"; limit = 2)[2]
        B1 = nm1 == "const" ? beta0 : (i = findfirst(x -> x[1] == nm1, btA.betas); i === nothing ? nothing : btA.betas[i][2])
        B1 === nothing && return _rfail("bound " * nm)
        base = (; M = W.tpA.M, K = W.tpA.K, A2 = W.tpA.A2)
        b = boot_tm(W.cpA, base, B1, btA.gap, W.wsA); b === nothing && return _rfail("bootTM")
        lev == 2 && (b = boot_tm(W.cpA, base, b, btA.gap, W.wsA); b === nothing && return _rfail("bootTM2"))
        betaA = b
    else
        btA = bound_tms(Arb, p0, rad, W.wsA)
        if mode === :LG
            i = findfirst(x -> x[1] == nm, btA.betas); i === nothing && return _rfail("bound " * nm)
            betaA = btA.betas[i][2]
        else
            i = findfirst(x -> x[1] == nm, btA.L0s); i === nothing && return _rfail("bound " * nm)
            L0A = btA.L0s[i][2]
        end
    end
    return (true, "", mode, betaA, L0A)
end

"""
    cert_with_kstar(rec, W) -> (ok, msg, kstar)

check_cert_t at κ = 5/2; for a record whose logged `kappa_cert` is the ladder fallback 2.5, additionally a certified
κ* by 14 bisection steps in [2, 5/2] (exact rationals; the criterion holds at every accepted value), else NaN.
"""
function cert_with_kstar(rec, W::TWorker)
    ok, msg, mode, betaA, L0A = rebuild_cert_t(rec, W)
    ok || return (false, msg, NaN)
    crit(k) = crit_eval!(W.cpA, W.tpA.M, W.tpA.K, W.tpA.A2, W.wsA; mode, beta = betaA, L0 = L0A, kappa = k)
    res = crit(5 // 2)
    res.ok || return (false, "criterion", NaN)
    which = String(res.which); kst = NaN
    if get(rec, "kappa_cert", 2.5) >= 2.5
        lo = Rational{BigInt}(2); hi = Rational{BigInt}(5, 2)
        for _ = 1:14
            mid = (lo + hi) / 2
            crit(mid).ok ? (hi = mid) : (lo = mid)
        end
        if hi == Rational{BigInt}(5, 2)                    # margin below the grid: keep bisecting towards 5/2
            for _ = 1:46
                mid = (lo + hi) / 2
                crit(mid).ok ? (hi = mid) : (lo = mid)
            end
        end
        kst = Float64(BigFloat(hi; precision = 256), RoundUp)
    end
    return (true, which, kst)
end

function check_cert_t(rec, W::TWorker; kappa = 5 // 2)
    ok, msg, mode, betaA, L0A = rebuild_cert_t(rec, W)
    ok || return (false, msg)
    res = crit_eval!(W.cpA, W.tpA.M, W.tpA.K, W.tpA.A2, W.wsA; mode, beta = betaA, L0 = L0A, kappa)
    res.ok || return (false, "criterion")
    return (true, String(res.which))
end

const VERIFIER_FILES_T = ["replay_t.jl", "cover_t.jl", "cover.jl", "geom.jl", "rules.jl", "fsolve.jl", "tm.jl", "assemble.jl",
                          "criterion.jl", "transport.jl", "bounds_tm.jl", "transport1.jl", "xcheck_cf.jl", "strictjson.jl",
                          "logcheck.jl", "config.toml"]

"read a results file; keep only complete lines bound to this log, this verifier and the current record (hash + kind)"
function load_bound_results(resf, recs, rawsha, src_sha, log_sha)
    done = Dict{String,Dict{String,Any}}(); rej = Dict{String,Int}()
    isfile(resf) || return done, rej
    for ln in eachline(resf)
        bump(k) = (rej[k] = get(rej, k, 0) + 1)
        endswith(ln, "\"end\":1}") || (bump("cut line"); continue)
        q = try
            strict_json_object(ln)
        catch
            bump("unparsable"); continue
        end
        path = get(q, "path", nothing); r = path isa String ? get(recs, path, nothing) : nothing
        r === nothing && (bump("unknown path"); continue)
        get(q, "rec_sha", "") == rawsha[path] || (bump("record hash"); continue)
        get(q, "src_sha", "") == src_sha || (bump("verifier hash"); continue)
        get(q, "log_sha", "") == log_sha || (bump("log hash"); continue)
        get(q, "kind", "") == r["status"] || (bump("kind"); continue)
        (get(q, "ok", nothing) isa Bool) || (bump("schema"); continue)
        done[path] = Dict{String,Any}(k => _tonative(v) for (k, v) in q)    # a later line for the same path wins
    end
    return done, rej
end

"write <results>_RECEIPT.json for an exhaustive, failure-free replay; returns (receipt path, results sha256)"
function write_receipt(resf, logf, log_sha, src_sha, src_files, ver_extra, recs, leaves, cnt, res, stat, ncert,
                       shells, invf, ninv, tot_inv, minmargin, minladder, worstk, nofine)
    res_sha = file_sha(resf)
    rcpt = replace(resf, r"\.jsonl$" => "") * "_RECEIPT.json"
    perkind = join([@sprintf("\"%s\":{\"checked\":%d,\"failed\":%d,\"tm_only_leaves\":%d,\"sub_leaves\":%d,\"thread_s\":%.1f}",
                             k, s[1], s[2], s[3], s[6], s[5]) for (k, s) in sort(collect(stat))], ",")
    files = join([@sprintf("\"%s\":\"%s\"", k, v) for (k, v) in sort(collect(src_files))], ",")
    nex = sum(Int(s[1]) for (k, s) in stat if k != "certified"; init = 0)
    open(rcpt, "w") do f
        @printf(f, "{\"verdict\":\"REPLAY OK\",\"scope\":\"exhaustive: all %d non-SHELL leaves (%d certified, %d exclusion/compound) re-verified under one verifier identity; %d SHELL leaf by exact containment\",\n",
                ncert + nex, ncert, nex, length(shells))
        @printf(f, " \"finished\":\"%s\",\"log\":\"%s\",\"log_sha256\":\"%s\",\"log_records\":%d,\"leaves\":%d,\n",
                Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), logf, log_sha, length(recs), length(leaves))
        @printf(f, " \"leaves_by_status\":{%s},\n", join([@sprintf("\"%s\":%d", k, v) for (k, v) in sort(collect(cnt))], ","))
        @printf(f, " \"results\":\"%s\",\"results_sha256\":\"%s\",\"results_bound_records\":%d,\"invocations\":\"%s\",\"invocations_under_identity\":%d,\"wall_s_total\":%.1f,\n",
                resf, res_sha, length(res), invf, ninv, tot_inv)
        @printf(f, " \"verifier_sha256\":\"%s\",\"verifier_extra\":\"%s\",\"verifier_files\":{%s},\n", src_sha, join(ver_extra, "; "), files)
        @printf(f, " \"per_kind\":{%s},\n", perkind)
        @printf(f, " \"smallest_certified_margin\":%.6e,\"largest_bisected_kappa_star\":%.15f,\"at\":\"%s\",\"ladder_min_margin\":%.6e,\"unresolved_kstar\":%d}\n",
                min(minmargin, minladder), worstk[1], worstk[2], minladder, nofine)
    end
    return rcpt, res_sha
end

function replay_main()
    isempty(ARGS) && (println("usage: julia -t N replay_t.jl <log> [--results f] [--certs-sample K] [--max-leaves K] [--region ...]"); exit(2))
    logf = ARGS[1]
    exhaustive = true; scope = String[]
    roots = root_boxes_t()
    if (i = findfirst(==("--region"), ARGS)) !== nothing
        v = parse.(Float64, split(ARGS[i+1], ","))
        roots = root_boxes_t([v[1], v[3], v[5], v[7]], [v[2], v[4], v[6], v[8]]); exhaustive = false; push!(scope, "region " * ARGS[i+1])
    end
    certsample = (i = findfirst(==("--certs-sample"), ARGS)) === nothing ? typemax(Int) : parse(Int, ARGS[i+1])
    certsample < typemax(Int) && (exhaustive = false; push!(scope, "certificate sample $certsample"))
    maxleaves = (i = findfirst(==("--max-leaves"), ARGS)) === nothing ? typemax(Int) : parse(Int, ARGS[i+1])
    maxleaves < typemax(Int) && (exhaustive = false; push!(scope, "at most $maxleaves leaves this invocation"))
    resf = (i = findfirst(==("--results"), ARGS)) === nothing ? replace(logf, r"\.jsonl$" => "") * "_replay_results.jsonl" : ARGS[i+1]
    invf = replace(resf, r"\.jsonl$" => "") * "_invocations.jsonl"
    twall = time(); tstarted = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    # ---- identities: the log and the verifier (sources, config, Julia/Arblib versions, precision)
    log_sha = file_sha(logf)
    vfiles = [joinpath(@__DIR__, f) for f in VERIFIER_FILES_T]
    ver_extra = ["julia " * string(VERSION), "Arblib " * string(pkgversion(Arblib)), "PREC " * string(PREC)]
    src_sha, src_files = verifier_identity(vfiles; extra = ver_extra)
    @printf("log %s sha256 %s\nverifier sha256 %s (%s; %s)\n", logf, log_sha, src_sha, join(ver_extra, ", "), join(VERIFIER_FILES_T, " "))
    flush(stdout)
    # ---- fail-closed loading: strict JSON, status whitelist, per-status schema, duplicates
    recs, rawsha, problems = load_log_strict(logf, schema_transport)
    if !isempty(problems)
        for p in problems[1:min(end, 50)]; println("  LOG PROBLEM ", p); end
        @printf("%d log problems\nREPLAY FOUND PROBLEMS (log rejected by the strict loader)\n", length(problems)); flush(stdout)
        exit(2)
    end
    # ---- coverage: tree rebuilt from the roots; every record must be reached exactly once
    leaves = Dict{String,Any}[]; uncovered = String[]; mismatch = String[]; reached = 0
    stack = copy(roots)
    while !isempty(stack)
        path, lo, hi = pop!(stack)
        r = get(recs, path, nothing)
        if r === nothing; push!(uncovered, path * " (no record: pending/uncovered)"); continue; end
        reached += 1
        if r["status"] == "split"
            c1, c2 = children(path, lo, hi, Int(r["dim"])); push!(stack, c1, c2)
        else
            (r["lo"] == lo && r["hi"] == hi) || push!(mismatch, path)
            push!(leaves, r)
        end
    end
    unreached = length(recs) - reached
    cnt = Dict{String,Int}(); for r in leaves; cnt[r["status"]] = get(cnt, r["status"], 0) + 1; end
    @printf("coverage: %d records, %d leaves %s, %d uncovered, %d bound mismatches, %d unreached records\n",
            length(recs), length(leaves), string(cnt), length(uncovered), length(mismatch), unreached); flush(stdout)
    if !isempty(uncovered) || !isempty(mismatch) || (unreached != 0 && exhaustive)
        for u in uncovered[1:min(end, 30)]; println("  UNCOVERED ", u); end
        for m in mismatch[1:min(end, 30)]; println("  MISMATCH ", m); end
        println("REPLAY FOUND PROBLEMS (coverage)"); flush(stdout); exit(3)
    end
    # ---- hash-bound cached results (resumable)
    done, rej = load_bound_results(resf, recs, rawsha, src_sha, log_sha)
    @printf("results file %s: %d bound results accepted, rejected %s\n", resf, length(done), string(rej)); flush(stdout)
    bad = String[]; lk = ReentrantLock()
    shells = filter(r -> r["status"] == "SHELL", leaves)
    for r in shells
        is_shell_box(r["lo"], r["hi"]) || push!(bad, r["path"] * " SHELL outside the certified shell region")
    end
    ex = filter(r -> r["status"] in ("E1", "E2", "E3", "E3T", "E4", "COMPOUND"), leaves)
    certs = filter(r -> r["status"] == "certified", leaves)
    length(certs) > certsample && (certs = certs[sort(randperm(MersenneTwister(1), length(certs))[1:certsample])])
    todo = [r for r in vcat(ex, certs) if !haskey(done, r["path"])]
    todo = todo[randperm(MersenneTwister(2), length(todo))]
    length(todo) > maxleaves && (todo = todo[1:maxleaves])
    @printf("to check now: %d (%d exclusion/compound, %d certificates)\n", length(todo),
            count(r -> r["status"] != "certified", todo), count(r -> r["status"] == "certified", todo)); flush(stdout)
    san(s) = replace(String(s), '"' => '\'', '\\' => '/', '\n' => ' ', '\r' => ' ')
    if isfile(resf) && filesize(resf) > 0                    # terminate a line cut by a crash
        lastb = open(f -> (seekend(f); skip(f, -1); read(f, UInt8)), resf)
        lastb == UInt8('\n') || open(f -> write(f, "\n"), resf, "a")
    end
    io = open(resf, "a")
    ndone = Threads.Atomic{Int}(0); nfail = Threads.Atomic{Int}(0); t0 = time(); lastrep = Ref(time())
    chunks = collect(Iterators.partition(todo, 32))
    Threads.@threads :dynamic for ch in chunks
        W = TWorker()
        for r in ch
            t1 = time(); ntm = 0; nsub = 1; kst = NaN
            if r["status"] == "certified"
                ok, why, kst = try
                    cert_with_kstar(r, W)
                catch err
                    (false, "error: " * first(sprint(showerror, err), 120), NaN)
                end
            else
                ok, ntm, nsub = try
                    r["status"] == "COMPOUND" ? check_compound2(r["tree"], r["lo"], r["hi"]) : (check_E_or_tm2(r["status"], r["lo"], r["hi"])..., 1)
                catch err
                    (false, 0, 0)
                end
                why = ok ? "" : "exclusion not re-verified"
            end
            line = @sprintf("{\"path\":\"%s\",\"kind\":\"%s\",\"ok\":%s,\"why\":\"%s\",\"t\":%.4f,\"tm\":%d,\"nsub\":%d,\"kstar\":%s,\"rec_sha\":\"%s\",\"src_sha\":\"%s\",\"log_sha\":\"%s\",\"end\":1}",
                            r["path"], r["status"], ok ? "true" : "false", san(why), time() - t1, ntm, nsub,
                            isnan(kst) ? "null" : repr(kst), rawsha[r["path"]], src_sha, log_sha)
            lock(lk) do
                println(io, line); flush(io)
                if time() - lastrep[] > 60
                    lastrep[] = time(); el = time() - t0; nd = ndone[] + 1
                    @printf("[%s] checked %d / %d, failures %d, elapsed %.0fs, ETA %.0fs\n", Dates.format(now(), "HH:MM:SS"),
                            nd, length(todo), nfail[], el, el / nd * (length(todo) - nd)); flush(stdout)
                end
            end
            Threads.atomic_add!(ndone, 1); ok || Threads.atomic_add!(nfail, 1)
        end
    end
    close(io)
    tinv = time() - twall
    open(invf, "a") do f
        @printf(f, "{\"start\":\"%s\",\"end\":\"%s\",\"wall\":%.1f,\"checked\":%d,\"threads\":%d,\"src_sha\":\"%s\",\"log_sha\":\"%s\"}\n",
                tstarted, Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), tinv, ndone[], Threads.nthreads(), src_sha, log_sha)
    end
    @printf("this invocation: checked %d leaves in %.0fs wall (%d threads)\n", ndone[], time() - t0, Threads.nthreads()); flush(stdout)
    # ---- aggregate over all bound results (this invocation and earlier ones under the same identities)
    res, rej2 = load_bound_results(resf, recs, rawsha, src_sha, log_sha)
    stat = Dict{String,Vector{Float64}}()                         # kind => [n, failed, tm-only leaves, tm sub-checks, sum t, sub-leaves]
    notchecked = 0; minmargin = Inf; worstk = Any[-Inf, ""]; nofine = 0; minladder = Inf
    for r in vcat(ex, certs)
        q = get(res, r["path"], nothing)
        if q === nothing                                   # missing results are a problem only for an exhaustive run
            notchecked += 1; exhaustive && push!(bad, r["path"] * " " * r["status"] * " not checked"); continue
        end
        s = get!(stat, r["status"], zeros(6))
        s[1] += 1; s[5] += q["t"]; s[6] += q["nsub"]
        q["tm"] > 0 && (s[3] += 1; s[4] += q["tm"])
        if !q["ok"]
            s[2] += 1; push!(bad, r["path"] * " " * r["status"] * ": " * q["why"])
        elseif r["status"] == "certified"
            kc = r["kappa_cert"]; ks = q["kstar"]
            if kc < 2.5
                minladder = min(minladder, 2.5 - kc)
            elseif ks >= 0
                minmargin = min(minmargin, 2.5 - ks)
                ks > worstk[1] && (worstk[1] = ks; worstk[2] = r["path"])
            else
                nofine += 1
            end
        end
    end
    for k in sort(collect(keys(stat)))
        s = stat[k]
        @printf("  %-10s checked %6d  failed %d  TM-only leaves %d (TM-only sub-checks %d)  sub-leaves %d  check time %.0f thread-s (%.3f s each)\n",
                k, s[1], s[2], s[3], s[4], s[6], s[5], s[5] / max(s[1], 1))
    end
    ncert = haskey(stat, "certified") ? Int(stat["certified"][1]) : 0
    nlogc = count(r -> r["status"] == "certified", leaves)
    @printf("certificates re-verified: %d of %d; exclusion/compound leaves re-checked: %d of %d; SHELL leaves %d (exact containment)\n",
            ncert, nlogc, sum(Int(s[1]) for (k, s) in stat if k != "certified"; init = 0), length(ex), length(shells))
    @printf("smallest certified margin 5/2 − κ*: %.4e (replay bisection over the %d ladder-fallback certificates; largest κ* %.12f at %s); ladder-logged certificates: margin ≥ %.4e; unresolved %d\n",
            min(minmargin, minladder), count(r -> r["status"] == "certified" && r["kappa_cert"] >= 2.5, certs), worstk[1], worstk[2], minladder, nofine)
    for b in bad[1:min(end, 50)]; println("  BAD ", b); end
    nbad = length(bad)
    tot_inv = 0.0; ninv = 0
    if isfile(invf)
        for ln in eachline(invf)
            q = try strict_json_object(ln) catch; continue end
            (get(q, "src_sha", "") == src_sha && get(q, "log_sha", "") == log_sha) || continue
            tot_inv += q["wall"]; ninv += 1
        end
    end
    @printf("wall time: this invocation %.0f s; all %d recorded invocations under this identity %.0f s\n", tinv, ninv, tot_inv)
    if nbad > 0
        println("REPLAY FOUND PROBLEMS"); flush(stdout); exit(1)
    elseif !exhaustive
        println("PARTIAL REPLAY PASSED (not exhaustive: ", join(scope, "; "), ")"); flush(stdout); exit(0)
    end
    # ---- exhaustive and clean: write the receipt
    rcpt, res_sha = write_receipt(resf, logf, log_sha, src_sha, src_files, ver_extra, recs, leaves, cnt, res, stat, ncert,
                                  shells, invf, ninv, tot_inv, minmargin, minladder, worstk, nofine)
    println("receipt written: ", rcpt, " (results sha256 ", res_sha, ")")
    println("REPLAY OK"); flush(stdout); exit(0)
end

if abspath(PROGRAM_FILE) == @__FILE__
    replay_main()
end
