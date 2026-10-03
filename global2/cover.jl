# cover.jl — adaptive bisection cover of the D4 fundamental domain ∩ a-priori box.
#
# usage:  julia -t 8 cover.jl <logfile.jsonl> [--pilot] [--maxboxes N] [--region a0,a1,b0,b1,c0,c1,d0,d1]
#   --pilot : Float64 only (no Arb): mimics the rigorous test (same Taylor-model formulas in Float64) to predict
#             box counts and time.  Leaves are logged with status "pilot_pass" instead of "certified".
#   without --pilot: every leaf is either discarded by a rigorous exclusion test (E1–E4, Arb) or certified by the
#             Arb certificate; boxes that reach the minimum size are logged as "FAIL" (uncovered).
# The log is append-only JSONL: one line per leaf (certified / E1–E4 / FAIL / ERROR) and one line per split
# ("status":"split","dim":k).  Restarting with the same log resumes: logged splits are replayed, logged leaves
# skipped, everything else (re)processed.  Split decisions are deterministic functions of the box.

include("geom.jl"); include("rules.jl"); include("fsolve.jl"); include("tm.jl"); include("assemble.jl"); include("criterion.jl")
using Printf, LinearAlgebra, Dates
BLAS.set_num_threads(1)
set_order!(3)

const CFG = joinpath(@__DIR__, "config.toml")
const CUBE = config_rat(CFG, "cube")                  # exact rational (e.g. 49//1000)
@assert config_rat(CFG, "eta_neartri") == RAT_ETA && config_rat(CFG, "thin_s0") == RAT_S0
const CUBE_F = Float64(CUBE, RoundDown)               # split point inside the cube
const A_PRIORI = ([0.0, 0.0, 0.0, -Float64(905 // 1000, RoundUp)],
                  [Float64(905 // 1000, RoundUp), Float64(436 // 100, RoundUp), Float64(436 // 100, RoundUp), Float64(905 // 1000, RoundUp)])
const MINW = 2.0^-12        # boxes narrower than this in every coordinate are logged FAIL
const WTRY = 0.3            # do not attempt a certificate while some width exceeds this
const DIMNAME = "abcd"
"rigorous lower bound of a + b for Float64 lower bounds a, b (round-to-nearest addition could round up)"
addlo(a::Float64, b::Float64) = Float64(Arblib.lbound(Arb(a; prec = PREC) + Arb(b; prec = PREC)), RoundDown)

function root_boxes(lo0 = A_PRIORI[1], hi0 = A_PRIORI[2])
    s = CUBE_F
    cut(l, h, pts) = (xs = vcat(l, [x for x in pts if l < x < h], h); [(xs[k], xs[k+1]) for k = 1:length(xs)-1])
    ivs = [cut(lo0[i], hi0[i], [s]) for i = 1:3]
    push!(ivs, cut(lo0[4], hi0[4], [-s, s]))
    roots = Tuple{String,Vector{Float64},Vector{Float64}}[]
    k = 0
    for ia in ivs[1], ib in ivs[2], ic in ivs[3], id in ivs[4]
        k += 1
        push!(roots, ("R$k", [ia[1], ib[1], ic[1], id[1]], [ia[2], ib[2], ic[2], id[2]]))
    end
    roots
end

children(path, lo, hi, dim) = begin
    mid = (lo[dim] + hi[dim]) / 2
    @assert lo[dim] < mid < hi[dim]
    hi1 = copy(hi); hi1[dim] = mid; lo2 = copy(lo); lo2[dim] = mid
    ((path * "." * DIMNAME[dim] * "0", copy(lo), hi1), (path * "." * DIMNAME[dim] * "1", lo2, copy(hi)))
end

# ------------------------------------------------------------------ per-thread workspace
struct Worker
    awF::AsmWork{Float64,Float64}; awA::AsmWork{Arb,Mag}
    cpF::CritPool{Float64,Float64}; cpA::CritPool{Arb,Mag}
    wsF::TMWork{Float64,Float64}; wsA::TMWork{Arb,Mag}
end
Worker() = Worker(newasm(Float64), newasm(Arb), critpool(Float64), critpool(Arb), newwork(Float64), newwork(Arb))

# point criterion (Float64) for choosing n: smallest κ with θ2 < κ L at the point (L = max(LG, L0))
function point_kappa_f(s, beta, L0)
    LG = beta * s.A1 > s.A0 ? lg_bound(s.A0, s.A1, s.A2, beta) : 0.0
    return s.theta[2] / max(LG, L0)
end

function choose_n(p0, beta, L0)
    s8 = fsolve(p0, 8); k8 = point_kappa_f(s8, beta, L0)
    k8 < 2.0 && return 8, k8, k8
    s12 = fsolve(p0, 12); k12 = point_kappa_f(s12, beta, L0)
    (2.5 - k8) >= 0.75 * (2.5 - k12) && return 8, k8, k12
    s10 = fsolve(p0, 10); k10 = point_kappa_f(s10, beta, L0)
    (2.5 - k10) >= 0.75 * (2.5 - k12) && return 10, k10, k12
    return 12, k12, k12
end

# influence of each variable on a TM (sum of |coefficients| of monomials containing t_i)
function influence(f::TM{Float64})
    tab = TAB[]; infl = zeros(4)
    for k = 2:tab.nc, i = 1:4
        tab.mono[k][i] > 0 && (infl[i] += abs(f.c[k]))
    end
    infl
end

# κ* (certified bound on ξ over the box) by bisection on κ (each evaluation is a full certificate)
function kappa_star!(cp, base, ws, mode, beta, L0, kap0::Rational)
    lo = 1.0; hi = Float64(kap0)
    for _ = 1:14
        mid = (lo + hi) / 2
        rr = crit_eval!(cp, base.M, base.K, base.A2, ws; mode, beta, L0, kappa = rationalize(BigInt, mid; tol = 1e-12))
        rr.ok ? (hi = mid) : (lo = mid)
    end
    return hi
end

"""
    process(path, lo, hi, W, pilot) -> (leaf::Union{Nothing,Dict}, split::Int)

Either returns a leaf record (Dict) or a split dimension.
"""
function process(path, lo, hi, W::Worker, pilot::Bool)
    t0 = time()
    rec = Dict{String,Any}("path" => path, "lo" => lo, "hi" => hi)
    ex = exclusion(lo, hi; cube = CUBE, ws = W.wsA)
    if ex !== :none
        rec["status"] = String(ex); rec["t"] = time() - t0
        return rec, 0
    end
    w = hi .- lo
    if maximum(w) < MINW
        rec["status"] = "FAIL"; rec["why"] = "min size"; rec["t"] = time() - t0
        return rec, 0
    end
    if maximum(w) > WTRY
        return nothing, argmax(w)
    end
    p0, rad = box_center_radius(lo, hi)
    cls = classify_point(p0; cube = Float64(CUBE), eta = Float64(RAT_ETA), s0 = Float64(RAT_S0))
    if cls === :nonconvex
        return nothing, argmax(w)
    end
    sb = scalar_bounds(p0, rad; ws = W.wsA)
    L0 = max(sb.L_rect, sb.L_hersch)
    beta = max(sb.beta_rect_ks, addlo(L0, sb.gapAC))
    beta0 = beta
    n, kpt, kbest = choose_n(p0, beta, L0)
    if kpt >= 2.5
        rec["why"] = "point kappa >= 2.5"
        return maximum(w) < 4MINW ? (merge(rec, Dict("status" => "FAIL", "t" => time() - t0, "kpt" => kpt, "n" => n)), 0) : (nothing, argmax(w))
    end
    sol = fsolve_d(p0, n)
    baseF = assemble_box!(W.awF, p0, rad, sol)
    baseF.ok || return nothing, argmax(w)
    # bootstrap β with the LG lower bound itself (λ2 ≥ λ1 + 3π²/D², Andrews–Clutterbuck)
    resF = crit_eval!(W.cpF, baseF.M, baseF.K, baseF.A2, W.wsF; mode = :LG, beta, L0)
    mode = :LG
    if resF.g1ok
        Llo = resF.g2lo / tm_upper(W.cpF.g1)          # Float64 candidate only (recomputed rigorously in Arb)
        if Llo + sb.gapAC > beta
            beta2 = Llo + sb.gapAC
            res2 = crit_eval!(W.cpF, baseF.M, baseF.K, baseF.A2, W.wsF; mode = :LG, beta = beta2, L0)
            if res2.ok || !resF.ok
                beta = beta2; resF = res2
            end
        end
    end
    if !resF.ok
        resM = crit_eval!(W.cpF, baseF.M, baseF.K, baseF.A2, W.wsF; mode = :MONO, beta, L0)
        if resM.ok
            resF = resM; mode = :MONO
        end
    end
    if !resF.ok
        # split: variable with the largest influence on the most promising failing test
        crit_eval!(W.cpF, baseF.M, baseF.K, baseF.A2, W.wsF; mode = :LG, beta, L0)
        f = !resF.g1ok ? W.cpF.g2 : (resF.a11up >= 0 ? W.cpF.At[1] : W.cpF.s11)
        infl = influence(f)
        return nothing, argmax(infl)
    end
    rec["n"] = n; rec["mode"] = String(mode); rec["beta"] = beta; rec["L0"] = L0; rec["kpt"] = kpt
    rec["which_f64"] = String(resF.which)
    rec["tF"] = time() - t0
    if pilot
        rec["status"] = "pilot_pass"
        rec["kappa_f64"] = kappa_star!(W.cpF, baseF, W.wsF, mode, beta, L0, 5 // 2)
        rec["t"] = time() - t0
        return rec, 0
    end
    # ---- rigorous (Arb) certificate
    baseA = assemble_box!(W.awA, p0, rad, sol)
    if !baseA.ok
        return nothing, argmax(w)
    end
    # β used in the Arb certificate: recompute the Andrews–Clutterbuck bootstrap rigorously (the Float64 value
    # chosen above is only a candidate): β_A = lower bound of [ L_LG(β0) + 3π²/D² ], β0 = beta0 (rigorous).
    betaA = beta0
    if mode === :LG && beta > beta0
        r0 = crit_eval!(W.cpA, baseA.M, baseA.K, baseA.A2, W.wsA; mode = :LG, beta = beta0, L0)
        if r0.g1ok
            g1up = tm_upper(W.cpA.g1)
            bA = Arb(r0.g2lo; prec = PREC) / Arb(g1up; prec = PREC) + Arb(sb.gapAC; prec = PREC)
            betaA = max(beta0, Float64(Arblib.lbound(bA), RoundDown))
        end
    end
    beta = min(beta, betaA)
    rec["beta"] = beta
    resA = crit_eval!(W.cpA, baseA.M, baseA.K, baseA.A2, W.wsA; mode, beta, L0)
    if !resA.ok
        rec["arb_fail"] = true
        return nothing, argmax(influence(W.cpF.s11))
    end
    rec["status"] = "certified"; rec["which"] = String(resA.which)
    rec["a11up"] = resA.a11up; rec["g1lo"] = resA.g1lo; rec["g2lo"] = resA.g2lo
    rec["ups"] = Dict(String(k) => v for (k, v) in pairs(resA.ups))
    rec["kappa_cert"] = kappa_star!(W.cpA, baseA, W.wsA, mode, beta, L0, 5 // 2)
    rec["t"] = time() - t0
    return rec, 0
end

# ------------------------------------------------------------------ JSON (flat, numbers/strings/arrays/dicts)
jsonval(x::AbstractString) = "\"" * x * "\""
jsonval(x::Symbol) = jsonval(String(x))
jsonval(x::Bool) = x ? "true" : "false"
jsonval(x::Integer) = string(x)
jsonval(x::AbstractFloat) = isfinite(x) ? repr(Float64(x)) : "null"
jsonval(x::AbstractVector) = "[" * join(jsonval.(x), ",") * "]"
jsonval(x::AbstractDict) = "{" * join([jsonval(string(k)) * ":" * jsonval(v) for (k, v) in x], ",") * "}"

function read_log(logfile)
    recs = Dict{String,String}()          # path => status ("split:<dim>" for splits)
    isfile(logfile) || return recs
    for ln in eachline(logfile)
        mp = match(r"\"path\":\"([^\"]+)\"", ln); ms = match(r"\"status\":\"([^\"]+)\"", ln)
        (mp === nothing || ms === nothing) && continue
        st = ms.captures[1]
        if st == "split"
            md = match(r"\"dim\":(\d)", ln); st = "split:" * md.captures[1]
        end
        recs[mp.captures[1]] = st
    end
    recs
end

# ------------------------------------------------------------------ driver
function run_cover(logfile; pilot = false, maxboxes = typemax(Int), roots = root_boxes(), procfun = process, mkworker = Worker)
    prev = read_log(logfile)
    io = open(logfile, "a")
    lk = ReentrantLock()
    queue = Channel{Tuple{String,Vector{Float64},Vector{Float64}}}(Inf)
    pending = Threads.Atomic{Int}(0)
    nproc = Threads.Atomic{Int}(0)
    stopped = Threads.Atomic{Bool}(false)
    for r in roots
        Threads.atomic_add!(pending, 1); put!(queue, r)
    end
    stats = Dict{String,Int}()
    tstart = time()
    lastreport = Ref(time())
    STOP = ("STOP", Float64[], Float64[])
    function finish_item!()
        if Threads.atomic_sub!(pending, 1) == 1      # was 1, now 0: all work done
            for _ = 1:Threads.nthreads(); put!(queue, STOP); end
        end
    end
    function worker_loop()
        W = mkworker()
        while true
            item = take!(queue)
            item[1] == "STOP" && break
            stopped[] && continue
            path, lo, hi = item
            st = get(prev, path, "")
            if startswith(st, "split:")
                dim = parse(Int, st[7:end])
                c1, c2 = children(path, lo, hi, dim)
                Threads.atomic_add!(pending, 2); put!(queue, c1); put!(queue, c2)
                finish_item!(); continue
            elseif st != "" && st != "ERROR"
                finish_item!(); continue                       # finished leaf from a previous run
            end
            leaf, dim = try
                procfun(path, lo, hi, W, pilot)
            catch err
                (Dict{String,Any}("path" => path, "lo" => lo, "hi" => hi, "status" => "ERROR", "why" => first(sprint(showerror, err), 300)), 0)
            end
            Threads.atomic_add!(nproc, 1)
            if leaf !== nothing
                lock(lk) do
                    println(io, jsonval(leaf)); flush(io)
                    stats[leaf["status"]] = get(stats, leaf["status"], 0) + 1
                end
            else
                lock(lk) do
                    println(io, jsonval(Dict{String,Any}("path" => path, "status" => "split", "dim" => dim))); flush(io)
                    stats["split"] = get(stats, "split", 0) + 1
                end
                c1, c2 = children(path, lo, hi, dim)
                Threads.atomic_add!(pending, 2); put!(queue, c1); put!(queue, c2)
            end
            if time() - lastreport[] > 60
                lock(lk) do
                    if time() - lastreport[] > 60
                        lastreport[] = time()
                        @printf(stderr, "[%s] processed %d, pending %d, leaves %s, elapsed %.0fs\n", Dates.format(now(), "HH:MM:SS"), nproc[], pending[], string(stats), time() - tstart)
                        flush(stderr)
                    end
                end
            end
            if nproc[] >= maxboxes && !stopped[]
                stopped[] = true
                for _ = 1:Threads.nthreads(); put!(queue, STOP); end
            end
            finish_item!()
        end
    end
    tasks = [Threads.@spawn worker_loop() for _ = 1:Threads.nthreads()]
    foreach(wait, tasks)
    close(io)
    @printf(stderr, "DONE processed %d, pending %d, leaves %s, elapsed %.0fs\n", nproc[], pending[], string(stats), time() - tstart)
    return stats
end

if abspath(PROGRAM_FILE) == @__FILE__
    logfile = ARGS[1]
    pilot = "--pilot" in ARGS
    maxb = (i = findfirst(==("--maxboxes"), ARGS)) === nothing ? typemax(Int) : parse(Int, ARGS[i+1])
    roots = root_boxes()
    if (i = findfirst(==("--region"), ARGS)) !== nothing
        v = parse.(Float64, split(ARGS[i+1], ","))
        roots = root_boxes([v[1], v[3], v[5], v[7]], [v[2], v[4], v[6], v[8]])
    end
    run_cover(logfile; pilot, maxboxes = maxb, roots)
end
