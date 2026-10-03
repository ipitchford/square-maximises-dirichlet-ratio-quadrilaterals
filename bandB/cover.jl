# cover.jl — adaptive certified cover of one chart domain (BAND_SPEC §4–§5).  Author: Anonymous.
#
# usage:  julia -t T cover.jl <chart 1|2> <log.jsonl> [--depth0 8] [--floor 24] [--only PREFIX] [--maxsub N]
#
# The root box is pre-split (round-robin bisection, depth0 levels) into subtrees; each subtree is refined depth-first
# by exact bisection.  Every leaf gets ONE status: OUTSIDE_CHART | LEMMA_A | LEMMA_B | OUT_OF_BAND (exact rational
# exclusion on the whole box, geom.jl) | CERT (certify(), cert.jl) | FAIL_FLOOR (a box that would have to be bisected
# below width 2^-floor: a LOUD failure, never dropped).  Records of a subtree are written together with a "done"
# marker, so a crashed run resumes by skipping finished subtrees.
include(joinpath(@__DIR__, "cert.jl"))
using Printf, Dates

qs(x::QQ) = string(numerator(x), "/", denominator(x))
boxjson(B::Box) = "\"lo\":[" * join(("\"" * qs(x) * "\"" for x in B.lo), ",") * "],\"hi\":[" *
                  join(("\"" * qs(x) * "\"" for x in B.hi), ",") * "]"
esc(s) = replace(string(s), "\\" => "\\\\", "\"" => "\\\"")

"bisection direction: largest float loss when that coordinate is collapsed to its midpoint (fallback: widest)"
function choose_dim(B::Box, prm::Params, r0::Float64)
    R = root_box(B.chart)
    relw = [Float64((B.hi[d] - B.lo[d]) / (R.hi[d] - R.lo[d])) for d = 1:4]
    widest = argmax(relw)
    (isfinite(r0) && r0 > 0) || return widest, "widest"
    mid = ntuple(d -> (B.lo[d] + B.hi[d]) / 2, 4)
    rc = try
        μ, R2, _, _ = float_eval(Box(B.chart, mid, mid), prm); R2 / μ
    catch
        Inf
    end
    (isfinite(rc) && rc < 2.45) || return widest, "widest(center $(round(rc, digits=4)))"
    loss = zeros(4)
    for d = 1:4
        relw[d] > 0 || continue
        lo = ntuple(k -> k == d ? mid[d] : B.lo[k], 4); hi = ntuple(k -> k == d ? mid[d] : B.hi[k], 4)
        rd = try
            μ, R2, _, _ = float_eval(Box(B.chart, lo, hi), prm); R2 / μ
        catch
            Inf
        end
        loss[d] = isfinite(rd) ? r0 - rd : -Inf
    end
    d = argmax(loss)
    loss[d] > 1e-9 || return widest, "widest(no loss)"
    return d, "loss"
end

function process_subtree(chart::Int, prefix::String, prm::Params, floorw::QQ, buf::Vector{String}, stats)
    stack = [(prefix, box_from_path(chart, prefix))]
    while !isempty(stack)
        path, B = pop!(stack)
        head = "{\"chart\":$chart,\"path\":\"$path\","
        ex = exclusion(B)
        if ex !== nothing
            push!(buf, head * "\"status\":\"$(ex[1])\",\"rule\":\"$(esc(ex[2]))\"," * boxjson(B) * "}")
            continue
        end
        t0 = time()
        r = try
            certify(B, prm)
        catch err
            (ok = false, E = NaN, R2lo = NaN, R2hi = NaN, nS = 0, nR = 0, why = "exception $(typeof(err))")
        end
        dt = time() - t0
        if r.ok
            margin = r.R2A / (Arb(r.E) * 5 / 2)
            push!(buf, head * "\"status\":\"CERT\"," * boxjson(B) * ",\"E\":\"$(repr(prevfloat(r.E)))\",\"R2\":\"$(updec(r.R2A))\"," *
                       "\"margin\":\"$(updec(margin))\",\"nS\":$(r.nS),\"nR\":$(r.nR),\"sec\":$(round(dt, digits=4))}")
            stats[:worst] = max(stats[:worst], up(margin))
            continue
        end
        r0 = (isfinite(r.R2hi) && isfinite(r.E) && r.E > 0) ? r.R2hi / r.E : Inf
        d, how = choose_dim(B, prm, r0)
        if B.hi[d] - B.lo[d] <= floorw || length(path) >= 120
            push!(buf, head * "\"status\":\"FAIL_FLOOR\"," * boxjson(B) * ",\"why\":\"$(esc(r.why))\",\"dim\":$d}")
            continue
        end
        L, U = bisect(B, d)
        push!(stack, (path * pathchar(d, true), U))
        push!(stack, (path * pathchar(d, false), L))
    end
end

function prefixes(depth0::Int)
    P = [""]
    for lev = 1:depth0
        d = mod1(lev, 4)
        P = [p * c for p in P for c in (pathchar(d, false), pathchar(d, true))]
    end
    return P
end

function main(args)
    chart = parse(Int, args[1]); logf = args[2]
    depth0 = 8; floorbits = 24; only = ""; maxsub = typemax(Int)
    i = 3
    while i <= length(args)
        a = args[i]
        if a == "--depth0"; depth0 = parse(Int, args[i+1]); i += 2
        elseif a == "--floor"; floorbits = parse(Int, args[i+1]); i += 2
        elseif a == "--only"; only = args[i+1]; i += 2
        elseif a == "--maxsub"; maxsub = parse(Int, args[i+1]); i += 2
        else error("unknown arg $a") end
    end
    prm = DEFAULT_PARAMS
    floorw = QQ(1, big(2)^floorbits)
    subs = [p for p in prefixes(depth0) if startswith(p, only) || startswith(only, p)]
    if only != "" && length(only) > depth0
        subs = [only]
    end
    subs = subs[1:min(end, maxsub)]
    # resume: keep only records of finished subtrees
    done = Set{String}()
    if isfile(logf)
        lines = readlines(logf)
        keep = String[]; recs = Dict{String,Vector{String}}()
        for ln in lines
            m = match(r"^\{\"subtree\":\"([A-Da-d]*)\",\"done\":true", ln)
            if m !== nothing
                push!(done, m.captures[1])
            end
        end
        for ln in lines
            m = match(r"\"path\":\"([A-Da-d]*)\"", ln)
            ms = match(r"^\{\"subtree\":\"([A-Da-d]*)\"", ln)
            key = ms !== nothing ? ms.captures[1] : (m === nothing ? nothing : m.captures[1])
            key === nothing && continue
            any(p -> startswith(key, p) && p in done, (key[1:min(end, k)] for k = 0:length(key))) && push!(keep, ln)
        end
        open(logf, "w") do io; foreach(l -> println(io, l), keep); end
        println("resume: $(length(done)) finished subtrees kept ($(length(keep)) lines)")
    end
    todo = [p for p in subs if !(p in done)]
    println("chart $chart: $(length(subs)) subtrees, $(length(todo)) to do, threads $(Threads.nthreads()), floor 2^-$floorbits")
    flush(stdout)
    lk = ReentrantLock()
    io = open(logf, "a")
    t0 = time()
    ch = Channel{String}(length(todo)); foreach(p -> put!(ch, p), todo); close(ch)
    nleaf = Threads.Atomic{Int}(0); nsub = Threads.Atomic{Int}(0)
    worst = Threads.Atomic{Float64}(0.0)
    @sync for _ = 1:Threads.nthreads()
        Threads.@spawn begin
            for p in ch
                buf = String[]; stats = Dict(:worst => 0.0)
                ts = time()
                process_subtree(chart, p, prm, floorw, buf, stats)
                lock(lk) do
                    for l in buf; println(io, l); end
                    println(io, "{\"subtree\":\"$p\",\"done\":true,\"leaves\":$(length(buf)),\"sec\":$(round(time()-ts, digits=2))}")
                    flush(io)
                end
                Threads.atomic_add!(nleaf, length(buf)); Threads.atomic_add!(nsub, 1)
                Threads.atomic_max!(worst, stats[:worst])
                @printf("[%s] subtree %-10s leaves %6d  (done %d/%d, leaves %d, %.0f s, worst margin %.5f)\n",
                        Dates.format(now(), "HH:MM:SS"), p, length(buf), nsub[], length(todo), nleaf[], time() - t0, worst[])
                flush(stdout)
            end
        end
    end
    close(io)
    wall = time() - t0
    println("finished chart $chart in $(round(wall, digits=1)) s")
    summarize(logf, chart, wall)
end

"append a summary line computed from ALL leaf records in the log (including resumed subtrees)"
function summarize(logf, chart, wall)
    cnt = Dict{String,Int}(); worst = 0.0; certsec = 0.0; maxdepth = 0; nsub = 0; subsec = 0.0; nS = 0; nR = 0
    for ln in eachline(logf)
        if startswith(ln, "{\"subtree\"")
            nsub += 1; subsec += parse(Float64, match(r"\"sec\":([0-9.e+-]+)", ln).captures[1]); continue
        end
        m = match(r"\"status\":\"([A-Z_]+)\"", ln); m === nothing && continue
        st = m.captures[1]; cnt[st] = get(cnt, st, 0) + 1
        maxdepth = max(maxdepth, length(match(r"\"path\":\"([A-Da-d]*)\"", ln).captures[1]))
        if st == "CERT"
            worst = max(worst, parse(Float64, match(r"\"margin\":\"([0-9.e+-]+)\"", ln).captures[1]))
            certsec += parse(Float64, match(r"\"sec\":([0-9.e+-]+)", ln).captures[1])
            nS = max(nS, parse(Int, match(r"\"nS\":(\d+)", ln).captures[1])); nR = max(nR, parse(Int, match(r"\"nR\":(\d+)", ln).captures[1]))
        end
    end
    s = "{\"summary\":{\"chart\":$chart,\"leaves\":$(sum(values(cnt))),\"counts\":{" *
        join(("\"$k\":$v" for (k, v) in sort(collect(cnt))), ",") * "},\"worst_margin_R2_over_2.5E\":\"$worst\"," *
        "\"max_depth\":$maxdepth,\"max_nS\":$nS,\"max_nR\":$nR,\"cert_seconds\":$(round(certsec, digits=1))," *
        "\"subtrees\":$nsub,\"subtree_thread_seconds\":$(round(subsec, digits=1)),\"wall_seconds_this_session\":$(round(wall, digits=1))," *
        "\"threads\":$(Threads.nthreads()),\"finished\":\"$(Dates.now())\"}}"
    open(logf, "a") do io; println(io, s); end
    println(s)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
