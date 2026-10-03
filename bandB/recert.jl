# recert.jl — re-certification of an existing cover log with the CURRENT certificate code.  Author: Anonymous.
#
# usage:  julia -t T recert.jl <chart> <old_log.jsonl> <new_log.jsonl>
#
# Written after the 2026-10-01 review found that the original cover's Gram assembly rounded nodal products in Float64
# (cert.jl, fixed).  Every leaf path of the old log is re-processed from scratch by cover.jl's process_subtree: the
# exclusion is re-derived exactly, a CERT leaf is re-certified with the fixed cert.jl, and any leaf that no longer
# passes is bisected further exactly as in the cover (FAIL_FLOOR below 2^-24).  The new log has the cover format (leaf
# records, one completion marker per depth-8 subtree, one summary line) and is checked by replay.jl.  Re-subdivided
# old leaves are listed in the .out file.  The old log is only read.
include(joinpath(@__DIR__, "cover.jl"))

function recert_main(args)
    chart = parse(Int, args[1]); oldf = args[2]; newf = args[3]
    isfile(newf) && error("refusing to overwrite existing $newf")
    groups = Dict{String,Vector{Tuple{String,String}}}()
    for ln in eachline(oldf)
        m = match(r"^\{\"chart\":(\d),\"path\":\"([A-Da-d]*)\",\"status\":\"([A-Z_]+)\"", ln)
        m === nothing && continue
        parse(Int, m.captures[1]) == chart || error("chart mismatch in $oldf")
        p = String(m.captures[2]); st = String(m.captures[3])
        st == "FAIL_FLOOR" && error("old log contains FAIL_FLOOR")
        length(p) >= 8 || error("old leaf $p shorter than the pre-split depth")
        push!(get!(groups, p[1:8], Tuple{String,String}[]), (p, st))
    end
    P8 = prefixes(8)
    Set(keys(groups)) == Set(P8) || error("old log does not cover all 256 depth-8 subtrees")
    prm = DEFAULT_PARAMS; floorw = QQ(1, big(2)^24)
    lk = ReentrantLock(); io = open(newf, "w")
    t0 = time()
    ch = Channel{String}(length(P8)); foreach(p -> put!(ch, p), P8); close(ch)
    nold = Threads.Atomic{Int}(0); nnew = Threads.Atomic{Int}(0); nres = Threads.Atomic{Int}(0); nsub = Threads.Atomic{Int}(0)
    worst = Threads.Atomic{Float64}(0.0)
    resub = String[]
    @sync for _ = 1:Threads.nthreads()
        Threads.@spawn for P in ch
            buf = String[]; stats = Dict(:worst => 0.0); ts = time(); myres = String[]
            for (p, st) in groups[P]
                n0 = length(buf)
                process_subtree(chart, p, prm, floorw, buf, stats)
                k = length(buf) - n0
                newst = match(r"\"status\":\"([A-Z_]+)\"", buf[n0+1]).captures[1]
                if k != 1 || !occursin("\"path\":\"$p\"", buf[n0+1])
                    push!(myres, "$p ($st) -> $k leaves")
                elseif newst != st
                    push!(myres, "$p status $st -> $newst")
                end
            end
            lock(lk) do
                for l in buf; println(io, l); end
                println(io, "{\"subtree\":\"$P\",\"done\":true,\"leaves\":$(length(buf)),\"sec\":$(round(time()-ts, digits=2))}")
                flush(io)
                append!(resub, myres)
            end
            Threads.atomic_add!(nold, length(groups[P])); Threads.atomic_add!(nnew, length(buf))
            Threads.atomic_add!(nres, length(myres)); Threads.atomic_add!(nsub, 1); Threads.atomic_max!(worst, stats[:worst])
            @printf("[%s] subtree %s: %d old leaves -> %d new leaves, %d changed  (done %d/256, %.0f s, worst %.7f)\n",
                    Dates.format(now(), "HH:MM:SS"), P, length(groups[P]), length(buf), length(myres), nsub[], time() - t0, worst[])
            flush(stdout)
        end
    end
    close(io)
    wall = time() - t0
    println("recert chart $chart: old leaves $(nold[]) -> new leaves $(nnew[]); changed (re-subdivided or status changed): $(nres[])")
    for r in sort(resub); println("  CHANGED: ", r); end
    summarize(newf, chart, wall)
end

if abspath(PROGRAM_FILE) == @__FILE__
    recert_main(ARGS)
end
