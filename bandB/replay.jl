# replay.jl — FAIL-CLOSED replay of the band cover logs (BAND_SPEC §5).  Author: Anonymous.
#
# usage:  julia -t T replay.jl [--cert-sample N] [--partial] [--out DIR] runs/cover_chart1.jsonl runs/cover_chart2.jsonl
# exit codes: 0 = PASSED (full replay of BOTH charts, every record verified);  1 = FAILED (anything wrong);
#             3 = a test-mode run (--partial and/or --cert-sample) with no failure found — NOT a proof.
#
# Fail-closed rules (any violation ⇒ FAILED, exit 1):
#   * every line of a log must match exactly one strict full-line schema: a leaf record of a WHITELISTED status
#     (CERT, LEMMA_A, LEMMA_B, OUT_OF_BAND, OUTSIDE_CHART — each with its own field list/order and formats; FAIL_FLOOR
#     is recognised but is itself a failure), a subtree-completion marker, or the single final summary line;
#     unknown statuses, unknown rules, malformed/extra/missing fields, non-canonical rationals ⇒ FAILED;
#   * the logged box must equal the box recomputed from the path (exact rationals);
#   * duplicate paths and duplicate lines ⇒ FAILED; the leaf paths must tile the chart root exactly (tiling.jl);
#   * "pending" work ⇒ FAILED: every leaf must lie in a subtree whose completion marker is present, the markers must
#     be exactly the 256 depth-8 prefixes, each once, with the right leaf counts; the summary counts must match;
#   * every exclusion is re-derived by the separate implementation ex_verify below (closed-form c_i and distances,
#     corner extrema of multilinear functions) and every CERT box is re-certified from scratch (cert.jl);
#   * both charts must be present (unless --partial, which forces exit 3).
# Binding: runs/replay_chart<c>.jsonl gets a header with the SHA-256 of every verifier source file and of the input
# log, one result line per leaf record carrying the SHA-256 of that exact log line, and a footer with the SHA-256
# digest of all result lines.
include(joinpath(@__DIR__, "cert.jl"))
include(joinpath(@__DIR__, "tiling.jl"))
using Printf, Dates, SHA

const SOURCES = ["replay.jl", "tiling.jl", "cert.jl", "cert_float.jl", "geom.jl"]
sha_file(f) = bytes2hex(open(sha256, f))
sha_str(s::AbstractString) = bytes2hex(sha256(s))

# ---------------- strict schemas ----------------
const RQ = raw"(0|[1-9][0-9]*)/([1-9][0-9]*)"
const RBOX = "\"lo\":\\[\"$RQ\",\"$RQ\",\"$RQ\",\"$RQ\"\\],\"hi\":\\[\"$RQ\",\"$RQ\",\"$RQ\",\"$RQ\"\\]"
const RFL = raw"-?[0-9]+\.[0-9]+(?:e-?[0-9]+)?"
const RE_EXCL = Regex("^\\{\"chart\":([12]),\"path\":\"([A-Da-d]*)\",\"status\":\"(OUTSIDE_CHART|LEMMA_A|LEMMA_B|OUT_OF_BAND)\",\"rule\":\"([^\"\\\\]*)\",$RBOX\\}\$")
const RE_CERT = Regex("^\\{\"chart\":([12]),\"path\":\"([A-Da-d]*)\",\"status\":\"CERT\",$RBOX,\"E\":\"($RFL)\",\"R2\":\"($RFL)\",\"margin\":\"($RFL)\",\"nS\":([0-9]+),\"nR\":([0-9]+),\"sec\":($RFL)\\}\$")
const RE_FAIL = Regex("^\\{\"chart\":([12]),\"path\":\"([A-Da-d]*)\",\"status\":\"FAIL_FLOOR\",$RBOX,\"why\":\"([^\"\\\\]*)\",\"dim\":([1-4])\\}\$")
const RE_MARK = Regex("^\\{\"subtree\":\"([A-Da-d]*)\",\"done\":true,\"leaves\":([0-9]+),\"sec\":($RFL)\\}\$")
# summary: every field enumerated in the order written by cover.jl: summarize(); no extra, missing or duplicate keys
const RSTAT = "(?:CERT|LEMMA_A|LEMMA_B|OUT_OF_BAND|OUTSIDE_CHART|FAIL_FLOOR)"
const RE_SUMM = Regex("^\\{\"summary\":\\{\"chart\":([12]),\"leaves\":([0-9]+),\"counts\":\\{(\"$RSTAT\":[0-9]+(?:,\"$RSTAT\":[0-9]+)*)\\}," *
    "\"worst_margin_R2_over_2\\.5E\":\"($RFL)\",\"max_depth\":([0-9]+),\"max_nS\":([0-9]+),\"max_nR\":([0-9]+),\"cert_seconds\":($RFL)," *
    "\"subtrees\":([0-9]+),\"subtree_thread_seconds\":($RFL),\"wall_seconds_this_session\":($RFL),\"threads\":([0-9]+)," *
    "\"finished\":\"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}(?:\\.[0-9]{1,3})?\"\\}\\}\$")

const RULES = Dict(
    (1, "OUTSIDE_CHART") => Set(["|P$(j)P$(k)| > 1" for j = 0:3 for k = j+1:3]),
    (2, "OUTSIDE_CHART") => Set(vcat(["|P$(j)P$(k)| > 1" for j = 0:3 for k = j+1:3], ["c$i <= 0" for i = 0:3], ["x3 >= x2"])),
    (1, "LEMMA_A") => Set(["c$i <= 7/200 c$(mod(i + 2, 4))" for i = 0:3]),
    (2, "LEMMA_A") => Set(["c$i <= 7/200 c$(mod(i + 2, 4))" for i = 0:3]),
    (1, "LEMMA_B") => Set(vcat(["y1+y3 <= 1/10"], ["width normal to edge $i <= 1/10" for i = 0:3])),
    (2, "LEMMA_B") => Set(["y2 <= 1/10"]),
    (1, "OUT_OF_BAND") => Set(["y3 > 2/5", "width lower bound > 2/5"]),
    (2, "OUT_OF_BAND") => Set(["y2 > 2/5"]))

"parse 8 canonical rationals from regex captures starting at index k0 (pairs num, den); nothing if not canonical"
function parse_box(caps, k0)
    vals = QQ[]
    for k = 0:7
        n = parse(BigInt, caps[k0+2k]); d = parse(BigInt, caps[k0+2k+1])
        gcd(n, d) == 1 || return nothing
        push!(vals, QQ(n, d))
    end
    return Tuple(vals[1:4]), Tuple(vals[5:8])
end

# ---------------- independent exclusion verification (direct formulas; corner extrema of multilinear functions) ----
function corner_pts(B)
    [(m & 1 == 1 ? B.hi[1] : B.lo[1], m & 2 == 2 ? B.hi[2] : B.lo[2], m & 4 == 4 ? B.hi[3] : B.lo[3], m & 8 == 8 ? B.hi[4] : B.lo[4]) for m = 0:15]
end
# closed forms (sympy_check.py 4a–4g):  chart I  p = (y3,t,x1,x3): c = (t y3, y3(1+t−x1−t x3), y3, y3(x1 + t x3));
#                                       chart II p = (y2,t,x2,x3): c = (y2, y2(1−x3) − t y2(1−x2), t y2 x2 − x3 y2, t y2)
cI(p) = (p[2] * p[1], p[1] * (1 + p[2] - p[3] - p[2] * p[4]), p[1], p[1] * (p[3] + p[2] * p[4]))
cII(p) = (p[1], p[1] * (1 - p[4]) - p[2] * p[1] * (1 - p[3]), p[2] * p[1] * p[3] - p[4] * p[1], p[2] * p[1])
vI(p) = ((QQ(0), QQ(0)), (p[3], -p[2] * p[1]), (QQ(1), QQ(0)), (p[4], p[1]))
vII(p) = ((QQ(0), QQ(0)), (QQ(1), QQ(0)), (p[3], p[1]), (p[4], p[2] * p[1]))
function ex_verify(chart, B, status, rule)
    P = corner_pts(B)
    cf = chart == 1 ? cI : cII; vf = chart == 1 ? vI : vII
    C = [cf(p) for p in P]; Vs = [vf(p) for p in P]
    cmax(i) = maximum(c[i+1] for c in C); cmin(i) = minimum(c[i+1] for c in C)
    function d2(j, k)            # (min, max) of |P_j − P_k|² over B; components multilinear
        dx = [V[j+1][1] - V[k+1][1] for V in Vs]; dy = [V[j+1][2] - V[k+1][2] for V in Vs]
        lo(a, b) = a <= 0 <= b ? QQ(0) : min(a * a, b * b)
        return lo(minimum(dx), maximum(dx)) + lo(minimum(dy), maximum(dy)), max(minimum(dx)^2, maximum(dx)^2) + max(minimum(dy)^2, maximum(dy)^2)
    end
    if status == "LEMMA_A"
        m = match(r"^c(\d) <= 7/200 c(\d)$", rule); m === nothing && return false
        i = parse(Int, m.captures[1]); j = parse(Int, m.captures[2]); j == mod(i + 2, 4) || return false
        return maximum(c[i+1] - QQ(7, 200) * c[j+1] for c in C) <= 0
    elseif status == "OUTSIDE_CHART"
        m = match(r"^\|P(\d)P(\d)\| > 1$", rule)
        m !== nothing && return d2(parse(Int, m.captures[1]), parse(Int, m.captures[2]))[1] > 1
        m = match(r"^c(\d) <= 0$", rule)
        m !== nothing && return chart == 2 && cmax(parse(Int, m.captures[1])) <= 0
        rule == "x3 >= x2" && return chart == 2 && B.lo[4] >= B.hi[3]
        return false
    elseif status == "LEMMA_B"
        chart == 1 && rule == "y1+y3 <= 1/10" && return B.hi[1] * (1 + B.hi[2]) <= QQ(1, 10)
        chart == 2 && rule == "y2 <= 1/10" && return B.hi[1] <= QQ(1, 10)
        m = match(r"^width normal to edge (\d) <= 1/10$", rule)
        if chart == 1 && m !== nothing
            i = parse(Int, m.captures[1]); mc = max(cmax(i), cmax(mod(i + 3, 4)))
            return mc >= 0 && mc^2 <= QQ(1, 100) * d2(i, mod(i + 1, 4))[1]
        end
        return false
    elseif status == "OUT_OF_BAND"
        chart == 1 && rule == "y3 > 2/5" && return B.lo[1] > QQ(2, 5)
        chart == 2 && rule == "y2 > 2/5" && return B.lo[1] > QQ(2, 5)
        if chart == 1 && rule == "width lower bound > 2/5"
            return all(i -> (m = max(cmin(i), cmin(mod(i + 3, 4))); m > 0 && m^2 > QQ(4, 25) * d2(i, mod(i + 1, 4))[2]), 0:3)
        end
        return false
    end
    return false
end

function depth8_prefixes()
    P = [""]
    for lev = 1:8
        d = mod1(lev, 4)
        P = [p * c for p in P for c in (pathchar(d, false), pathchar(d, true))]
    end
    return P
end

"replay one log; returns (ok::Bool, chart, report line)"
function replay_log(f::String, cert_sample::Int, sources_sha::Vector{Pair{String,String}}, outdir::String)
    errs = String[]
    err!(m) = (push!(errs, m); length(errs) <= 20 && println("FAIL: ", m))
    lines = readlines(f)
    logsha = sha_file(f)
    println("\n== $f  ($(length(lines)) lines, sha256 $logsha)")
    isempty(lines) && (err!("empty log"); return false, 0, "$f: empty log")
    # ---------- phase 1: schemas, uniqueness, boxes, markers, summary (sequential; no long-lived BigInts) ----------
    kinds = fill(:none, length(lines))
    leafidx = Int[]; paths = String[]; stats = String[]; rules = String[]
    marks = Dict{String,Int}(); charts = Set{Int}()
    seenline = Set{String}(); seenpath = Set{String}()
    summ = nothing
    for (n, ln) in enumerate(lines)
        ln in seenline && (err!("line $n: duplicate line"); continue)
        push!(seenline, ln)
        m = match(RE_CERT, ln); st = "CERT"; rule = ""; k0 = 3
        if m === nothing
            m = match(RE_EXCL, ln)
            if m !== nothing
                st = m.captures[3]; rule = m.captures[4]; k0 = 5
            end
        end
        if m === nothing
            mf = match(RE_FAIL, ln)
            if mf !== nothing
                err!("line $n: FAIL_FLOOR leaf $(mf.captures[2]) — refined below the floor, NOT PROVED"); kinds[n] = :failfloor; continue
            end
            mk = match(RE_MARK, ln)
            if mk !== nothing
                p = mk.captures[1]
                haskey(marks, p) && err!("line $n: duplicate subtree marker $p")
                marks[p] = parse(Int, mk.captures[2]); kinds[n] = :marker; continue
            end
            ms = match(RE_SUMM, ln)
            if ms !== nothing
                n == length(lines) || err!("line $n: summary line is not the last line")
                summ === nothing || err!("line $n: second summary line")
                summ = ms; kinds[n] = :summary; continue
            end
            err!("line $n: does not match any whitelisted record schema: $(first(ln, 160))"); continue
        end
        chart = parse(Int, m.captures[1]); path = String(m.captures[2])
        push!(charts, chart)
        if st != "CERT"
            rule in get(RULES, (chart, st), Set{String}()) || (err!("line $n: rule '$rule' not whitelisted for $st in chart $chart"); continue)
        end
        pb = parse_box(m.captures, k0)
        pb === nothing && (err!("line $n: non-canonical rational"); continue)
        path in seenpath && (err!("line $n: duplicate path $path"); continue)
        push!(seenpath, path)
        B = box_from_path(chart, path)
        (B.lo == pb[1] && B.hi == pb[2]) || (err!("line $n: logged box differs from the box of path $path"); continue)
        kinds[n] = :leaf
        push!(leafidx, n); push!(paths, path); push!(stats, st); push!(rules, rule)
    end
    length(charts) == 1 || (err!("leaf records of $(length(charts)) charts in one log"); return false, 0, "$f: chart mixture")
    chart = first(charts)
    # completion markers: exactly the 256 depth-8 prefixes, every leaf inside a completed subtree, counts match
    P8 = depth8_prefixes()
    Set(keys(marks)) == Set(P8) || err!("subtree markers are not exactly the 256 depth-8 prefixes ($(length(marks)) markers)")
    cntp = Dict{String,Int}()
    for p in paths
        length(p) >= 8 || (err!("leaf path $p shorter than the pre-split depth"); continue)
        haskey(marks, p[1:8]) || (err!("leaf $p belongs to a subtree without completion marker (pending)"); continue)
        cntp[p[1:8]] = get(cntp, p[1:8], 0) + 1
    end
    for (p, c) in marks
        get(cntp, p, 0) == c || err!("subtree $p: marker says $c leaves, log has $(get(cntp, p, 0))")
    end
    # summary must exist and agree
    cnt = Dict{String,Int}(); for s in stats; cnt[s] = get(cnt, s, 0) + 1; end
    if summ === nothing
        err!("missing summary line")
    else
        parse(Int, summ.captures[1]) == chart || err!("summary chart mismatch")
        parse(Int, summ.captures[2]) == length(paths) || err!("summary leaf count mismatch")
        pairs = [(String(x.captures[1]), parse(Int, x.captures[2])) for x in eachmatch(r"\"([A-Z_]+)\":([0-9]+)", summ.captures[3])]
        keysq = first.(pairs)
        all(k -> keysq[k] < keysq[k+1], 1:length(keysq)-1) || err!("summary count keys not strictly sorted (duplicate or reordered keys)")
        sc = Dict(pairs)
        sc == cnt || err!("summary counts $sc differ from the records $cnt")
        parse(Int, summ.captures[9]) == length(marks) || err!("summary subtree count mismatch")
        parse(Int, summ.captures[5]) == maximum(length, paths; init = 0) || err!("summary max_depth mismatch")
    end
    # exact tiling
    okt, msgt = check_tiling(chart, paths)
    println("tiling: ", okt ? "PASS " : "FAIL ", msgt); okt || err!("tiling: $msgt")
    println("phase 1 (schemas, boxes, markers, summary, tiling): $(isempty(errs) ? "ok" : "$(length(errs)) errors"); statuses $cnt")
    # ---------- phase 2: re-derive every exclusion, re-run every certificate (threaded) ----------
    nleaf = length(paths)
    docert = trues(nleaf)
    if cert_sample > 0
        ic = findall(==("CERT"), stats)
        keep = Set(ic[unique(rand(1:length(ic), min(cert_sample, length(ic))))])
        for k in ic; docert[k] = k in keep; end
    end
    results = Vector{String}(undef, nleaf)
    verdict = fill(false, nleaf)
    worst = zeros(nleaf)
    t2 = time()
    Threads.@threads :dynamic for k = 1:nleaf
        n = leafidx[k]; path = paths[k]; st = stats[k]
        rec = sha_str(lines[n])
        B = box_from_path(chart, path)
        if st == "CERT"
            if !docert[k]
                results[k] = "{\"line\":$n,\"rec_sha256\":\"$rec\",\"path\":\"$path\",\"status\":\"CERT\",\"verdict\":\"NOT_REPLAYED\"}"
                verdict[k] = true
                continue
            end
            c = try
                certify(B)
            catch e
                (ok = false, why = "exception $(typeof(e))", E = NaN, R2A = Arb(NaN))
            end
            if c.ok
                mg = c.R2A / (Arb(c.E) * 5 / 2)
                worst[k] = up(mg); verdict[k] = true
                results[k] = "{\"line\":$n,\"rec_sha256\":\"$rec\",\"path\":\"$path\",\"status\":\"CERT\",\"verdict\":\"PASS\",\"E\":\"$(repr(c.E))\",\"R2_hi\":\"$(updec(c.R2A))\",\"margin_hi\":\"$(updec(mg))\"}"
            else
                results[k] = "{\"line\":$n,\"rec_sha256\":\"$rec\",\"path\":\"$path\",\"status\":\"CERT\",\"verdict\":\"FAIL\",\"why\":\"$(c.why)\"}"
            end
        else
            ok = try ex_verify(chart, B, st, rules[k]) catch; false end
            verdict[k] = ok
            results[k] = "{\"line\":$n,\"rec_sha256\":\"$rec\",\"path\":\"$path\",\"status\":\"$st\",\"rule\":\"$(rules[k])\",\"verdict\":\"$(ok ? "PASS" : "FAIL")\"}"
        end
    end
    nfail = count(!, verdict)
    for k in findall(!, verdict)
        err!("re-verification failed: line $(leafidx[k]) path $(paths[k]) $(stats[k])")
    end
    nex = count(s -> s != "CERT", stats); ncr = count(k -> stats[k] == "CERT" && docert[k], 1:nleaf)
    wm = maximum(worst; init = 0.0)
    @printf("phase 2: re-derived exclusions %d, re-run certificates %d (worst margin R2/(2.5E) ≤ %.7f), failures %d, %.1f s\n",
            nex, ncr, wm, nfail, time() - t2)
    # ---------- bound result file ----------
    outf = joinpath(outdir, "replay_chart$(chart).jsonl")
    digest = SHA.SHA256_CTX()
    open(outf, "w") do io
        hdr = "{\"replay_header\":{\"input_log\":\"$(basename(f))\",\"input_log_sha256\":\"$logsha\",\"verifier_sources\":{" *
              join(("\"$k\":\"$v\"" for (k, v) in sources_sha), ",") * "},\"threads\":$(Threads.nthreads()),\"cert_sample\":$cert_sample,\"started\":\"$(Dates.now())\"}}"
        println(io, hdr); SHA.update!(digest, codeunits(hdr * "\n"))
        for r in results
            println(io, r); SHA.update!(digest, codeunits(r * "\n"))
        end
        ok = isempty(errs)
        ftr = "{\"replay_result\":{\"chart\":$chart,\"verdict\":\"$(ok ? (cert_sample > 0 ? "NO_FAILURE_SAMPLED" : "PASSED") : "FAILED")\",\"leaves\":$nleaf," *
              "\"counts\":{" * join(("\"$k\":$v" for (k, v) in sort(collect(cnt))), ",") * "},\"exclusions_rederived\":$nex,\"certificates_rerun\":$ncr," *
              "\"worst_margin_R2_over_2.5E_hi\":\"$(repr(wm))\",\"errors\":$(length(errs)),\"results_sha256\":\"$(bytes2hex(SHA.digest!(digest)))\"}}"
        println(io, ftr); println(ftr)
    end
    rep = @sprintf("chart %d: %d leaves; tiling %s; CERT %d (re-run %d), LEMMA_A %d, LEMMA_B %d, OUT_OF_BAND %d, OUTSIDE_CHART %d; errors %d; worst margin ≤ %.7f",
                   chart, nleaf, okt ? "exact" : "FAILED", get(cnt, "CERT", 0), ncr, get(cnt, "LEMMA_A", 0), get(cnt, "LEMMA_B", 0),
                   get(cnt, "OUT_OF_BAND", 0), get(cnt, "OUTSIDE_CHART", 0), length(errs), wm)
    return isempty(errs), chart, rep
end

function main(args)
    cert_sample = 0; partial = false; files = String[]; outdir = joinpath(@__DIR__, "runs")
    i = 1
    while i <= length(args)
        if args[i] == "--cert-sample"; cert_sample = parse(Int, args[i+1]); i += 2
        elseif args[i] == "--partial"; partial = true; i += 1
        elseif args[i] == "--out"; outdir = args[i+1]; i += 2
        else push!(files, args[i]); i += 1 end
    end
    t0 = time()
    sources_sha = [s => sha_file(joinpath(@__DIR__, s)) for s in SOURCES]
    println("replay started ", Dates.now(), "  threads ", Threads.nthreads(), "  cert_sample ", cert_sample, "  partial ", partial)
    for (k, v) in sources_sha; println("  source sha256 $k $v"); end
    allok = true; seen = Int[]; reps = String[]
    for f in files
        ok, chart, rep = try
            replay_log(f, cert_sample, sources_sha, outdir)
        catch e
            println("FAIL: exception while replaying $f: ", sprint(showerror, e)); (false, 0, "$f: exception")
        end
        allok &= ok; push!(seen, chart); push!(reps, rep)
    end
    if !partial
        sort(seen) == [1, 2] || (println("FAIL: need exactly one log for chart 1 and one for chart 2, got charts $seen"); allok = false)
    end
    println("\n", join(reps, "\n"))
    code = !allok ? 1 : ((partial || cert_sample > 0) ? 3 : 0)
    @printf("\nREPLAY %s  (%.1f s)\n", code == 0 ? "PASSED" : (code == 1 ? "FAILED" : "INCOMPLETE TEST RUN (no failure found; NOT a proof)"), time() - t0)
    return code
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
