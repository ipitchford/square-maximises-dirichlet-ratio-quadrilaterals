include("../replay_t.jl"); include("../xcheck_cf.jl")
using Printf, Random
rng = MersenneTwister(3)
recs = Dict{String,Dict{String,Any}}()
for ln in eachline("runs/transport_cert.jsonl"); r = parse_rec(ln); recs[r["path"]] = r; end
res = [parse_rec(ln) for ln in eachline("runs/transport_cert_replay_results.jsonl") if endswith(ln, "\"end\":1}")]
cvals_f(p) = cvals(p)
function holds_pt(kind, cond, p)
    c = cvals(p); V = verts(p)
    minimum(c) <= 0 && kind != "E4" && return true          # non-convex points are excluded anyway (E3/E2 tests are polynomial; just report)
    if kind == "E4"; return c[cond[1]] <= 1e-12; end
    if kind == "E2"; i = cond[1]; return c[i] - 0.035 * c[mod1(i + 2, 4)] <= 1e-12; end
    i, kk, ll = cond; s = kind == "E3" ? 0.1 : 0.4
    ex = V[mod1(i + 1, 4), 1] - V[i, 1]; ey = V[mod1(i + 1, 4), 2] - V[i, 2]
    e2 = ex^2 + ey^2; d2 = (V[kk, 1] - V[ll, 1])^2 + (V[kk, 2] - V[ll, 2])^2
    all(j -> (ex * (V[j, 2] - V[i, 2]) - ey * (V[j, 1] - V[i, 1]))^2 - s^2 * e2 * d2 <= 1e-12 * (1 + e2 * d2), (mod1(i + 2, 4), mod1(i + 3, 4)))
end
function subleaves(tree, lo, hi, out)
    tree = strip(tree)
    if startswith(tree, "E"); push!(out, (String(tree), lo, hi)); return out; end
    dim = findfirst(==(tree[1]), DIMNAME); inner = tree[3:end-1]; depth = 0; cut = 0
    for (i, ch) in enumerate(inner); ch == '(' && (depth += 1); ch == ')' && (depth -= 1); if ch == ',' && depth == 0; cut = i; break; end; end
    mid = (lo[dim] + hi[dim]) / 2; h1 = copy(hi); h1[dim] = mid; l2 = copy(lo); l2[dim] = mid
    subleaves(inner[1:cut-1], lo, h1, out); subleaves(inner[cut+1:end], l2, hi, out)
end
check_E_cf_which("E2", [0.3, 0.3, 0.3, 0.3], [0.31, 0.31, 0.31, 0.31])   # compile
stats = Dict{String,Vector{Float64}}()   # kind => [n, cf ok, old tm-only, cf ok among tm-only, time, unsound]
for q in res[1:min(end, 1500)]
    q["kind"] == "certified" && continue
    r = recs[q["path"]]
    items = r["status"] == "COMPOUND" ? subleaves(r["tree"], r["lo"], r["hi"], Tuple{String,Vector{Float64},Vector{Float64}}[]) : [(r["status"], r["lo"], r["hi"])]
    for (kind, lo, hi) in items
        s = get!(stats, kind, zeros(6))
        t = @elapsed (w = check_E_cf_which(kind, lo, hi))
        s[1] += 1; s[5] += t
        if w !== nothing
            s[2] += 1
            for _ = 1:100
                p = lo .+ (hi .- lo) .* rand(rng, 4)
                holds_pt(kind, w, p) || (s[6] += 1; println("UNSOUND? ", kind, " ", w, " p=", p); break)
            end
        end
    end
    if q["tm"] > 0
        s = stats[r["status"] == "COMPOUND" ? "E2" : r["status"]]
    end
end
for (k, s) in stats
    @printf("%-5s sub-leaves %5d  proved by centred form %5d (%.1f%%)  time %.4f s each  unsound samples %d\n", k, s[1], s[2], 100s[2]/s[1], s[5]/s[1], s[6])
end
# soundness on certified boxes (contain points of R): any proved condition must hold at sampled points
nproved = 0; nbad = 0
for q in res[1:min(end, 1500)]
    q["kind"] == "certified" || continue
    r = recs[q["path"]]
    for kind in ("E2", "E3T", "E4")
        w = check_E_cf_which(kind, r["lo"], r["hi"]; depth = 2, budget = 20)
        w === nothing && continue
        global nproved += 1
        for _ = 1:200
            p = r["lo"] .+ (r["hi"] .- r["lo"]) .* rand(rng, 4)
            holds_pt(kind, w, p) || (global nbad += 1; println("UNSOUND on certified box ", r["path"], " ", kind, w); break)
        end
    end
end
@printf("certified boxes: conditions proved by cf %d, unsound %d\n", nproved, nbad)
