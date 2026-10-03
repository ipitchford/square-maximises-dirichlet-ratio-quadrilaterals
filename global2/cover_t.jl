# cover_t.jl — the transport cover of the region outside the certified shell (TRANSPORT_SPEC.md).
#
# usage:  julia -t 8 cover_t.jl <logfile.jsonl> [--pilot] [--maxboxes N] [--region a0,a1,b0,b1,c0,c1,d0,d1]
#
# Region: the a-priori box (a ≤ 0.905, b, c ≤ 4.36, |d| ≤ 0.905) of the fundamental domain, split at a, b, c = S and
# d = ±S with S = Float64(1/5); the inner root [0,S]³×[−S,S] lies inside the shell certified by cover.jl
# (runs/shell0225_cert.jsonl, 0.049 < |p|∞ ≤ 0.225, plus the local cube) and is logged with status "SHELL".
# Every other box is (i) discarded by one exclusion condition E1–E4 (Arb Taylor models), or (ii) discarded by the
# COMPOUND test (an internal bisection of depth ≤ CDEPTH whose every leaf satisfies one exclusion condition; the
# leaf tree is logged), or (iii) certified by the affine-transport certificate (transport.jl), or (iv) bisected.
include("cover.jl")
include("transport.jl")
include("bounds_tm.jl")
include("transport1.jl")
set_order!(4)          # order-4 Taylor models (order 3 folds too much at radius 0.04)

const S_SHELL = Float64(1 // 5)        # inner root [0,S]^3 × [−S,S] ⊂ certified shell region (S ≤ 0.225)
# DEFERRED thin band (NOT an established lemma): boxes on which width/diam ≤ S1 at every point are logged "E3T" and
# left to a strengthened thin lemma (1-D reduction, see README §11).  S1 = 2/5 exactly.
const RAT_S1 = 2 // 5

"thin test with an arbitrary threshold s (exact rational): width/diam ≤ s on the whole box (same polynomial test as E3)"
function thin_test(lo, hi, srat::Rational, ws)
    p0, rad = box_center_radius(lo, hi)
    G = box_geom(Arb, p0, rad, ws)
    s2 = tsc(Arb, srat^2)
    Vf = verts(p0)
    pairs = [(k, l) for k = 1:4 for l = k+1:4]
    sort!(pairs, by = kl -> -hypot(Vf[kl[1], 1] - Vf[kl[2], 1], Vf[kl[1], 2] - Vf[kl[2], 2]))
    e2 = newtm(Arb); d2 = newtm(Arb); dx = newtm(Arb); dy = newtm(Arb); rhs = newtm(Arb); lhs = newtm(Arb)
    for (k, l) in pairs[1:2]
        tm_sub!(dx, G.V[k, 1], G.V[l, 1]); tm_sub!(dy, G.V[k, 2], G.V[l, 2]); tm_sq2!(d2, dx, dy, ws)
        for i = 1:4
            tm_sq2!(e2, G.E[i, 1], G.E[i, 2], ws)
            tm_mul_full!(rhs, e2, d2, ws); tm_scale!(rhs, rhs, s2, ws)
            okall = true
            for j in (mod1(i + 2, 4), mod1(i + 3, 4))
                tm_mul_full!(lhs, G.A[i, j], G.A[i, j], ws); tm_sub!(lhs, lhs, rhs)
                if !(tm_upper(lhs) <= 0); okall = false; break; end
            end
            okall && return true
        end
    end
    return false
end
"single exclusion conditions of the transport cover: E1–E4, then the deferred thin band E3T"
function exclusion_t(lo, hi; ws = newwork(Arb))
    ex = exclusion(lo, hi; cube = CUBE, ws)
    ex !== :none && return ex
    thin_test(lo, hi, RAT_S1, ws) && return :E3T
    return :none
end
@assert Rational{BigInt}(S_SHELL) <= 225 // 1000
const WTRY_T = 0.15
const MINW_T = 2.0^-11
const CDEPTH = 6

function root_boxes_t(lo0 = A_PRIORI[1], hi0 = A_PRIORI[2])
    s = S_SHELL
    cut(l, h, pts) = (xs = vcat(l, [x for x in pts if l < x < h], h); [(xs[k], xs[k+1]) for k = 1:length(xs)-1])
    ivs = [cut(lo0[i], hi0[i], [s]) for i = 1:3]
    push!(ivs, cut(lo0[4], hi0[4], [-s, s]))
    roots = Tuple{String,Vector{Float64},Vector{Float64}}[]
    k = 0
    for ia in ivs[1], ib in ivs[2], ic in ivs[3], id in ivs[4]
        k += 1
        push!(roots, ("T$k", [ia[1], ib[1], ic[1], id[1]], [ia[2], ib[2], ic[2], id[2]]))
    end
    roots
end
"the inner root (inside the certified shell region)"
is_shell_box(lo, hi) = all(lo[i] >= 0 && Rational{BigInt}(hi[i]) <= 225 // 1000 for i = 1:3) &&
                       Rational{BigInt}(lo[4]) >= -225 // 1000 && Rational{BigInt}(hi[4]) <= 225 // 1000

# ------------------------------------------------------------------ compound exclusion
"""
    compound(lo, hi, depth, ws) -> Union{Nothing,String}

Bisect (widest side) up to `depth` levels; succeed if every leaf satisfies one single exclusion condition.
Returns the leaf tree as a string: a leaf is its condition (E1..E4), a node is "<dim>(<left>,<right>)".
"""
function compound(lo, hi, depth, ws)
    ex = exclusion_t(lo, hi; ws)
    ex !== :none && return String(ex)
    depth == 0 && return nothing
    dim = argmax(hi .- lo)
    mid = (lo[dim] + hi[dim]) / 2
    (lo[dim] < mid < hi[dim]) || return nothing
    hi1 = copy(hi); hi1[dim] = mid; lo2 = copy(lo); lo2[dim] = mid
    l = compound(lo, hi1, depth - 1, ws); l === nothing && return nothing
    r = compound(lo2, hi, depth - 1, ws); r === nothing && return nothing
    return string(DIMNAME[dim], "(", l, ",", r, ")")
end

# ------------------------------------------------------------------ workspace
struct TWorker
    tpF::TPool{Float64,Float64}; tpA::TPool{Arb,Mag}
    cpF::CritPool{Float64,Float64}; cpA::CritPool{Arb,Mag}
    wsF::TMWork{Float64,Float64}; wsA::TMWork{Arb,Mag}
end
TWorker() = TWorker(tpool(Float64), tpool(Arb), critpool(Float64), critpool(Arb), newwork(Float64), newwork(Arb))

function tpoint_kappa(R::TRefs{Float64}, beta, L0)
    K, M, A2 = tpoint(R)
    th = eigvals(Symmetric(K), Symmetric(M))
    LG = beta * M[1, 1] > K[1, 1] ? (beta * M[1, 1] - K[1, 1]) / (beta * A2 - M[1, 1]) : 0.0
    return th[2] / max(LG, L0)
end

function kappa_star_t!(cp, tp, ws, mode, beta, L0; steps = 12)
    lo = 1.0; hi = 2.5
    for _ = 1:steps
        mid = (lo + hi) / 2
        rr = crit_eval!(cp, tp.M, tp.K, tp.A2, ws; mode, beta, L0, kappa = rationalize(BigInt, mid; tol = 1e-12))
        rr.ok ? (hi = mid) : (lo = mid)
    end
    hi
end

# Float64 screen of one (n, shift) choice.  Lower-bound candidates: the box constants (β0, L0, and the
# Andrews–Clutterbuck bootstrap from the LG bound) and the pointwise Taylor-model bounds of bounds_tm.jl.
# Returns the best (smallest certified κ*) candidate: (; ok, kappa, mode, bname, beta, L0, sol)
function screen(p0, rad, n, shift, sb, bt, W::TWorker; variants = (true, false))
    L0c = max(sb.L_rect, sb.L_hersch)
    beta0 = max(sb.beta_rect_ks, addlo(L0c, sb.gapAC))
    sol1 = tsolve1(p0, n; shift)                     # first-order corrected choice (Float64)
    best = (; ok = false, kappa = Inf, mode = :LG, bname = "", beta = beta0, L0 = L0c, sol = sol1, fo = true, dim = 0)
    for fo in variants
        sol = fo ? sol1 : merge(sol1, (; dC = zero(sol1.dC), ddpsi = zero(sol1.ddpsi)))   # plain transport
        R = transport1_refs(Float64, sol)
        b = transport1_base!(W.tpF, R, p0, rad, W.wsF, shift)
        b.ok || continue
        base = (; M = W.tpF.M, K = W.tpF.K, A2 = W.tpF.A2)
        # candidate order = observed usage: constants (LG), Krahn–Szegő, constants (MONO), strips, bootstrap, rest
        val0(x) = x[2].c[1]
        cands = Tuple{Symbol,String,Any,Any}[(:LG, "const", beta0, L0c)]
        for (nm, B) in bt.betas; nm == "KS" && push!(cands, (:LG, nm, B, L0c)); end
        push!(cands, (:MONO, "const", beta0, L0c))
        for (nm, L) in sort(filter(x -> startswith(x[1], "strip"), bt.L0s), by = x -> -val0(x)); push!(cands, (:MONO, nm, beta0, L)); end
        r0 = crit_eval!(W.cpF, base.M, base.K, base.A2, W.wsF; mode = :LG, beta = beta0, L0 = L0c)
        if best.dim == 0                               # split suggestion from the first combination
            f = !r0.g1ok ? W.cpF.g2 : (r0.a11up >= 0 ? W.cpF.At[1] : W.cpF.s11)
            best = merge(best, (; dim = argmax(influence(f))))
        end
        if r0.g1ok
            Llo = r0.g2lo / tm_upper(W.cpF.g1)
            Llo + sb.gapAC > beta0 && push!(cands, (:LG, "boot", Llo + sb.gapAC, L0c))
        end
        for (nm, B) in sort(filter(x -> x[1] != "KS", bt.betas), by = x -> -val0(x)); push!(cands, (:LG, nm, B, L0c)); end
        for (nm, L) in sort(filter(x -> !startswith(x[1], "strip"), bt.L0s), by = x -> -val0(x)); push!(cands, (:MONO, nm, beta0, L)); end
        for (mode, nm, beta, L0) in cands
            rr = crit_eval!(W.cpF, base.M, base.K, base.A2, W.wsF; mode, beta, L0)
            rr.ok || continue
            k = kappa_star_t!(W.cpF, base, W.wsF, mode, beta, L0; steps = 5)
            return (; ok = true, kappa = k, mode, bname = nm, beta, L0, sol, fo, dim = best.dim)
        end
        # pointwise bootstrap β2(p) = L_LG(p; β1) + 3π²/D̂(p)², β1 ∈ {KS, const, other β candidates}, iterated twice
        for (nm1, B1) in vcat([("const", beta0)], [(x[1], x[2]) for x in bt.betas])
            b2 = boot_tm(W.cpF, base, B1, bt.gap, W.wsF); b2 === nothing && continue
            nmb = "bootTM:" * nm1
            rr = crit_eval!(W.cpF, base.M, base.K, base.A2, W.wsF; mode = :LG, beta = b2, L0 = L0c)
            if !rr.ok
                b3 = boot_tm(W.cpF, base, b2, bt.gap, W.wsF); b3 === nothing && continue
                rr = crit_eval!(W.cpF, base.M, base.K, base.A2, W.wsF; mode = :LG, beta = b3, L0 = L0c)
                rr.ok || continue
                b2 = b3; nmb = "bootTM2:" * nm1
            end
            k = kappa_star_t!(W.cpF, base, W.wsF, :LG, b2, L0c; steps = 5)
            return (; ok = true, kappa = k, mode = :LG, bname = nmb, beta = b2, L0 = L0c, sol, fo, dim = best.dim)
        end
    end
    return best
end

"""
    arb_bound(sc, sb, beta0, L0, p0, rad, W) -> (betaA, L0A, ok)

Rigorous (Arb) version of the lower-bound candidate chosen in the Float64 screen: the box constants, the
Andrews–Clutterbuck bootstrap recomputed from the Arb LG bound, or the named Taylor-model candidate rebuilt in Arb.
"""
function arb_bound(sc, sb, beta0, L0, p0, rad, W::TWorker)
    nm = sc.bname
    nm == "const" && return (beta0, L0, true)
    if nm == "boot"
        r0 = crit_eval!(W.cpA, W.tpA.M, W.tpA.K, W.tpA.A2, W.wsA; mode = :LG, beta = beta0, L0)
        r0.g1ok || return (beta0, L0, false)
        g1up = tm_upper(W.cpA.g1)
        bb = Arb(r0.g2lo; prec = PREC) / Arb(g1up; prec = PREC) + Arb(sb.gapAC; prec = PREC)
        return (max(beta0, min(sc.beta, Float64(Arblib.lbound(bb), RoundDown))), L0, true)
    end
    btA = bound_tms(Arb, p0, rad, W.wsA)
    if startswith(nm, "bootTM")
        lev = startswith(nm, "bootTM2:") ? 2 : 1
        nm1 = split(nm, ":"; limit = 2)[2]
        B1 = nm1 == "const" ? beta0 : (i = findfirst(x -> x[1] == nm1, btA.betas); i === nothing ? nothing : btA.betas[i][2])
        B1 === nothing && return (beta0, L0, false)
        base = (; M = W.tpA.M, K = W.tpA.K, A2 = W.tpA.A2)
        b = boot_tm(W.cpA, base, B1, btA.gap, W.wsA); b === nothing && return (beta0, L0, false)
        if lev == 2
            b = boot_tm(W.cpA, base, b, btA.gap, W.wsA); b === nothing && return (beta0, L0, false)
        end
        return (b, L0, true)
    end
    if sc.mode === :LG
        i = findfirst(x -> x[1] == nm, btA.betas); i === nothing && return (beta0, L0, false)
        return (btA.betas[i][2], L0, true)
    else
        i = findfirst(x -> x[1] == nm, btA.L0s); i === nothing && return (beta0, L0, false)
        return (beta0, btA.L0s[i][2], true)
    end
end

function process_transport(path, lo, hi, W::TWorker, pilot::Bool)
    t0 = time()
    rec = Dict{String,Any}("path" => path, "lo" => lo, "hi" => hi)
    if is_shell_box(lo, hi)
        rec["status"] = "SHELL"; rec["t"] = time() - t0
        return rec, 0
    end
    ex = exclusion_t(lo, hi; ws = W.wsA)
    if ex !== :none
        rec["status"] = String(ex); rec["t"] = time() - t0
        return rec, 0
    end
    w = hi .- lo
    if maximum(w) < MINW_T
        rec["status"] = "FAIL"; rec["why"] = "min size"; rec["t"] = time() - t0
        return rec, 0
    end
    p0, rad = box_center_radius(lo, hi)
    cls = classify_point(p0; cube = Float64(CUBE), eta = Float64(RAT_ETA), s0 = Float64(RAT_S1))   # thin band excluded
    if cls !== :R
        ctree = compound(lo, hi, CDEPTH, W.wsA)
        if ctree !== nothing
            rec["status"] = "COMPOUND"; rec["tree"] = ctree; rec["t"] = time() - t0
            return rec, 0
        end
        cls === :nonconvex && return nothing, argmax(w)
    end
    maximum(w) > WTRY_T && return nothing, argmax(w)
    sb = scalar_bounds(p0, rad; ws = W.wsA)
    L0 = max(sb.L_rect, sb.L_hersch)
    beta0 = max(sb.beta_rect_ks, addlo(L0, sb.gapAC))
    # degree: point κ at n = 8 (first diagonal choice), compare with n = 12 when the margin is small
    sh0 = choose_shift(p0)
    k8 = tpoint_kappa(transport_refs(Float64, tsolve(p0, 8; shift = sh0)), beta0, L0)
    n = 8
    if k8 >= 2.0
        k12 = tpoint_kappa(transport_refs(Float64, tsolve(p0, 12; shift = sh0)), beta0, L0)
        n = (2.5 - k8) >= 0.75 * (2.5 - k12) ? 8 : 12
        k8 = min(k8, k12)
    end
    if k8 >= 2.5
        return maximum(w) < 4MINW_T ? (merge(rec, Dict("status" => "FAIL", "why" => "point kappa", "kpt" => k8, "t" => time() - t0)), 0) : (nothing, argmax(w))
    end
    # Float64 screen: preferred diagonal, then the other one; first-order variant, then plain transport;
    # lower-bound candidates in order of observed usage; the first certified combination is kept.
    bt = bound_tms(Float64, p0, rad, W.wsF)
    best = nothing; bestk = Inf; besth = 0; dimsug = 0
    for (sh, vars) in ((sh0, (true,)), (1 - sh0, (true,)), (sh0, (false,)), (1 - sh0, (false,)))
        sc = screen(p0, rad, n, sh, sb, bt, W; variants = vars)
        dimsug == 0 && (dimsug = sc.dim)
        if sc.ok
            best = sc; bestk = sc.kappa; besth = sh; break
        end
    end
    if best === nothing
        # influence-guided split, with a progress guarantee: never split a side narrower than 1/4 of the widest
        d = (dimsug == 0 || w[dimsug] < 0.25 * maximum(w)) ? argmax(w) : dimsug
        return nothing, d
    end
    sc = best; sh = besth
    rec["cert"] = "transport1"; rec["fo"] = sc.fo; rec["n"] = n; rec["shift"] = sh; rec["mode"] = String(sc.mode); rec["bound"] = sc.bname
    rec["beta"] = sc.beta isa TM ? NaN : sc.beta; rec["L0"] = sc.L0 isa TM ? NaN : sc.L0
    rec["kappa_f64"] = bestk; rec["tF"] = time() - t0
    if pilot
        rec["status"] = "pilot_pass"; rec["t"] = time() - t0
        return rec, 0
    end
    # ---- rigorous (Arb): same Float64 choices (sol), reference integrals, Taylor models and bounds in ball arithmetic
    RA = transport1_refs(Arb, sc.sol)
    bA = transport1_base!(W.tpA, RA, p0, rad, W.wsA, sh)
    bA.ok || return nothing, argmax(w)
    betaA, L0A, okb = arb_bound(sc, sb, beta0, L0, p0, rad, W)
    okb || return nothing, argmax(w)
    resA = crit_eval!(W.cpA, W.tpA.M, W.tpA.K, W.tpA.A2, W.wsA; mode = sc.mode, beta = betaA, L0 = L0A)
    if !resA.ok
        rec["arb_fail"] = true
        return nothing, argmax(w)
    end
    betaA isa TM || (rec["beta"] = betaA)
    rec["status"] = "certified"; rec["which"] = String(resA.which)
    rec["dlo"] = bA.dlo; rec["a11up"] = resA.a11up; rec["g1lo"] = resA.g1lo; rec["g2lo"] = resA.g2lo
    rec["ups"] = Dict(String(k) => v for (k, v) in pairs(resA.ups))
    # certified bound κ on ξ over the box: first value of a short ladder above the Float64 estimate that the Arb
    # certificate accepts (the criterion is monotone in κ; κ = 5/2 itself was certified above)
    # (changed 2026-10-01 after the run: logs the upward-rounded Float64 of the rational actually tested; the frozen run log
    #  records the Float64 rung kk itself, which lies within 1e-12 of that rational — README §10.2)
    kc = 2.5
    for kk in (bestk + 0.002, bestk + 0.01, bestk + 0.04)
        kk < 2.5 || break
        kr = rationalize(BigInt, kk; tol = 1e-12)
        rk = crit_eval!(W.cpA, W.tpA.M, W.tpA.K, W.tpA.A2, W.wsA; mode = sc.mode, beta = betaA, L0 = L0A, kappa = kr)
        if rk.ok; kc = Float64(BigFloat(kr; precision = 256), RoundUp); break; end
    end
    rec["kappa_cert"] = kc
    rec["t"] = time() - t0
    return rec, 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    logfile = ARGS[1]
    pilot = "--pilot" in ARGS
    maxb = (i = findfirst(==("--maxboxes"), ARGS)) === nothing ? typemax(Int) : parse(Int, ARGS[i+1])
    roots = root_boxes_t()
    if (i = findfirst(==("--region"), ARGS)) !== nothing
        v = parse.(Float64, split(ARGS[i+1], ","))
        roots = root_boxes_t([v[1], v[3], v[5], v[7]], [v[2], v[4], v[6], v[8]])
    end
    run_cover(logfile; pilot, maxboxes = maxb, roots, procfun = process_transport, mkworker = TWorker)
end
