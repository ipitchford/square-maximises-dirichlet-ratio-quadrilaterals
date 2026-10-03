# tests_t.jl — self-tests of the affine-transport certificate and semantic checks of a transport cover log.
#
# usage: julia -t 1 tests_t.jl [logfile] [--quick]       (logfile default runs/transport_cert.jsonl)
#
# 0. Exact-rational two-triangle Piola conformity: unequal determinants d1 ≠ d2, both affine maps agree on the
#    diagonal, equal physical normal traces of w_K = A_K(d_K φ̂ e + curl ψ̂)/d_K, and div_x w = −u (rational arithmetic).
# 1. Point identity: at δ = 0 the certified 3×3 pencil (constant terms of the order-4 Taylor models; D(p0) = 1) has
#    the Rayleigh–Ritz eigenvalues of fsolve (independent physical-coordinate quadrature; validation_points.csv,
#    which also carries the FEM comparison of README §4), for both diagonals; Arb constant terms ⊇ Float64.
# 2. Taylor-model enclosure: at random points of boxes (thin and fat shapes, r = 0.01 and 0.04, first-order and plain
#    transport, both diagonals) the directly assembled Float64 pencil D·K, D·M, D·A2 (full-basis transported matrices
#    of tsolve1, a separate Float64 code path) lies in the Arb Taylor-model enclosure.
# 3. D4 invariance of the transported Rayleigh–Ritz values.
# 4. Pointwise lower bounds: every β candidate (incl. the bootstrap bootTM) ≤ θ2(p) and every L0 candidate ≤ θ1(p)
#    at random p (θ = Rayleigh–Ritz upper bounds, n = 12): necessary conditions for validity.
# 5. Log semantics: random points of the a-priori box ∩ F are located in the logged tree; excluded leaves must satisfy
#    their condition at the point (Float64), certified leaves must have θ2/θ1 ≤ kappa_cert (+5e-3, the Rayleigh–Ritz error of thin shapes); plus random points
#    in random certified / excluded leaves.
include("replay_t.jl")
using Random, Printf, LinearAlgebra
const NP = Ref(0); const NF = Ref(0)
check(c::Bool, msg) = c ? (NP[] += 1; true) : (NF[] += 1; println("FAIL: ", msg); false)
quick = "--quick" in ARGS
logf = (a = filter(x -> !startswith(x, "--"), ARGS); isempty(a) ? joinpath(@__DIR__, "runs", "transport_cert.jsonl") : a[1])
W = TWorker(); rng = MersenneTwister(20261001)

function tm_eval_arb(f::TM{Arb}, t)
    tab = TAB[]; s = Arb(0; prec = PREC); ta = [x isa Arb ? x : Arb(x; prec = PREC) for x in t]
    for k = 1:tab.nc
        m = tab.mono[k]
        s += f.c[k] * (ta[1]^m[1] * ta[2]^m[2] * ta[3]^m[3] * ta[4]^m[4])
    end
    s + Arb(f.r[1]; prec = PREC) * UNIT
end
inball(b::Arb, v, tol) = Float64(Arblib.lbound(b), RoundDown) - tol <= v <= Float64(Arblib.ubound(b), RoundUp) + tol
midf(b::Arb) = Float64(Arblib.midref(b))
pen3(tp) = (K = zeros(3, 3); M = zeros(3, 3);
            for (q, (i, j)) in enumerate(KL); K[i, j] = K[j, i] = tp.K[q].c[1]; M[i, j] = M[j, i] = tp.M[q].c[1]; end; (K, M))

# ---------------------------------------------------------------- 0
println("== 0. exact-rational two-triangle Piola conformity (unequal determinants)")
const QQ = Rational{BigInt}
# polynomials in (x, y) with exact rational coefficients: Dict((i, j) => c) for c·x^i·y^j
peval(P, x, y) = sum(c * x^i * y^j for ((i, j), c) in P; init = zero(QQ))
pdx(P) = Dict((i - 1, j) => c * i for ((i, j), c) in P if i > 0)
pdy(P) = Dict((i, j - 1) => c * j for ((i, j), c) in P if j > 0)
vertsQ(p) = [QQ(V0[j, k]) + sum(QQ(DV[j, k, i]) * p[i] for i = 1:4) for j = 1:4, k = 1:2]   # DV entries are dyadic: exact
rot(v) = [-v[2], v[1]]
let ok = true, ntr = 0
    phi = Dict((0, 0) => QQ(3, 7), (1, 0) => QQ(1), (1, 1) => QQ(-2, 5), (0, 2) => QQ(3), (3, 0) => QQ(1, 3))
    psi = Dict((2, 0) => QQ(1, 2), (0, 1) => QQ(-1), (1, 2) => QQ(2, 9), (2, 1) => QQ(1))
    pref = QQ.([1 // 5, 1 // 4, 1 // 3, -1 // 10]); ptgt = pref .+ QQ.([1 // 50, -1 // 40, 1 // 30, 1 // 20])
    for sh in (0, 1)
        Vr = vertsQ(pref)[vperm(sh), :]; Vt = vertsQ(ptgt)[vperm(sh), :]
        e = Vr[3, :] - Vr[1, :]                                       # reference diagonal v1 -> v3 (τ̂ = φ̂ e)
        A = Matrix{QQ}[]; d = QQ[]; Psi = Function[]
        for (i0, i1, i2, _) in TRIS
            Mr = [Vr[i1, 1]-Vr[i0, 1] Vr[i2, 1]-Vr[i0, 1]; Vr[i1, 2]-Vr[i0, 2] Vr[i2, 2]-Vr[i0, 2]]
            Mt = [Vt[i1, 1]-Vt[i0, 1] Vt[i2, 1]-Vt[i0, 1]; Vt[i1, 2]-Vt[i0, 2] Vt[i2, 2]-Vt[i0, 2]]
            Ak = Mt * inv(Mr); push!(A, Ak); push!(d, Ak[1, 1] * Ak[2, 2] - Ak[1, 2] * Ak[2, 1])
            push!(Psi, xh -> Vt[i0, :] + Ak * (xh - Vr[i0, :]))
        end
        ok &= check(d[1] > 0 && d[2] > 0 && d[1] != d[2], "Piola test: determinants $(d) must be positive and unequal")
        ok &= check(A[1] * e == A[2] * e, "Piola test: the two affine maps disagree on the diagonal direction")
        nphys = rot(Vt[3, :] - Vt[1, :])                                 # normal of the physical diagonal
        what(k, xh) = d[k] * peval(phi, xh...) .* e .+ [peval(pdy(psi), xh...), -peval(pdx(psi), xh...)]   # d_K φ̂ e + curl ψ̂
        wphys(k, xh) = A[k] * what(k, xh) / d[k]                          # contravariant Piola
        for s in QQ.([0, 1 // 3, 1 // 2, 4 // 5, 1])
            xh = Vr[1, :] + s * e
            ok &= check(Psi[1](xh) == Psi[2](xh), "Piola test: maps differ on the diagonal at s=$s")
            w1 = wphys(1, xh); w2 = wphys(2, xh)
            ok &= check(nphys' * w1 == nphys' * w2, "Piola test: normal traces differ at s=$s (shift $sh)")
            ntr += (w1 != w2)                                             # the tangential parts do differ (non-trivial test)
        end
        # divergence: div_x w_K = tr(A_K Dŵ A_K⁻¹)/d_K = div_x̂ ŵ_K / d_K = e·∇φ̂ = −û with û := −e·∇φ̂
        for (k, xh) in ((1, (Vr[1, :] + Vr[2, :] + Vr[3, :]) / 3), (2, (Vr[1, :] + Vr[4, :] + Vr[3, :]) / 3), (1, Vr[2, :] / 2 + Vr[3, :] / 2))
            Dw = [d[k] * e[1] * peval(pdx(phi), xh...) + peval(pdx(pdy(psi)), xh...)  d[k] * e[1] * peval(pdy(phi), xh...) + peval(pdy(pdy(psi)), xh...);
                  d[k] * e[2] * peval(pdx(phi), xh...) - peval(pdx(pdx(psi)), xh...)  d[k] * e[2] * peval(pdy(phi), xh...) - peval(pdy(pdx(psi)), xh...)]
            divx = sum(diag(A[k] * Dw * inv(A[k]))) / d[k]
            uhat = -(e[1] * peval(pdx(phi), xh...) + e[2] * peval(pdy(phi), xh...))
            ok &= check(divx == -uhat, "Piola test: div_x w ≠ −u on triangle $k (shift $sh)")
        end
    end
    check(ntr > 0, "Piola test: tangential traces never differed (test not discriminating)")
    println(ok ? "   exact: d1 ≠ d2, maps agree on the diagonal, normal traces equal, tangential traces differ ($ntr points), div_x w = −u" : "   FAILED")
end

# ---------------------------------------------------------------- 1
println("== 1. point identity at δ = 0 (validation_points.csv), both diagonals")
rows = [split(l, ",") for l in readlines(joinpath(@__DIR__, "validation_points.csv"))[2:end]]
quick && (rows = rows[1:4:end])
let maxrel = 0.0, nA = 0
    for r in rows
        p = parse.(Float64, r[2:5]); n = parse(Int, r[14]); th = parse.(Float64, r[15:17])
        for sh in (0, 1)
            sol = tsolve1(p, n; shift = sh)
            RF = transport1_refs(Float64, sol); transport1_base!(W.tpF, RF, p, fill(1e-3, 4), W.wsF, sh)
            K3, M3 = pen3(W.tpF)
            θ = eigvals(Symmetric(K3), Symmetric(M3))
            rel = maximum(abs.(θ .- th) ./ th); maxrel = max(maxrel, rel)
            check(rel < 1e-6, "point identity p=$p n=$n shift=$sh rel=$rel")
            if n == 8
                RA = transport1_refs(Arb, sol); transport1_base!(W.tpA, RA, p, fill(1e-3, 4), W.wsA, sh)
                for q = 1:6, (fA, fF) in ((W.tpA.K[q], W.tpF.K[q]), (W.tpA.M[q], W.tpF.M[q]), (W.tpA.A2, W.tpF.A2))
                    check(inball(fA.c[1], fF.c[1], 1e-9 * max(1.0, abs(fF.c[1]))), "Arb ⊇ F64 constant term p=$p shift=$sh")
                end
                nA += 1
            end
        end
    end
    @printf("   %d rows × 2 diagonals: max relative deviation of θ1..θ3 from fsolve %.2e; Arb ⊇ Float64 checked on %d\n", length(rows), maxrel, nA)
end

# ---------------------------------------------------------------- 2
println("== 2. Taylor-model enclosure vs direct Float64 assembly at random points")
function direct_pencil(sol, p)
    Kp, Mp, TT, TC, CC = sol.pencil_at(p)
    δ = p .- sol.p0
    C = sol.C .+ sum(δ[i] .* sol.dC[:, :, i] for i = 1:4)
    ψ = sol.dpsi .+ sol.ddpsi * δ
    c1 = C[:, 1]
    A2 = dot(c1, TT * c1) + 2 * dot(c1, TC * ψ) + dot(ψ, CC * ψ)
    V0p = verts(sol.p0)[vperm(sol.shift), :]; Vp = verts(p)[vperm(sol.shift), :]
    D = 1.0
    for (i0, i1, i2, _) in TRIS
        Mk(V) = [V[i1, 1]-V[i0, 1] V[i2, 1]-V[i0, 1]; V[i1, 2]-V[i0, 2] V[i2, 2]-V[i0, 2]]
        D *= det(Mk(Vp) / Mk(V0p))
    end
    return D .* (C' * Kp * C), D .* (C' * Mp * C), D * A2
end
"the transported pencil D·K, D·M, D·A2 at the point p in ball arithmetic from the Arb reference integrals (no TMs)"
function direct_pencil_arb(R::TRefs1{Arb}, p0, p, shift)
    A_(x) = Arb(x; prec = PREC)
    Vp = verts_T(Arb, p)[vperm(shift), :]
    e = [A_(1.0); [A_(p[i]) - A_(p0[i]) for i = 1:4]]
    quad(X, rows, cols) = sum(e[m] * e[n] * X[rows[m], cols[n]] for m = 1:5, n = 1:5)
    rowsk(k) = [5(k - 1) + mu + 1 for mu = 0:4]
    K = [A_(0) for _ = 1:3, _ = 1:3]; M = [A_(0) for _ = 1:3, _ = 1:3]; A2 = A_(0); ds = Arb[]
    for (t, (i0, i1, i2, _)) in enumerate(TRIS)
        Mk = [Vp[i1, 1]-Vp[i0, 1] Vp[i2, 1]-Vp[i0, 1]; Vp[i1, 2]-Vp[i0, 2] Vp[i2, 2]-Vp[i0, 2]]
        A = Mk * R.Minv0[t]
        d = A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]; push!(ds, d)
        adj = [A[2, 2] -A[1, 2]; -A[2, 1] A[1, 1]]; Ga = adj * permutedims(adj); GA = permutedims(A) * A
        for i = 1:3, j = 1:3
            K[i, j] += sum(Ga[a, b] * quad(R.I[t][a, b], rowsk(i), rowsk(j)) for a = 1:2, b = 1:2) / d
            M[i, j] += d * quad(R.M[t], rowsk(i), rowsk(j))
        end
        r5 = collect(1:5)
        A2 += d * sum(GA[a, b] * R.e[a] * R.e[b] for a = 1:2, b = 1:2) * quad(R.TT[t], r5, r5)
        A2 += 2 * sum(GA[a, b] * R.e[a] * quad(R.TC[t], r5, [5(b - 1) + nu for nu = 1:5]) for a = 1:2, b = 1:2)
        A2 += sum(GA[a, b] * quad(R.CC[t], [5(a - 1) + mu for mu = 1:5], [5(b - 1) + nu for nu = 1:5]) for a = 1:2, b = 1:2) / d
    end
    D = ds[1] * ds[2]
    return D .* K, D .* M, D * A2
end
boxes = [([0.3, 0.5, 0.8, 0.5], 0.01), ([0.5, 0.2, 0.4, -0.2], 0.01), ([0.3, 0.2, 0.4, 0.1], 0.04),
         ([0.171, 0.374, 0.092, -0.438], 0.02), ([0.6, 0.1, 0.2, 0.3], 0.01), ([0.25, 0.3, 0.3, 0.0], 0.04),
         ([0.597, 0.15, 0.15, -0.15], 0.005), ([0.35, 0.45, 0.25, 0.15], 0.02)]
for (p0, r) in boxes; @assert minimum(cvals(p0)) > 0.05 "test box centre not convex: $p0"; end
quick && (boxes = boxes[1:3])
let worst = 0.0, npt = 0, f64dev = 0.0
    for (p0, r) in boxes, sh in (choose_shift(p0), 1 - choose_shift(p0)), fo in (true, false)
        rad = fill(r, 4)
        sol = tsolve1(p0, 8; shift = sh)
        fo || (sol = merge(sol, (; dC = zero(sol.dC), ddpsi = zero(sol.ddpsi))))
        RA = transport1_refs(Arb, sol)
        b = transport1_base!(W.tpA, RA, p0, rad, W.wsA, sh)
        b.ok || (println("   (fold on box p0=$p0 r=$r shift=$sh: skipped)"); continue)
        for _ = 1:(quick ? 3 : 6)
            t = 0.999 .* (2 .* rand(rng, 4) .- 1); p = p0 .+ rad .* t
            # the exact Taylor-model variable of the Float64 point p (p0 + rad·t is rounded; |tA| < 1 strictly)
            tA = [(Arb(p[i]; prec = PREC) - Arb(p0[i]; prec = PREC)) / Arb(rad[i]; prec = PREC) for i = 1:4]
            K3, M3, A2 = direct_pencil(sol, p)                 # Float64, from tsolve1's own reference matrices
            KA, MA, A2A = direct_pencil_arb(RA, p0, p, sh)       # ball arithmetic, from the Arb reference integrals
            for (q, (i, j)) in enumerate(KL)
                for (f, v, vA) in ((W.tpA.K[q], K3[i, j], KA[i, j]), (W.tpA.M[q], M3[i, j], MA[i, j]))
                    e = tm_eval_arb(f, tA)
                    # strict: the TM enclosure must overlap the rigorous point value (both are balls)
                    check(Arblib.overlaps(e, vA), "TM enclosure vs Arb point value p0=$p0 r=$r sh=$sh fo=$fo q=$q enc=$(e) point=$(vA)")
                    worst = max(worst, abs(midf(e) - midf(vA)) / (Float64(Arblib.radref(e)) + Float64(Arblib.radref(vA)) + 1e-300))
                    f64dev = max(f64dev, abs(v - midf(vA)) / max(1.0, abs(midf(vA))))
                end
            end
            e = tm_eval_arb(W.tpA.A2, tA)
            check(Arblib.overlaps(e, A2A), "TM enclosure A2 vs Arb point value p0=$p0 r=$r sh=$sh fo=$fo")
            f64dev = max(f64dev, abs(A2 - midf(A2A)) / max(1.0, abs(midf(A2A))))
            npt += 1
        end
    end
    @printf("   %d points: TM enclosure overlaps the Arb point value everywhere required; max |mid(TM) − point| / (sum of radii) = %.3f\n", npt, worst)
    @printf("   Float64 assembly from tsolve1's matrices vs Arb point value: max relative deviation %.2e (Float64 conditioning; informative)\n", f64dev)
end

# ---------------------------------------------------------------- 3
println("== 3. D4 invariance of the transported Rayleigh–Ritz values")
function d4_images(p)
    G = [p]
    for _ = 1:3, q in copy(G), g in (sigma_h, sigma_v, rho_rot)
        x = g(q); any(y -> maximum(abs.(y .- x)) < 1e-15, G) || push!(G, x)
    end
    G
end
let worst = 0.0
    for p in ([0.3, 0.5, 0.8, 0.5], [0.5, 0.2, 0.4, -0.2], [0.3, 0.2, 0.4, 0.1], [0.6, 0.1, 0.2, 0.3])
        imgs = d4_images(p)
        check(length(imgs) == 8, "D4 orbit size $(length(imgs))")
        ths = [tsolve1(q, 8).theta[1:3] for q in imgs]
        spread = maximum(maximum(abs.(t .- ths[1]) ./ ths[1]) for t in ths)
        worst = max(worst, spread)
        check(spread < 1e-7, "D4 spread $spread at p=$p")
    end
    @printf("   max relative spread over the 8 images: %.2e\n", worst)
end

# ---------------------------------------------------------------- 4
println("== 4. pointwise lower bounds below the Rayleigh–Ritz values")
let nchk = 0, minslack = Inf
    for (p0, r) in boxes
        rad = fill(r, 4)
        btA = bound_tms(Arb, p0, rad, W.wsA)
        sh = choose_shift(p0); sol = tsolve1(p0, 8; shift = sh)
        RA = transport1_refs(Arb, sol); b = transport1_base!(W.tpA, RA, p0, rad, W.wsA, sh)
        base = (; M = W.tpA.M, K = W.tpA.K, A2 = W.tpA.A2)
        boots = Tuple{String,Any}[]
        if b.ok
            for (nm, B) in btA.betas
                bb = boot_tm(W.cpA, base, B, btA.gap, W.wsA); bb === nothing || push!(boots, ("bootTM:" * nm, bb))
            end
        end
        for _ = 1:(quick ? 2 : 4)
            t = 2 .* rand(rng, 4) .- 1; p = p0 .+ rad .* t
            cls = classify_point(p; cube = 0.0, eta = 0.0, s0 = 0.0)
            cls === :nonconvex && continue
            th = tsolve1(p, 12).theta
            for (nm, B) in vcat(btA.betas, boots)
                v = midf(tm_eval_arb(B, t)); minslack = min(minslack, (th[2] - v) / th[2])
                check(v <= th[2] * (1 + 1e-9), "beta candidate $nm = $v > θ2 = $(th[2]) at p=$p"); nchk += 1
            end
            for (nm, L) in btA.L0s
                v = midf(tm_eval_arb(L, t)); minslack = min(minslack, (th[1] - v) / th[1])
                check(v <= th[1] * (1 + 1e-9), "L0 candidate $nm = $v > θ1 = $(th[1]) at p=$p"); nchk += 1
            end
        end
    end
    @printf("   %d candidate evaluations; smallest relative slack below θ: %.3e\n", nchk, minslack)
end

# ---------------------------------------------------------------- 5
function logtests(logf, rng, quick)
    println("== 5. log semantics: ", logf)
    recs = Dict{String,Dict{String,Any}}()
    for ln in eachline(logf)
        r = parse_rec(ln); haskey(r, "path") && (recs[r["path"]] = r)
    end
    roots = root_boxes_t()
    function locate(p)
        for (path, lo, hi) in roots
            all(lo .<= p .<= hi) || continue
            while true
                r = get(recs, path, nothing); r === nothing && return nothing
                r["status"] == "split" || return r
                d = Int(r["dim"]); mid = (lo[d] + hi[d]) / 2
                c = children(path, lo, hi, d)
                (path, lo, hi) = p[d] <= mid ? c[1] : c[2]
            end
        end
        return nothing
    end
    w2d(p) = (wd = width_diam(verts(p)); wd[1] / wd[2])
    function semantic(r, p)
        st = r["status"]; c = cvals(p)
        minimum(c) <= 0 && return true                                     # non-convex: excluded by definition
        st == "E1" && return maximum(abs, p) <= 0.049
        st == "E4" && return false
        st == "E2" && return any(c[i] <= 0.035 * c[mod1(i + 2, 4)] for i = 1:4)
        st == "E3" && return w2d(p) <= 0.1 + 1e-12
        st == "E3T" && return w2d(p) <= 0.4 + 1e-12
        st == "COMPOUND" && return maximum(abs, p) <= 0.049 || any(c[i] <= 0.035 * c[mod1(i + 2, 4)] for i = 1:4) || w2d(p) <= 0.4 + 1e-12
        st == "SHELL" && return all(0 .<= p[1:3] .<= 0.225) && abs(p[4]) <= 0.225
        if st == "certified"
            th = tsolve1(p, 8).theta
            return th[2] / th[1] <= r["kappa_cert"] + 5e-3
        end
        return false
    end
    cnt = Dict{String,Int}(); nmiss = 0
    lo0, hi0 = A_PRIORI
    for _ = 1:(quick ? 500 : 4000)
        p = lo0 .+ (hi0 .- lo0) .* rand(rng, 4)
        r = locate(p)
        if r === nothing; nmiss += 1; check(false, "point $p not located in the tree"); continue; end
        cnt[r["status"]] = get(cnt, r["status"], 0) + 1
        check(semantic(r, p), "semantic $(r["status"]) at p=$p (leaf $(r["path"]))")
    end
    @printf("   uniform points in the a-priori box: %s, not located %d\n", string(cnt), nmiss)
    leaves = [r for r in values(recs) if r["status"] != "split"]
    bykind = Dict{String,Vector{Dict{String,Any}}}()
    for r in leaves; push!(get!(bykind, r["status"], Dict{String,Any}[]), r); end
    worstk = -Inf
    for (st, rs) in bykind
        nsamp = min(length(rs), st == "certified" ? (quick ? 100 : 1500) : (quick ? 50 : 300))
        for r in rs[randperm(rng, length(rs))[1:nsamp]]
            p = r["lo"] .+ (r["hi"] .- r["lo"]) .* rand(rng, 4)
            check(semantic(r, p), "leaf semantic $st at p=$p (leaf $(r["path"]))")
            if st == "certified" && minimum(cvals(p)) > 0
                th = tsolve1(p, 8).theta; worstk = max(worstk, th[2] / th[1] - r["kappa_cert"])
            end
        end
        @printf("   %-10s %6d leaves, %5d sampled\n", st, length(rs), nsamp)
    end
    @printf("   certified leaves: max over samples of θ2/θ1 − kappa_cert = %.4f (must be ≤ 5e-3)\n", worstk)
end
isfile(logf) && logtests(logf, rng, quick)
@printf("\ntests_t: %d passed, %d failed\n", NP[], NF[])
println(NF[] == 0 ? "ALL TESTS PASSED" : "SOME TESTS FAILED")
