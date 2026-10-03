# tests.jl — unit tests (BAND_SPEC §5 (c) and more).  Author: Anonymous.
# usage: julia -t 1 tests.jl [nboxes]          (exit code ≠ 0 on any failure)
#
#  T1  containment, exact rationals, at random points of random boxes (both charts):
#        h_Q(x) ≤ h_out(x);  a_Q(x) ≤ a_in(x) < b_in(x) ≤ b_Q(x);  the FLOAT-node Ritz profile actually used by the
#        certificate (inward rounded, P1 interpolated) lies inside [a_in, b_in];  h_Q(x) ≤ H_k on every Sturm cell;
#        proj_x Q ⊂ [0, 1].
#  T2  exclusion soundness at random points of excluded boxes (exact).
#  T3  Iints_arb encloses BigFloat Gauss-Legendre values (random ha, hb incl. ρ ≈ 1/4, Δ = 0, tiny Δ).
#  T4  Sturm check on a constant potential (exact μ1 = V + π²) and gram_arb vs float matrices.
#  T5  tiling checker (replay.jl's) on random trees, and on corrupted trees.
#  T6  Julia's mixed Float64/Rational comparisons are exact (used to place breakpoints in cells).
#  T7  EXACT-DYADIC ENCLOSURE tests (added after the 2026-10-01 review): element_forms, Iints_arb and the full gram_arb
#      assembly (random meshes AND real certificate data) against exact Rational{BigInt} references (polynomial parts,
#      constant-h and end elements) and 2048-bit references for the log/π parts; STRICT containment is required (no
#      tolerance).  The pre-review assembly (Float64 nodal products) is kept as gram_arb_v1 and must FAIL this test.
#  (T4 is a float-vs-Arb sanity check only; it is not the rigour test.)
include(joinpath(@__DIR__, "cert.jl"))
include(joinpath(@__DIR__, "tiling.jl"))
using Random, Printf
Random.seed!(20261001)
nbox = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 400
nfail = 0
fail(msg) = (global nfail += 1; println("FAIL: ", msg))
randq(lo::QQ, hi::QQ) = lo + (hi - lo) * QQ(rand(0:2^20), 2^20)

"random box: a random bisection path of depth 0..28 from the root (the boxes the cover produces)"
function random_box(chart)
    path = join(pathchar(rand(1:4), rand(Bool)) for _ = 1:rand(0:28))
    return box_from_path(chart, path), path
end

function T1(nbox)
    ntest = 0
    for chart in (1, 2), _ = 1:nbox
        B, path = random_box(chart)
        bd = try meshes(B, DEFAULT_PARAMS) catch; continue end
        H = sturm_cells(bd.xs, PWLA(bd.hout))
        a, b = inner_nodes(bd)
        xr = bd.xr; xs = bd.xs
        for _ = 1:6
            p = ntuple(k -> randq(B.lo[k], B.hi[k]), 4)
            V = verts(chart, p)
            all(c -> c > 0, cvals(V)) || continue                   # certificate claims are for convex Q only
            minimum(v[1] for v in V) >= 0 && maximum(v[1] for v in V) <= 1 || fail("proj_x Q not in [0,1] $path")
            xsamp = vcat([QQ(rand(1:2^20 - 1), 2^20) for _ = 1:20], bd.hout.X, bd.bin.X, QQ.(xr[2:4]), QQ.(xs[2:4]))
            for x in xsamp
                (0 < x < 1) || continue
                sl = slice_exact(V, x); sl === nothing && (fail("empty slice"); continue)
                lo, hi = sl
                ntest += 1
                hi - lo <= pwl_eval(bd.hout, x) || fail("h_Q > h_out at x=$x, chart $chart, path $path")
                ai = pwl_eval(bd.ain, x); bi = pwl_eval(bd.bin, x)
                (lo <= ai && ai < bi && bi <= hi) || fail("Q_in not in Q at x=$x chart $chart path $path: $lo $ai $bi $hi")
                # float-node Ritz profile actually used
                j = clamp(searchsortedlast(xr, x), 1, length(xr) - 1)
                xa = QQ(xr[j]); xb = QQ(xr[j+1]); s = (x - xa) / (xb - xa)
                aint = QQ(a[j]) + (QQ(a[j+1]) - QQ(a[j])) * s; bint = QQ(b[j]) + (QQ(b[j+1]) - QQ(b[j])) * s
                (ai <= aint && aint < bint && bint <= bi) || fail("Ritz profile not inside Q_in at x=$x chart $chart path $path")
                # Sturm cell bound
                k = clamp(searchsortedlast(xs, x), 1, length(xs) - 1)
                hi - lo <= QQ(H[k]) || fail("h_Q > H_k at x=$x chart $chart path $path")
                pwl_eval(bd.hout, x) <= QQ(H[k]) || fail("h_out > H_k at x=$x")
            end
        end
    end
    println("T1 containment: $ntest exact slice comparisons, failures so far $nfail")
end

function T2(nbox)
    counts = Dict{Symbol,Int}(); ntest = 0
    for chart in (1, 2), _ = 1:8nbox
        B, path = random_box(chart)
        ex = exclusion(B); ex === nothing && continue
        st = ex[1]; counts[st] = get(counts, st, 0) + 1
        for _ = 1:10
            p = ntuple(k -> randq(B.lo[k], B.hi[k]), 4)
            V = verts(chart, p); c = cvals(V); ntest += 1
            if st == :OUTSIDE_CHART
                valid_point(chart, p) && fail("OUTSIDE_CHART box contains a valid point $p ($path)")
            elseif st == :LEMMA_A
                any(c[i] <= ETA * c[mod1(i + 2, 4)] for i = 1:4) || fail("LEMMA_A violated at $p ($path)")
            elseif st == :LEMMA_B
                valid_point(chart, p) && width2_exact(V) > S_LO^2 && fail("LEMMA_B violated at $p ($path)")
            elseif st == :OUT_OF_BAND
                valid_point(chart, p) && width2_exact(V) <= S_HI^2 && fail("OUT_OF_BAND violated at $p ($path)")
            end
        end
    end
    println("T2 exclusions: $ntest point checks, statuses $counts, failures so far $nfail")
end

function T3()
    setprecision(BigFloat, 200)
    # 60-point Gauss-Legendre in BigFloat on [0,1] for smooth integrands (exact up to ~1e-40 for these rational ones)
    function gl(n)
        x = zeros(BigFloat, n); w = zeros(BigFloat, n)
        for i = 1:n
            z = BigFloat(cos(pi * (i - 0.25) / (n + 0.5)))
            for _ = 1:100
                p1 = one(BigFloat); p2 = zero(BigFloat)
                for j = 1:n; p3 = p2; p2 = p1; p1 = ((2j - 1) * z * p2 - (j - 1) * p3) / j; end
                pp = n * (z * p1 - p2) / (z^2 - 1); z1 = z; z = z1 - p1 / pp
                abs(z - z1) < BigFloat(10)^-60 && break
            end
            p1 = one(BigFloat); p2 = zero(BigFloat)
            for j = 1:n; p3 = p2; p2 = p1; p1 = ((2j - 1) * z * p2 - (j - 1) * p3) / j; end
            pp = n * (z * p1 - p2) / (z^2 - 1)
            x[i] = (1 - z) / 2; w[i] = 1 / ((1 - z^2) * pp^2)
        end
        return x, w
    end
    X, W = gl(80)
    nt = 0
    cases = [(1.0, 1.0), (1.0, 1.0 + 1e-14), (0.3, 0.5), (1.0, 5 / 3 - 1e-9), (1.0, 5 / 3 + 1e-9), (0.02, 0.025), (1e-3, 2e-3)]
    for _ = 1:200; push!(cases, (exp(randn()), exp(randn()))); end
    for (ha, hb) in cases
        I = Iints_arb(Arb(ha), Arb(hb))
        # reference: for strongly varying h use the exact log formula in BigFloat; else Gauss
        Ha = BigFloat(ha); Hb = BigFloat(hb)
        if abs(Hb - Ha) / (Ha + Hb) > 0.05
            Δ = Hb - Ha; lg = log(Hb / Ha)
            r11 = (Δ * (Hb - 3Ha) / 2 + Ha^2 * lg) / Δ^3; r00 = (Δ * (Ha - 3Hb) / 2 + Hb^2 * lg) / Δ^3
            r01 = (Δ - Ha * lg) / Δ^2 - r11
        else
            r00 = sum(W[i] * (1 - X[i])^2 / (Ha * (1 - X[i]) + Hb * X[i]) for i in eachindex(X))
            r01 = sum(W[i] * X[i] * (1 - X[i]) / (Ha * (1 - X[i]) + Hb * X[i]) for i in eachindex(X))
            r11 = sum(W[i] * X[i]^2 / (Ha * (1 - X[i]) + Hb * X[i]) for i in eachindex(X))
        end
        for (Ik, r) in zip(I, (r00, r01, r11))
            nt += 1
            (Arblib.lbound(Ik) <= r + BigFloat(10)^-45 && r - BigFloat(10)^-45 <= Arblib.ubound(Ik)) ||
                fail("Iints_arb($ha,$hb) = $Ik does not contain $r")
            Float64(Arblib.radref(Ik)) < 1e-25 * abs(Float64(r)) || fail("Iints_arb($ha,$hb) too wide: $Ik")
        end
    end
    println("T3 element integrals: $nt enclosure checks, failures so far $nfail")
end

function T4()
    # constant potential V on [0,1]: μ1 = V + π² exactly; sturm_ok must accept E < μ1 and reject E > μ1 (mesh-free)
    for h in (0.05, 0.2, 0.5, 1.0)
        xs = collect(range(0.0, 1.0, length = 401)); H = fill(h, 400)
        μ = pi^2 / h^2 + pi^2
        sturm_ok(xs, H, μ * (1 - 1e-9)) == :ok || fail("sturm rejects E < μ1 (h=$h)")
        sturm_ok(xs, H, μ * (1 + 1e-9)) == :ok && fail("sturm accepts E > μ1 (h=$h)")
    end
    # gram_arb vs float assembly on random boxes
    for chart in (1, 2), _ = 1:20
        B, _ = random_box(chart)
        μ1, R2, bd, F = float_eval(B)
        a, b = inner_nodes(bd)
        G = gram_arb(bd.xr, a, b, F); G === nothing && continue
        N, D = G
        Kd, Ke, Md, Me = ritz_matrices_f(bd.xr, a, b)
        for i = 1:2, j = 1:2
            fi = F[2:end-1, i]; fj = F[2:end-1, j]
            nf = dot(fi, tri_mul(Kd, Ke, fj)); df = dot(fi, tri_mul(Md, Me, fj))
            abs(Float64(N[i, j]) - nf) <= 1e-9 * abs(nf) + 1e-9 || fail("gram N mismatch $(Float64(N[i,j])) vs $nf")
            abs(Float64(D[i, j]) - df) <= 1e-9 * abs(df) + 1e-12 || fail("gram D mismatch")
        end
    end
    println("T4 Sturm/Gram sanity done, failures so far $nfail")
end

function T5()
    for chart in (1, 2), trial = 1:20
        # random tree
        leaves = String[]; stack = [""]
        while !isempty(stack)
            p = pop!(stack)
            if length(p) < 10 && rand() < 0.6
                d = rand(1:4); push!(stack, p * pathchar(d, false)); push!(stack, p * pathchar(d, true))
            else
                push!(leaves, p)
            end
        end
        ok, msg = check_tiling(chart, leaves)
        ok || fail("tiling of a valid tree rejected: $msg")
        if length(leaves) > 1
            bad = copy(leaves); deleteat!(bad, rand(1:length(bad)))
            check_tiling(chart, bad)[1] && fail("tiling with a missing leaf accepted")
            bad = vcat(leaves, [leaves[rand(1:end)] * "a"])
            check_tiling(chart, bad)[1] && fail("overlapping leaf accepted")
            bad = copy(leaves); bad[end] = bad[end] * "b"
            check_tiling(chart, bad)[1] && fail("half-split leaf accepted")
        end
    end
    println("T5 tiling checker done, failures so far $nfail")
end

function T6()
    # sturm_cells / eval_arb locate breakpoints by Float64-vs-Rational comparisons: they must be exact
    for _ = 1:10000
        a = rand(); b = nextfloat(a); m = (QQ(a) + QQ(b)) / 2
        (a < m && m < b && !(m < a) && !(b < m) && searchsortedlast([a, b], m) == 1) || fail("inexact Float64/Rational comparison at $a")
        q = QQ(rand(1:10^9), rand(1:10^9))
        ((Float64(q) < q) == (QQ(Float64(q)) < q)) || fail("inexact comparison for $q")
        ((Float64(q) == q) == (QQ(Float64(q)) == q)) || fail("inexact equality for $q")
    end
    println("T6 exact Float64/Rational comparisons done, failures so far $nfail")
end

# ---------------- T7: exact-dyadic enclosure ----------------
const BPREC = 2048
"strict containment of an exact rational in a ball"
inball(X::Arb, q::QQ) = QQ(BigFloat(Arblib.lbound(X); precision = BPREC)) <= q <= QQ(BigFloat(Arblib.ubound(X); precision = BPREC))
"strict containment of the interval [r − δ, r + δ] (r a 2048-bit reference with error ≤ δ) in a ball"
function inball(X::Arb, r::BigFloat)
    setprecision(BigFloat, BPREC) do
        δ = BigFloat(2)^-1900 * (1 + abs(r))
        lb = BigFloat(Arblib.lbound(X)); ub = BigFloat(Arblib.ubound(X))
        lb <= r - δ && r + δ <= ub
    end
end
"exact / 2048-bit reference of one element: (N_e(f,g)::BigFloat, D_e(f,g)::QQ, I00, I01, I11 refs)"
function ref_element(xa, xb, aa, ab, ba, bb, fa, fb, ga, gb, isfirst, islast)
    setprecision(BigFloat, BPREC) do
        L = QQ(xb) - QQ(xa); ha = QQ(ba) - QQ(aa); hb = QQ(bb) - QQ(ab)
        ap = (QQ(ab) - QQ(aa)) / L; bp = (QQ(bb) - QQ(ba)) / L; hp = bp - ap
        Fa, Fb, Ga, Gb = QQ(fa), QQ(fb), QQ(ga), QQ(gb)
        piB = BigFloat(pi); C6 = inv(BigFloat(6)) + inv(4 * piB^2)
        ce = QQ(1, 2) + ap * bp / 2 + hp^2 * C6
        K1 = (Fb - Fa) * (Gb - Ga) * (ha + hb) / (4 * L)                 # exact
        K2 = hp / 4 * (Fb * Gb - Fa * Ga)                                 # exact
        Dq = L / 24 * (Fa * Ga * (3ha + hb) + (Fa * Gb + Fb * Ga) * (ha + hb) + Fb * Gb * (ha + 3hb))   # exact
        if isfirst
            S = L * Fb * Gb / (2hb); I = (nothing, nothing, nothing)
        elseif islast
            S = L * Fa * Ga / (2ha); I = (nothing, nothing, nothing)
        elseif ha == hb
            I = (1 / (3ha), 1 / (6ha), 1 / (3ha))                         # exact
            S = L * (Fa * Ga * I[1] + (Fa * Gb + Fb * Ga) * I[2] + Fb * Gb * I[3])
        else
            Ha = BigFloat(ha); Hb = BigFloat(hb); Δ = Hb - Ha; lg = log(Hb / Ha)
            I11 = (Δ * (Hb - 3Ha) / 2 + Ha^2 * lg) / Δ^3
            I00 = (Δ * (Ha - 3Hb) / 2 + Hb^2 * lg) / Δ^3
            I01 = (Δ - Ha * lg) / Δ^2 - I11
            I = (I00, I01, I11)
            S = BigFloat(L) * (BigFloat(Fa * Ga) * I00 + BigFloat(Fa * Gb + Fb * Ga) * I01 + BigFloat(Fb * Gb) * I11)
        end
        N = BigFloat(K1 + K2) + piB^2 * ce * S
        return N, Dq, I
    end
end
"the pre-review (DEFECTIVE) assembly: nodal products in Float64 before conversion — must fail T7"
function gram_arb_v1(xs, a, b, F)
    n = length(xs); m = size(F, 2)
    N = [Arb(0) for _ = 1:m, _ = 1:m]; D = [Arb(0) for _ = 1:m, _ = 1:m]
    for e = 1:n-1
        loc = element_forms(xs[e], xs[e+1], a[e], a[e+1], b[e], b[e+1], e == 1, e == n - 1)
        Kaa, Kab, Kbb, Maa, Mab, Mbb = loc
        for i = 1:m, j = i:m
            fa = F[e, i]; fb = F[e+1, i]; ga = F[e, j]; gb = F[e+1, j]
            N[i, j] += Kaa * (fa * ga) + Kab * (fa * gb + fb * ga) + Kbb * (fb * gb)
            D[i, j] += Maa * (fa * ga) + Mab * (fa * gb + fb * ga) + Mbb * (fb * gb)
        end
    end
    return N, D
end
"reference Gram matrices of a mesh (sum of element references)"
function ref_gram(xs, a, b, F)
    n = length(xs); m = size(F, 2)
    N = [BigFloat(0; precision = BPREC) for _ = 1:m, _ = 1:m]; D = [QQ(0) for _ = 1:m, _ = 1:m]
    for e = 1:n-1, i = 1:m, j = i:m
        Ne, De, _ = ref_element(xs[e], xs[e+1], a[e], a[e+1], b[e], b[e+1], F[e, i], F[e+1, i], F[e, j], F[e+1, j], e == 1, e == n - 1)
        setprecision(BigFloat, BPREC) do; N[i, j] += Ne; end
        D[i, j] += De
    end
    return N, D
end
function random_mesh(n)
    xs = sort(vcat(0.0, rand(n - 2), 1.0))
    a = vcat(0.0, -0.3 .* rand(n - 2), 0.0); b = vcat(0.0, 0.05 .+ 0.4 .* rand(n - 2), 0.0)
    if rand() < 0.3                                   # some elements with exactly constant h
        k = rand(2:n-2); b[k+1] = b[k] + (a[k+1] - a[k])
    end
    F = hcat(vcat(0.0, randn(n - 2), 0.0), vcat(0.0, randn(n - 2), 0.0))
    return xs, a, b, F
end
function T7()
    nt = 0; nfail0 = nfail
    # (a) single elements: local forms combined with random nodal values, all element types
    for trial = 1:400
        kind = rand(1:4)          # 1 interior generic, 2 interior constant h, 3 first, 4 last
        xa = rand(); xb = xa + 10.0^(-rand(1:4)) * rand()
        xb > xa || continue
        aa, ab = -0.3rand(), -0.3rand(); ba, bb = 0.05 + 0.4rand(), 0.05 + 0.4rand()
        kind == 2 && (bb = ba + (ab - aa))
        kind == 3 && (aa = 0.0; ba = 0.0)
        kind == 4 && (ab = 0.0; bb = 0.0)
        (kind == 2 && QQ(bb) - QQ(ab) != QQ(ba) - QQ(aa)) && continue     # rounding broke constancy: skip
        fa, fb, ga, gb = randn(), randn(), randn(), randn()
        kind == 3 && (fa = 0.0; ga = 0.0); kind == 4 && (fb = 0.0; gb = 0.0)
        loc = element_forms(xa, xb, aa, ab, ba, bb, kind == 3, kind == 4)
        loc === nothing && (fail("element_forms returned nothing on a valid element"); continue)
        Kaa, Kab, Kbb, Maa, Mab, Mbb = loc
        FA, FB, GA, GB = Arb(fa), Arb(fb), Arb(ga), Arb(gb)
        Nc = Kaa * (FA * GA) + Kab * (FA * GB + FB * GA) + Kbb * (FB * GB)
        Dc = Maa * (FA * GA) + Mab * (FA * GB + FB * GA) + Mbb * (FB * GB)
        Nr, Dr, _ = ref_element(xa, xb, aa, ab, ba, bb, fa, fb, ga, gb, kind == 3, kind == 4)
        nt += 2
        inball(Nc, Nr) || fail("T7a: N_e ball does not contain the reference (kind $kind)")
        inball(Dc, Dr) || fail("T7a: D_e ball does not contain the exact rational (kind $kind)")
    end
    # (b) Iints_arb on both branches, near the switch ρ = 1/4, tiny and zero Δ
    for trial = 1:600
        ha = exp(randn())
        # target r = |hb − ha|/(ha + hb) (the code's ρ = |ε|/2): straddling the switch r = 1/4, tiny, and generic
        r = trial % 3 == 0 ? 0.25 * (1 + 10.0^(-rand(3:12)) * randn()) : (trial % 3 == 1 ? 10.0^(-rand(0:12)) * rand() : 0.999 * rand())
        hb = rand(Bool) ? ha * (1 + r) / (1 - r) : ha * (1 - r) / (1 + r)
        trial % 50 == 0 && (hb = ha)
        (hb > 0 && isfinite(hb)) || continue
        I = Iints_arb(Arb(ha), Arb(hb))
        _, _, R = ref_element(0.0, 1.0, 0.0, 0.0, ha, hb, 1.0, 1.0, 1.0, 1.0, false, false)
        for k = 1:3
            nt += 1
            r = R[k]
            ok = r isa QQ ? inball(I[k], r) : inball(I[k], r)
            ok || fail("T7b: Iints_arb($ha, $hb)[$k] does not contain the reference")
        end
    end
    # (c) full gram_arb on random meshes (incl. constant-h elements), and the old assembly must fail
    nold = 0; ncases = 0
    for trial = 1:60
        xs, a, b, F = random_mesh(rand(4:12))
        G = gram_arb(xs, a, b, F); G === nothing && continue
        Nr, Dr = ref_gram(xs, a, b, F); ncases += 1
        for i = 1:2, j = i:2
            nt += 2
            inball(G[1][i, j], Nr[i, j]) || fail("T7c: gram_arb N[$i,$j] misses the reference")
            inball(G[2][i, j], Dr[i, j]) || fail("T7c: gram_arb D[$i,$j] misses the exact rational")
        end
        Go = gram_arb_v1(xs, a, b, F)
        nold += !all(inball(Go[1][i, j], Nr[i, j]) && inball(Go[2][i, j], Dr[i, j]) for i = 1:2 for j = i:2)
    end
    # (d) real certificate data: gram_arb exactly as certify() calls it, on random cover-shaped boxes
    nreal = 0
    for chart in (1, 2), _ = 1:6
        B, _ = random_box(chart)
        μ1, R2, bd, F = float_eval(B)
        a, b = inner_nodes(bd)
        G = gram_arb(bd.xr, a, b, F); G === nothing && continue
        Nr, Dr = ref_gram(bd.xr, a, b, F); nreal += 1
        for i = 1:2, j = i:2
            nt += 2
            inball(G[1][i, j], Nr[i, j]) || fail("T7d: gram_arb N[$i,$j] misses the reference on real data (chart $chart)")
            inball(G[2][i, j], Dr[i, j]) || fail("T7d: gram_arb D[$i,$j] misses the exact rational on real data (chart $chart)")
        end
        Go = gram_arb_v1(bd.xr, a, b, F)
        nold += !all(inball(Go[1][i, j], Nr[i, j]) && inball(Go[2][i, j], Dr[i, j]) for i = 1:2 for j = i:2)
        ncases += 1
    end
    nold >= ncases ÷ 2 || fail("T7: the defective Float64-product assembly was NOT detected ($nold/$ncases cases)")
    println("T7 exact-dyadic enclosure: $nt strict containment checks ($nreal real certificate meshes); the defective ",
            "pre-review assembly fails containment in $nold/$ncases meshes (as it must); failures so far $nfail")
end

t0 = time()
T1(nbox); T2(nbox); T3(); T4(); T5(); T6(); T7()
@printf("ALL TESTS: %s (%d failures, %.1f s)\n", nfail == 0 ? "PASS" : "FAIL", nfail, time() - t0)
exit(nfail == 0 ? 0 : 1)
