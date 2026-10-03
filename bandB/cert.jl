# cert.jl — the box certificate of BAND_SPEC §3 (rigorous layer in Arb ball arithmetic).  Author: Anonymous.
#
# For a box B:  λ1(Q) > E  (Sturm shooting on a step minorant of π²/h_out², rigorous)  and
#               λ2(Q) ≤ R2 < (5/2)E  (Rayleigh–Ritz on Q_in with u = f(x) sin(π(y−a)/h), Gram matrices rigorous)
# for every convex Q in B.  Floating point is used only to CHOOSE E, the meshes and the test functions f1, f2.
using Arblib, LinearAlgebra
include(joinpath(@__DIR__, "geom.jl"))
include(joinpath(@__DIR__, "cert_float.jl"))

const PREC = 128
setprecision(Arb, PREC)
const PI_A = Arb(π)
const PI2_A = PI_A^2
const C6_A = inv(Arb(6)) + inv(4 * PI2_A)
const HALF_A = inv(Arb(2))

# exact series coefficients (σ ∈ [−1/2,1/2]):  C_n = ∫σ^n(1/2+σ)²,  B_n = ∫σ^n(1/4−σ²),  M_n = ∫σ^n
const NSER = 160
mom(n) = iseven(n) ? QQ(1, big(2)^n * (n + 1)) : QQ(0)
const SERC_A = [Arb(mom(n) / 4 + mom(n + 1) + mom(n + 2)) for n = 0:NSER]
const SERB_A = [Arb(mom(n) / 4 - mom(n + 2)) for n = 0:NSER]

Base.@kwdef struct Params
    dS::Float64 = 0.002      # Sturm mesh: weight per cell
    pS::Float64 = 3.0
    dxS::Float64 = 0.02
    dR::Float64 = 0.02       # Ritz mesh
    pR::Float64 = 2.0
    dxR::Float64 = 0.008
    Erel::Float64 = 1e-7     # E = (1 − Erel)·(float μ1 of the step minorant)
end
const DEFAULT_PARAMS = Params()

up(x::Arb) = Float64(Arblib.ubound(x), RoundUp)
dn(x::Arb) = Float64(Arblib.lbound(x), RoundDown)
updec(x::Arb) = repr(nextfloat(up(x)))       # decimal string ≥ every point of the ball
dndec(x::Arb) = repr(prevfloat(dn(x)))       # decimal string ≤ every point of the ball

# --- PWL profiles in Arb (exact rational breakpoints; segment located by exact comparison) ---
struct PWLA
    X::Vector{QQ}
    XA::Vector{Arb}
    YA::Vector{Arb}
end
PWLA(p::PWL) = PWLA(p.X, [Arb(x) for x in p.X], [Arb(y) for y in p.Y])
"ball enclosing the exact value p(x) at the float x"
function eval_arb(p::PWLA, x::Float64)
    xq = QQ(x)
    j = clamp(searchsortedlast(p.X, xq), 1, length(p.X) - 1)
    return p.YA[j] + (p.YA[j+1] - p.YA[j]) * (Arb(x) - p.XA[j]) / (p.XA[j+1] - p.XA[j])
end
evalf(p::PWL, x::Float64) = (Xf = p.X; j = clamp(searchsortedlast(Xf, QQ(x)), 1, length(Xf) - 1);
    Float64(p.Y[j]) + (Float64(p.Y[j+1]) - Float64(p.Y[j])) * (x - Float64(Xf[j])) / (Float64(Xf[j+1]) - Float64(Xf[j])))

# ===================================== lower bound: Sturm =====================================
"""
    sturm_cells(xs, hout) -> H::Vector{Float64}

H[k] ≥ max_{[xs[k], xs[k+1]]} hout, rigorously: hout is linear between its breakpoints and concave, so its max over a
cell is attained at a cell end or at a breakpoint inside the cell; all are evaluated in ball arithmetic and rounded up.
"""
function sturm_cells(xs::Vector{Float64}, houtA::PWLA)
    n = length(xs)
    vals = [up(eval_arb(houtA, x)) for x in xs]
    H = [max(vals[k], vals[k+1]) for k = 1:n-1]
    for (j, X) in enumerate(houtA.X)              # breakpoints strictly inside a cell
        k = searchsortedlast(xs, X)               # exact Float64-vs-Rational comparisons
        if 1 <= k <= n - 1 && xs[k] < X < xs[k+1]
            H[k] = max(H[k], up(houtA.YA[j]))
        end
    end
    return H
end

"""
    sturm_ok(xs, H, E) -> :ok | :fail | :straddle

Rigorous check (ball arithmetic) that the solution of −v'' + (V − E)v = 0, V = π²/H_k² on [xs[k], xs[k+1]],
v(0) = 0, v'(0) = 1, has no zero in (0, 1].  Cell rule: v > 0 at both cell ends (v'(0) > 0 on the first cell) and,
if V_k < E, ω_k·|J_k| < π with ω_k = √(E − V_k).  (Proof: BAND_README.md, Lemma C.)
"""
function sturm_ok(xs::Vector{Float64}, H::Vector{Float64}, E::Float64)
    EA = Arb(E)
    v = Arb(0); dv = Arb(1)
    for k in eachindex(H)
        H[k] > 0 || return :fail
        L = Arb(xs[k+1]) - Arb(xs[k])               # exact (difference of doubles at 128 bits)
        q = PI2_A / Arb(H[k])^2 - EA
        if Arblib.ispositive(q)
            κ = sqrt(q); x = κ * L
            c = cosh(x); s = sinh(x)
            vn = v * c + dv * (s / κ); dvn = v * (κ * s) + dv * c
        elseif Arblib.isnegative(q)
            ω = sqrt(-q); x = ω * L
            Arblib.ispositive(PI_A - x) || return :fail
            c = cos(x); s = sin(x)
            vn = v * c + dv * (s / ω); dvn = dv * c - v * (ω * s)
        else
            return :straddle
        end
        if k == 1
            Arblib.ispositive(dv) || return :fail
        else
            Arblib.ispositive(v) || return :fail
        end
        Arblib.ispositive(vn) || return :fail
        v = vn; dv = dvn
    end
    return :ok
end

# ===================================== upper bound: Ritz =====================================
"""
    Iints_arb(ha, hb) -> (I00, I01, I11)

Enclosures of ∫₀¹(1−s)²/h, ∫₀¹ s(1−s)/h, ∫₀¹ s²/h with h = ha(1−s) + hb·s, ha, hb > 0 (balls).
Let ρ̄ := (upper bound of |ε|)/2 with ε = (hb − ha)/hm, hm = (ha + hb)/2, an exact Arf.  If ρ̄ ≤ 1/4 (exact Arf
comparison): 1/h = (1/hm)·Σ_{n<K}(−εσ)^n + R_K, σ = s − 1/2, |εσ| ≤ ρ̄, |R_K| ≤ (1/hm)·ρ̄^K/(1−ρ̄); the partial sums use
exact rational coefficients and the tail bound ρ̄^K/(1−ρ̄)·(1/hm)·∫w (∫w = 1/3, 1/6, 1/3) is computed in Arb and added
as an Arf upper bound.  K is fixed per bracket of ρ̄ (16/24/32/48 for ρ̄ ≤ 1/64, 1/16, 1/8, 1/4; any K is valid, the
tail is always bounded rigorously).  Otherwise: the closed forms with log(hb/ha), all in Arb.
No Float64 arithmetic enters the result.
"""
function Iints_arb(ha::Arb, hb::Arb)
    hm = (ha + hb) * HALF_A
    Δ = hb - ha
    ε = Δ / hm
    ρA = Arb(Arblib.ubound(abs(ε))) * HALF_A          # exact: Arf upper bound of |ε|, times 1/2
    ρf = Arblib.ubound(ρA)                             # = ρA (radius 0); exact Arf
    if ρf <= Arf(0.25)
        K = ρf <= Arf(0.015625) ? 16 : (ρf <= Arf(0.0625) ? 24 : (ρf <= Arf(0.125) ? 32 : 48))
        se = Arb(0); so = Arb(0); sb = Arb(0); pw = Arb(1)
        for n = 0:K-1
            if iseven(n)
                se += pw * SERC_A[n+1]; sb += pw * SERB_A[n+1]
            else
                so += pw * SERC_A[n+1]
            end
            pw *= ε
        end
        ihm = inv(hm)
        I00 = (se + so) * ihm; I11 = (se - so) * ihm; I01 = sb * ihm
        ρa = Arb(ρf)
        tail = ρa^K / (1 - ρa) * ihm                    # ball; every point ≥ 0 bounds the truncation error · 1/hm
        t3 = Arblib.ubound(abs(tail / 3)); t6 = Arblib.ubound(abs(tail / 6))   # Arf upper bounds
        Arblib.add_error!(I00, t3); Arblib.add_error!(I11, t3); Arblib.add_error!(I01, t6)
        return I00, I01, I11
    else
        lg = log(hb / ha)
        Δ3 = Δ^3
        I11 = (Δ * (hb - 3 * ha) * HALF_A + ha^2 * lg) / Δ3
        I00 = (Δ * (ha - 3 * hb) * HALF_A + hb^2 * lg) / Δ3
        I01 = (Δ - ha * lg) / Δ^2 - I11
        return I00, I01, I11
    end
end

"""
    element_forms(xa, xb, aa, ab, ba, bb, isfirst, islast) -> (Kaa, Kab, Kbb, Maa, Mab, Mbb) or nothing

Local matrices (Arb balls) of N and D on the P1 element [xa, xb] with linear a, b (node values aa, ab, ba, bb).
Every Float64 input is converted to Arb BEFORE any arithmetic (exact).  isfirst / islast: end element on which h and f
vanish at the left / right node (h(0) = f(0) = 0, resp. at 1): only the other node's entries are nonzero, and
∫fg/h = L f g/(2h) exactly.  nothing if h is not provably positive where required.
"""
function element_forms(xa::Float64, xb::Float64, aa_::Float64, ab_::Float64, ba_::Float64, bb_::Float64, isfirst::Bool, islast::Bool)
    L = Arb(xb) - Arb(xa)
    aa = Arb(aa_); ab = Arb(ab_); ba = Arb(ba_); bb = Arb(bb_)
    ha = ba - aa; hb = bb - ab
    ap = (ab - aa) / L; bp = (bb - ba) / L; hp = bp - ap
    ce = HALF_A + ap * bp * HALF_A + hp^2 * C6_A
    k1 = (ha + hb) / (4 * L)
    pc = PI2_A * ce * L
    z = Arb(0)
    if isfirst
        Arblib.ispositive(hb) || return nothing
        return (z, z, k1 + hp / 4 + pc / (2 * hb), z, z, L * 3 * hb / 24)
    elseif islast
        Arblib.ispositive(ha) || return nothing
        return (k1 - hp / 4 + pc / (2 * ha), z, z, L * 3 * ha / 24, z, z)
    else
        (Arblib.ispositive(ha) && Arblib.ispositive(hb)) || return nothing
        I00, I01, I11 = Iints_arb(ha, hb)
        return (k1 - hp / 4 + pc * I00, -k1 + pc * I01, k1 + hp / 4 + pc * I11,
                L * (3 * ha + hb) / 24, L * (ha + hb) / 24, L * (ha + 3 * hb) / 24)
    end
end

"""
    gram_arb(xs, a, b, F) -> (N::Matrix{Arb}, D::Matrix{Arb}) or nothing

Rigorous Gram matrices N_ij = N(f_i, f_j), D_ij = D(f_i, f_j) of the weighted form for the P1 functions with nodes xs
and values F[:, i] (F[1, :] = F[end, :] = 0) on the domain {a < y < b}, a, b the P1 interpolants of the node values
(exact floats, a[1] = b[1] = a[end] = b[end] = 0, b − a > 0 at interior nodes).  The nodal values are promoted to Arb
before any product (fixes the Float64-product defect found by the 2026-10-01 review).
"""
function gram_arb(xs::Vector{Float64}, a::Vector{Float64}, b::Vector{Float64}, F::Matrix{Float64})
    n = length(xs); m = size(F, 2)
    @assert n >= 3 && length(a) == n && length(b) == n && size(F, 1) == n
    @assert a[1] == 0 && b[1] == 0 && a[end] == 0 && b[end] == 0
    @assert all(F[1, :] .== 0) && all(F[end, :] .== 0)
    FA = [Arb(F[e, i]) for e = 1:n, i = 1:m]            # exact conversion, before any arithmetic
    N = [Arb(0) for _ = 1:m, _ = 1:m]; D = [Arb(0) for _ = 1:m, _ = 1:m]
    for e = 1:n-1
        loc = element_forms(xs[e], xs[e+1], a[e], a[e+1], b[e], b[e+1], e == 1, e == n - 1)
        loc === nothing && return nothing
        Kaa, Kab, Kbb, Maa, Mab, Mbb = loc
        for i = 1:m, j = i:m
            fa = FA[e, i]; fb = FA[e+1, i]; ga = FA[e, j]; gb = FA[e+1, j]
            paa = fa * ga; pab = fa * gb + fb * ga; pbb = fb * gb          # Arb (exact at 128 bits)
            N[i, j] += Kaa * paa + Kab * pab + Kbb * pbb
            D[i, j] += Maa * paa + Mab * pab + Mbb * pbb
        end
    end
    for i = 1:m, j = 1:i-1
        N[i, j] = N[j, i]; D[i, j] = D[j, i]
    end
    return N, D
end

# ===================================== the box evaluation =====================================
struct BoxData
    hout::PWL; ain::PWL; bin::PWL
    xs::Vector{Float64}       # Sturm nodes
    xr::Vector{Float64}       # Ritz nodes
end

function meshes(B::Box, prm::Params)
    hout, ain, bin = profiles(B)
    xs = graded_nodes(Float64.(hout.X), Float64.(hout.Y); delta = prm.dS, p = prm.pS, dxmax = prm.dxS)
    Xi = sort(unique(vcat(Float64.(ain.X), Float64.(bin.X))))
    Yi = [evalf(bin, x) - evalf(ain, x) for x in Xi]
    xr = graded_nodes(Xi, Yi; delta = prm.dR, p = prm.pR, dxmax = prm.dxR)
    return BoxData(hout, ain, bin, xs, xr)
end

"""
    inner_nodes(bd) -> (a, b)

Float node values of the Ritz profile, rounded INWARD: a[j] ≥ a_in(x_j), b[j] ≤ b_in(x_j) (ball upper/lower bounds),
exactly 0 at x = 0, 1.  Since a_in is convex and b_in concave, the P1 interpolants satisfy a_int ≥ a_in, b_int ≤ b_in
on [0,1], so {a_int < y < b_int} ⊂ Q_in.
"""
function inner_nodes(bd::BoxData)
    ainA = PWLA(bd.ain); binA = PWLA(bd.bin)
    a = [up(eval_arb(ainA, x)) for x in bd.xr]
    b = [dn(eval_arb(binA, x)) for x in bd.xr]
    @assert bd.ain.Y[1] == 0 && bd.ain.Y[end] == 0 && bd.bin.Y[1] == 0 && bd.bin.Y[end] == 0
    a[1] = a[end] = b[1] = b[end] = 0.0
    return a, b
end

"float estimate of (μ1 lower, R2 upper) for the box (non-rigorous; used for E, f and the bisection heuristic)"
function float_eval(B::Box, prm::Params = DEFAULT_PARAMS)
    bd = meshes(B, prm)
    xs = bd.xs
    vals = [evalf(bd.hout, x) for x in xs]
    H = [max(vals[k], vals[k+1]) for k = 1:length(xs)-1]
    for (j, X) in enumerate(bd.hout.X)
        k = searchsortedlast(xs, X)
        if 1 <= k <= length(xs) - 1 && xs[k] < X < xs[k+1]; H[k] = max(H[k], Float64(bd.hout.Y[j])); end
    end
    V = pi^2 ./ H .^ 2
    μ1 = mu1_step_f(xs, V)
    a = [evalf(bd.ain, x) for x in bd.xr]; b = [evalf(bd.bin, x) for x in bd.xr]
    a[1] = a[end] = b[1] = b[end] = 0.0
    R, F = ritz2_f(bd.xr, a, b)
    return μ1, R[2], bd, F
end

"""
    rigorous_bounds(B, prm; μ1f, bd, F) -> (E, R2A, D_pd, why)

E: a float with λ1(Q) > E for all Q in B (Sturm-verified; NaN if the check failed).  R2A: ball containing the larger
Ritz value R2 ≥ λ2(Q) for all convex Q in B (NaN ball if D is not positive definite).  No 5/2 gate here.
"""
function rigorous_bounds(B::Box, prm::Params, bd::BoxData, F::Matrix{Float64})
    houtA = PWLA(bd.hout)
    H = sturm_cells(bd.xs, houtA)
    V = pi^2 ./ H .^ 2
    μ1s = mu1_step_f(bd.xs, V)
    E = μ1s * (1 - prm.Erel)
    st = sturm_ok(bd.xs, H, E)
    tries = 0
    while st != :ok && tries < 3
        tries += 1
        E *= (1 - 1e-5 * 10.0^tries)
        st = sturm_ok(bd.xs, H, E)
    end
    st == :ok || return (E = NaN, R2A = Arb(NaN), N = nothing, D = nothing, why = "sturm $st")
    a, b = inner_nodes(bd)
    G = gram_arb(bd.xr, a, b, F)
    G === nothing && return (E = E, R2A = Arb(NaN), N = nothing, D = nothing, why = "inner profile not positive")
    N, D = G
    detD = D[1, 1] * D[2, 2] - D[1, 2]^2
    (Arblib.ispositive(D[1, 1]) && Arblib.ispositive(detD)) || return (E = E, R2A = Arb(NaN), N = N, D = D, why = "D not PD")
    bq = N[1, 1] * D[2, 2] + N[2, 2] * D[1, 1] - 2 * N[1, 2] * D[1, 2]
    detN = N[1, 1] * N[2, 2] - N[1, 2]^2
    disc = bq^2 - 4 * detD * detN
    R2 = Arblib.ispositive(disc) ? (bq + sqrt(disc)) / (2 * detD) : Arb(NaN)
    return (E = E, R2A = R2, N = N, D = D, why = "")
end

"""
    certify(B; prm) -> NamedTuple (ok, E, R2lo, R2hi, nS, nR, why, R2A)

Full certificate.  ok = true means (rigorously): λ1(Q) > E (Sturm, Lemma C) and (5/2)E·D − N ≻ 0 for the rigorous
2×2 Gram matrices with D ≻ 0, i.e. λ2(Q) ≤ R2 < (5/2)E, for every convex Q in B.
"""
function certify(B::Box, prm::Params = DEFAULT_PARAMS)
    μ1f, R2f, bd, Ffl = float_eval(B, prm)
    nS = length(bd.xs) - 1; nR = length(bd.xr) - 1
    if !(μ1f > 0) || !(R2f < 2.5 * μ1f * (1 - 1e-6))
        return (ok = false, E = μ1f, R2lo = R2f, R2hi = R2f, nS = nS, nR = nR, why = "float ratio $(R2f/μ1f)", R2A = Arb(NaN))
    end
    rb = rigorous_bounds(B, prm, bd, Ffl)
    isfinite(rb.E) && rb.D !== nothing || return (ok = false, E = rb.E, R2lo = R2f, R2hi = R2f, nS = nS, nR = nR, why = rb.why, R2A = Arb(NaN))
    N, D = rb.N, rb.D
    rb.why == "" || return (ok = false, E = rb.E, R2lo = R2f, R2hi = R2f, nS = nS, nR = nR, why = rb.why, R2A = Arb(NaN))
    Λ = Arb(rb.E) * 5 / 2
    A11 = Λ * D[1, 1] - N[1, 1]; A12 = Λ * D[1, 2] - N[1, 2]; A22 = Λ * D[2, 2] - N[2, 2]
    okR = Arblib.ispositive(A11) && Arblib.ispositive(A11 * A22 - A12^2)
    return (ok = okR, E = rb.E, R2lo = dn(rb.R2A), R2hi = up(rb.R2A), nS = nS, nR = nR, why = okR ? "" : "R2 >= 2.5E", R2A = rb.R2A)
end
