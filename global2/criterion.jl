# criterion.jl — the box criterion on the base TMs, rigorous scalar eigenvalue bounds over a box,
# and the exclusion tests E1–E4.
#
# Criterion (κ = 5/2 for the theorem; smaller κ gives a certified bound ξ < κ):
#   LG mode:    Ã = K·(βA2 − A1) − κ(βA1 − A0)·M,  A0 = K11, A1 = M11, A2 = ‖w‖²,  β ≤ λ2 (constant on the box)
#   MONO mode:  Ã = K − κ L0·M,  L0 ≤ λ1 (constant on the box)
#   Ã11 < 0 and S' = −Ã11 Ã_cc + Ã_c1 Ã_1c (c = {2,3}) has a negative eigenvalue (tr S' < 0, det S' < 0,
#   S'11 < 0, S'22 < 0, or zᵀS'z < 0 for a fixed z) ⇒ Ã has ≥ 2 negative eigenvalues (Haynsworth) ⇒ λ2 < κ·λ1.
#   Side conditions: LG: βA1 − A0 > 0 (then βA2 − A1 > 0 by Cauchy–Schwarz; also checked); M ≻ 0 (checked).

# ------------------------------------------------------------------ exact region constants
const RAT_ETA = 35 // 1000
const RAT_S0 = 1 // 10
"exact rational value of a decimal literal string such as \"0.039\""
function dec_to_rat(s::AbstractString)
    s = strip(s)
    neg = startswith(s, "-"); s = lstrip(s, ['-', '+'])
    if occursin('.', s)
        ip, fp = split(s, '.')
        r = parse(BigInt, ip * fp) // big(10)^length(fp)
    else
        r = parse(BigInt, s) // 1
    end
    return neg ? -r : r
end
"read a numeric key from config.toml as the exact decimal written there"
function config_rat(path::AbstractString, key::AbstractString)
    for ln in eachline(path)
        m = match(Regex("^\\s*" * key * "\\s*=\\s*([-+0-9.]+)"), ln)
        m === nothing || return dec_to_rat(m.captures[1])
    end
    error("key $key not found in $path")
end

# ------------------------------------------------------------------ geometry TMs for a box (Arb or Float64)
struct BoxGeom{T,R}
    V::Matrix{TM{T,R}}     # 4×2 vertex TMs
    E::Matrix{TM{T,R}}     # edges
    A::Matrix{TM{T,R}}     # A[k,j] = cross(e_k, v_j − v_k)  (c_k = A[k, k+2])
end
function box_geom(::Type{T}, p0::Vector{Float64}, rad::Vector{Float64}, ws) where {T}
    TMT = tmtype(T)
    V = [newtm(T) for _ = 1:4, _ = 1:2]; E = [newtm(T) for _ = 1:4, _ = 1:2]; A = [newtm(T) for _ = 1:4, _ = 1:4]
    for j = 1:4, k = 1:2
        if T == Arb
            c0 = Arb(V0[j, k]; prec = PREC)
            for i = 1:4; Arblib.addmul!(c0, Arb(DV[j, k, i]; prec = PREC), Arb(p0[i]; prec = PREC)); end
            tm_linear!(V[j, k], c0, [Arb(DV[j, k, i]; prec = PREC) * Arb(rad[i]; prec = PREC) for i = 1:4])
        else
            tm_linear!(V[j, k], V0[j, k] + sum(DV[j, k, i] * p0[i] for i = 1:4), [DV[j, k, i] * rad[i] for i = 1:4])
        end
    end
    for k = 1:4, c = 1:2; tm_sub!(E[k, c], V[mod1(k + 1, 4), c], V[k, c]); end
    t1 = newtm(T); t2 = newtm(T); t3 = newtm(T)
    for k = 1:4, j = 1:4
        (j == k || j == mod1(k + 1, 4)) && continue
        tm_sub!(t1, V[j, 1], V[k, 1]); tm_sub!(t2, V[j, 2], V[k, 2])
        tm_mul!(A[k, j], E[k, 1], t2, ws); tm_mul!(t3, E[k, 2], t1, ws); tm_sub!(A[k, j], A[k, j], t3)
    end
    BoxGeom{T,tmtype(T).parameters[2]}(V, E, A)
end
tsc(::Type{Float64}, x::Rational) = Float64(x)
tsc(::Type{Arb}, x::Rational) = Arb(numerator(x); prec = PREC) / Arb(denominator(x); prec = PREC)

# |u|² for a vector TM pair, into h
function tm_sq2!(h, ux, uy, ws)
    tm_mul!(h, ux, ux, ws); tm_addmul!(h, uy, uy, ws); h
end

# ------------------------------------------------------------------ exclusion tests (rigorous in Arb mode)
"""
    exclusion(lo, hi; cube) -> Symbol

:E1 (box inside the closed cube |p|∞ ≤ cube, exact rational test), :E4 (some c_i ≤ 0 on the whole box),
:E2 (some c_i ≤ η c_{i+2} on the whole box), :E3 (width/diam ≤ s0 on the whole box), or :none.
All tests are evaluated with Arb Taylor models over a box that contains [lo, hi].
"""
function exclusion(lo::Vector{Float64}, hi::Vector{Float64}; cube::Rational{BigInt}, ws = newwork(Arb))
    # E1: lo, hi are exact binary rationals
    if all(abs(Rational{BigInt}(lo[i])) <= cube && abs(Rational{BigInt}(hi[i])) <= cube for i = 1:4)
        return :E1
    end
    p0, rad = box_center_radius(lo, hi)
    G = box_geom(Arb, p0, rad, ws)
    # E4
    for i = 1:4
        tm_upper(G.A[i, mod1(i + 2, 4)]) <= 0 && return :E4
    end
    # E2: c_i − η c_{i+2} ≤ 0
    eta = tsc(Arb, RAT_ETA); h = newtm(Arb)
    for i = 1:4
        tm_set!(h, G.A[i, mod1(i + 2, 4)]); tm_axpy!(h, -eta, G.A[mod1(i + 2, 4), i], ws)
        tm_upper(h) <= 0 && return :E2
    end
    # E3: for some edge i and some vertex pair (k,l): cross(e_i, v_j − v_i)² ≤ s0² |e_i|² |v_k − v_l|² for both j ≠ i, i+1
    s02 = tsc(Arb, RAT_S0^2)
    Vf = verts(p0)
    pairs = [(k, l) for k = 1:4 for l = k+1:4]
    sort!(pairs, by = kl -> -hypot(Vf[kl[1], 1] - Vf[kl[2], 1], Vf[kl[1], 2] - Vf[kl[2], 2]))
    e2 = newtm(Arb); d2 = newtm(Arb); dx = newtm(Arb); dy = newtm(Arb); rhs = newtm(Arb); lhs = newtm(Arb)
    for (k, l) in pairs[1:2]
        tm_sub!(dx, G.V[k, 1], G.V[l, 1]); tm_sub!(dy, G.V[k, 2], G.V[l, 2]); tm_sq2!(d2, dx, dy, ws)
        for i = 1:4
            tm_sq2!(e2, G.E[i, 1], G.E[i, 2], ws)
            tm_mul!(rhs, e2, d2, ws); tm_scale!(rhs, rhs, s02, ws)
            okall = true
            for j in (mod1(i + 2, 4), mod1(i + 3, 4))
                tm_mul!(lhs, G.A[i, j], G.A[i, j], ws); tm_sub!(lhs, lhs, rhs)
                if !(tm_upper(lhs) <= 0); okall = false; break; end
            end
            okall && return :E3
        end
    end
    return :none
end

"centre and radius (Float64) of a TM box that provably contains [lo, hi] (checked in exact arithmetic)."
function box_center_radius(lo::Vector{Float64}, hi::Vector{Float64})
    p0 = (lo .+ hi) ./ 2
    rad = similar(p0)
    for i = 1:4
        r = max(p0[i] - lo[i], hi[i] - p0[i])
        r = nextfloat(r)
        @assert Rational{BigInt}(p0[i]) - Rational{BigInt}(r) <= Rational{BigInt}(lo[i]) &&
                Rational{BigInt}(p0[i]) + Rational{BigInt}(r) >= Rational{BigInt}(hi[i])
        rad[i] = r
    end
    return p0, rad
end

# ------------------------------------------------------------------ rigorous scalar bounds over a box (Arb)
const PI_A = Arb(π; prec = PREC)
const J01_A = Arb((Arb("2.4048255576957727"; prec = PREC), Arb("2.4048255576957729"; prec = PREC)); prec = PREC)

# sup over the box of an affine function c0 + Σ c_i δ_i  (exact: c0 + Σ |c_i| rad_i)
function affine_sup(c0::Arb, c::Vector{Arb}, rad::Vector{Float64})
    s = Arb(c0; prec = PREC)
    for i = 1:4; s += abs(c[i]) * Arb(rad[i]; prec = PREC); end
    s
end

"""
    scalar_bounds(p0, rad) -> NamedTuple (Float64 numbers, rigorous)

beta_rect_ks : lower bound of λ2(Q_p) on the box from Krahn–Szegő and enclosing rectangles,
L_rect       : lower bound of λ1 from enclosing rectangles (domain monotonicity),
L_hersch     : lower bound of λ1 from Hersch λ1 ≥ π²/(4 r²) (convex Q; r = inradius, upper-bounded),
gapAC        : 3π²/D_max² (Andrews–Clutterbuck: λ2 − λ1 ≥ 3π²/D² for convex Q).
Domain of validity (for p with Q_p convex): the rectangle, Hersch and Andrews–Clutterbuck parts are computed from
p0 and rad in ball arithmetic and hold on the closed box p0 ± rad; the Krahn–Szegő area bound uses the Float64
endpoints fl(p0[i] ∓ rad[i]) of the a- and d-ranges, so it holds on [fl(p0 − rad), fl(p0 + rad)] in those two
coordinates.  That box contains every logged leaf [lo, hi]: box_center_radius guarantees p0 − rad ≤ lo and
hi ≤ p0 + rad exactly, and rounding to nearest is monotone with lo, hi themselves Float64, so fl(p0 − rad) ≤ lo and
hi ≤ fl(p0 + rad).  Hence every number is a rigorous bound for all p in the leaf, which is the domain on which the
certificate uses it.
"""
function scalar_bounds(p0::Vector{Float64}, rad::Vector{Float64}; ws = newwork(Arb))
    a_lo = p0[1] - rad[1]; a_hi = p0[1] + rad[1]; d_lo = p0[4] - rad[4]; d_hi = p0[4] + rad[4]
    # |Q| = 1 − a² − d² ≤ 1 − min a² − min d²   (exact: ball arithmetic on the endpoints)
    amin2 = (a_lo <= 0 <= a_hi) ? Arb(0) : min(Arb(a_lo)^2, Arb(a_hi)^2)
    dmin2 = (d_lo <= 0 <= d_hi) ? Arb(0) : min(Arb(d_lo)^2, Arb(d_hi)^2)
    areamax = 1 - amin2 - dmin2
    ks = 2 * PI_A * J01_A^2 / areamax
    beta = Arblib.ispositive(areamax) ? ks : Arb(0)
    Vf = verts(p0)
    # vertex affine data in Arb: v_jk(p0+δ) = c0 + Σ DV δ
    V0A = [let c = Arb(V0[j, k]; prec = PREC); for i = 1:4; c += Arb(DV[j, k, i]) * Arb(p0[i]); end; c end for j = 1:4, k = 1:2]
    pi2 = PI_A^2
    Lrect = Arb(0)
    for th in candidate_dirs(Vf)
        cs = cos(th); sn = sin(th)              # any Float64 direction; normalised rigorously below
        nrm = sqrt(Arb(cs)^2 + Arb(sn)^2)
        Amax = Arb(-Inf); Bmax = Arb(-Inf)
        for j = 1:4, k = 1:4
            j == k && continue
            # (v_j − v_k)·(cs, sn) and (v_j − v_k)·(−sn, cs): affine in δ
            cA = (V0A[j, 1] - V0A[k, 1]) * cs + (V0A[j, 2] - V0A[k, 2]) * sn
            cB = -(V0A[j, 1] - V0A[k, 1]) * sn + (V0A[j, 2] - V0A[k, 2]) * cs
            gA = [Arb(DV[j, 1, i] - DV[k, 1, i]) * cs + Arb(DV[j, 2, i] - DV[k, 2, i]) * sn for i = 1:4]   # DV differences exact
            gB = [-(Arb(DV[j, 1, i] - DV[k, 1, i]) * sn) + Arb(DV[j, 2, i] - DV[k, 2, i]) * cs for i = 1:4]
            Amax = max(Amax, affine_sup(cA, gA, rad)); Bmax = max(Bmax, affine_sup(cB, gB, rad))
        end
        A = Arb(Arblib.ubound(Amax)) / nrm; B = Arb(Arblib.ubound(Bmax)) / nrm
        lam2 = pi2 * min(4 / A^2 + 1 / B^2, 1 / A^2 + 4 / B^2)
        lam1 = pi2 * (1 / A^2 + 1 / B^2)
        beta = max(beta, Arb(Arblib.lbound(lam2)))
        Lrect = max(Lrect, Arb(Arblib.lbound(lam1)))
    end
    # diameter upper bound: max over pairs of sup |v_j − v_k|
    G = box_geom(Arb, p0, rad, ws)
    d2 = newtm(Arb); dx = newtm(Arb); dy = newtm(Arb)
    D2max = 0.0
    for j = 1:4, k = j+1:4
        tm_sub!(dx, G.V[j, 1], G.V[k, 1]); tm_sub!(dy, G.V[j, 2], G.V[k, 2]); tm_sq2!(d2, dx, dy, ws)
        D2max = max(D2max, tm_upper(d2))
    end
    gapAC = 3 * pi2 / Arb(D2max)
    # inradius upper bound: strip bounds w_i/2 and certified-bounded drop-triangles
    rmax = Arb(Inf)
    e2 = newtm(Arb)
    for i = 1:4
        tm_sq2!(e2, G.E[i, 1], G.E[i, 2], ws)
        elo = tm_lower(e2)
        elo > 0 || continue
        wmax = max(tm_upper(G.A[i, mod1(i + 2, 4)]), tm_upper(G.A[i, mod1(i + 3, 4)]))
        rmax = min(rmax, Arb(wmax) / (2 * sqrt(Arb(elo))))
    end
    # drop-triangles via interval vertices
    VA = [Arb((Arblib.lbound(V0A[j, k] - sum(abs(Arb(DV[j, k, i])) * Arb(rad[i]) for i = 1:4)),
               Arblib.ubound(V0A[j, k] + sum(abs(Arb(DV[j, k, i])) * Arb(rad[i]) for i = 1:4)))) for j = 1:4, k = 1:2]
    cr = newtm(Arb)
    for l = 1:4
        lm = mod1(l - 1, 4); lp = mod1(l + 1, 4); lo2 = mod1(l + 2, 4)
        # bounded iff cross(e_{l−1}, e_{l+1}) > 0 on the whole box
        tm_mul!(cr, G.E[lm, 1], G.E[lp, 2], ws); tm_mul!(d2, G.E[lm, 2], G.E[lp, 1], ws); tm_sub!(cr, cr, d2)
        tm_lower(cr) > 0 || continue
        r = tri_inradius_interval(VA, (lm, lp, lo2))
        Arblib.isfinite(r) && (rmax = min(rmax, Arb(Arblib.ubound(r))))
    end
    Lh = Arblib.isfinite(rmax) && Arblib.ispositive(rmax) ? pi2 / (4 * rmax^2) : Arb(0)
    lb(x) = Float64(Arblib.lbound(x), RoundDown)
    return (; beta_rect_ks = lb(beta), L_rect = lb(Lrect), L_hersch = lb(Lh), gapAC = lb(gapAC), D2max)
end

# inradius of the triangle cut out by the lines of edges (k1,k2,k3) of a quadrilateral with interval vertices:
# triangle vertices are pairwise line intersections; r = 2|T|/perimeter.
function tri_inradius_interval(VA, ks)
    P = Vector{Tuple{Arb,Arb}}()
    lines = [(VA[k, 1], VA[k, 2], VA[mod1(k + 1, 4), 1] - VA[k, 1], VA[mod1(k + 1, 4), 2] - VA[k, 2]) for k in ks]
    for (a, b) in ((1, 2), (2, 3), (3, 1))
        (px, py, ex, ey) = lines[a]; (qx, qy, fx, fy) = lines[b]
        den = ex * fy - ey * fx
        Arblib.contains_zero(den) && return Arb(Inf)
        s = ((qx - px) * fy - (qy - py) * fx) / den
        push!(P, (px + s * ex, py + s * ey))
    end
    area2 = abs((P[2][1] - P[1][1]) * (P[3][2] - P[1][2]) - (P[2][2] - P[1][2]) * (P[3][1] - P[1][1]))
    per = sqrt((P[1][1] - P[2][1])^2 + (P[1][2] - P[2][2])^2) + sqrt((P[2][1] - P[3][1])^2 + (P[2][2] - P[3][2])^2) +
          sqrt((P[3][1] - P[1][1])^2 + (P[3][2] - P[1][2])^2)
    return area2 / per
end

# ------------------------------------------------------------------ the criterion on base TMs
struct CritPool{T,R}
    At::Vector{TM{T,R}}; g1::TM{T,R}; g2::TM{T,R}
    s11::TM{T,R}; s12::TM{T,R}; s22::TM{T,R}; tr::TM{T,R}; det::TM{T,R}; zsz::TM{T,R}
    t1::TM{T,R}; t2::TM{T,R}; t3::TM{T,R}; m2::TM{T,R}; m3::TM{T,R}
    q::TM{T,R}; Mn::Vector{TM{T,R}}; Kn::Vector{TM{T,R}}; A2n::TM{T,R}
    buf::Vector{T}
end
function critpool(::Type{T}) where {T}
    TMT = tmtype(T)
    CritPool{T,TMT.parameters[2]}([newtm(T) for _ = 1:6], (newtm(T) for _ = 1:14)...,
        [newtm(T) for _ = 1:6], [newtm(T) for _ = 1:6], newtm(T), [tconv(T, 0.0) for _ = 1:FTAB[].nc2])
end
_fmul!(cp::CritPool, h, f, g, ws) = tm_mul_full!(h, f, g, ws, cp.buf)
function _faddmul!(cp::CritPool, h, f, g, ws)       # h += f g  (exact product, folded)
    tm_mul_full!(cp.t3, f, g, ws, cp.buf); tm_add!(h, h, cp.t3)
end

"""
    choose_q(c) -> Vector{Float64}

Coefficients of a polynomial q(t) ≈ 1/M11(t) (truncated geometric series; any choice is admissible — it is
only used as a common positive factor, positivity is certified separately).  `c` = coefficients of M11.
"""
function choose_q(c::Vector{Float64})
    K = TAB[].order; F = c[1]; ws = newwork(Float64)
    u = newtm(Float64); copyto!(u.c, c ./ F); u.c[1] = 0.0
    q = newtm(Float64); tm_const!(q, 1.0); up = newtm(Float64); tm_const!(up, 1.0); tmp = newtm(Float64)
    for j = 1:K
        tm_mul!(tmp, up, u, ws); tm_set!(up, tmp); tm_axpy!(q, (-1.0)^j, up)
    end
    return q.c ./ F
end

"""
    normalize!(cp, M, K, A2, ws) -> Bool

Multiply all base TMs by the common polynomial q(t) ≈ 1/M11(t) (thin, exactly representable coefficients).
Every criterion polynomial is homogeneous in (K, M, A2), so its sign is unchanged when q > 0 on the box;
this is certified here (lower bound of q > 0).  Products are exact (full product, folded).
"""
function normalize!(cp::CritPool{T}, M, K, A2, ws) where {T}
    cM = T == Arb ? [Float64(Arblib.midref(x)) for x in M[1].c] : copy(M[1].c)
    qc = choose_q(cM)
    tm_zero!(cp.q)
    for k in eachindex(qc); _setc!(cp.q, k, tconv(T, qc[k])); end      # Float64 → exact thin ball
    tm_lower(cp.q) > 0 || return false
    for k = 1:6
        _fmul!(cp, cp.Mn[k], cp.q, M[k], ws); _fmul!(cp, cp.Kn[k], cp.q, K[k], ws)
    end
    _fmul!(cp, cp.A2n, cp.q, A2, ws)
    return true
end

"""
    crit_eval!(cp, M, K, A2, ws; mode, beta, L0, kappa, z, normalize) -> NamedTuple

Evaluates the criterion TMs and their bounds.  `ok` = every side condition and at least one of the
negativity tests certified.  Rigorous in Arb mode.
"""
function crit_eval!(cp::CritPool{T}, M, K, A2, ws; mode::Symbol, beta::Union{Float64,TM} = 0.0, L0::Union{Float64,TM} = 0.0,
                    kappa::Rational = 5 // 2, z = (1.0, 0.0), normalize::Bool = true) where {T}
    if normalize
        normalize!(cp, M, K, A2, ws) || return (; ok = false, which = :qneg, a11up = NaN, ups = (tr = NaN, det = NaN, s11 = NaN, s22 = NaN, zsz = NaN), g1lo = NaN, g2lo = NaN, mpos = false, g1ok = false)
        M = cp.Mn; K = cp.Kn; A2 = cp.A2n
    end
    kap = tsc(T, kappa)
    At = cp.At
    g1ok = true; g2lo = NaN; g1lo = NaN
    if mode === :LG
        # g1 = β A2 − A1,  g2 = β A1 − A0   (β a constant or a pointwise-valid Taylor model β(p) ≤ λ2(Q_p))
        if beta isa TM
            _fmul!(cp, cp.g1, beta, A2, ws); tm_axpy!(cp.g1, tconv(T, -1.0), M[1], ws)
            _fmul!(cp, cp.g2, beta, M[1], ws); tm_axpy!(cp.g2, tconv(T, -1.0), K[1], ws)
        else
            b = tconv(T, beta)
            tm_scale!(cp.g1, A2, b, ws); tm_axpy!(cp.g1, tconv(T, -1.0), M[1], ws)
            tm_scale!(cp.g2, M[1], b, ws); tm_axpy!(cp.g2, tconv(T, -1.0), K[1], ws)
        end
        g2lo = tm_lower(cp.g2); g1lo = tm_lower(cp.g1)
        g1ok = g2lo > 0 && g1lo > 0
        for q = 1:6
            _fmul!(cp, At[q], K[q], cp.g1, ws)
            _fmul!(cp, cp.t1, cp.g2, M[q], ws); tm_axpy!(At[q], -kap, cp.t1, ws)
        end
    else
        if L0 isa TM          # Ã = K − κ L0(p) M,  L0(p) ≤ λ1(Q_p) pointwise
            for q = 1:6
                _fmul!(cp, cp.t1, L0, M[q], ws); tm_set!(At[q], K[q]); tm_axpy!(At[q], -kap, cp.t1, ws)
            end
        else
            cL = T == Arb ? -(kap * Arb(L0; prec = PREC)) : -(Float64(kappa) * L0)
            for q = 1:6
                tm_set!(At[q], K[q]); tm_axpy!(At[q], cL, M[q], ws)
            end
        end
    end
    # S' entries: s11 = −Ã11 Ã22 + Ã12²,  s12 = −Ã11 Ã23 + Ã12 Ã13,  s22 = −Ã11 Ã33 + Ã13²
    A11, A12, A13, A22, A23, A33 = At
    m1 = tconv(T, -1.0)
    _fmul!(cp, cp.s11, A11, A22, ws); tm_scale!(cp.s11, cp.s11, m1, ws); _faddmul!(cp, cp.s11, A12, A12, ws)
    _fmul!(cp, cp.s12, A11, A23, ws); tm_scale!(cp.s12, cp.s12, m1, ws); _faddmul!(cp, cp.s12, A12, A13, ws)
    _fmul!(cp, cp.s22, A11, A33, ws); tm_scale!(cp.s22, cp.s22, m1, ws); _faddmul!(cp, cp.s22, A13, A13, ws)
    tm_add!(cp.tr, cp.s11, cp.s22)
    _fmul!(cp, cp.det, cp.s11, cp.s22, ws); _fmul!(cp, cp.t1, cp.s12, cp.s12, ws); tm_sub!(cp.det, cp.det, cp.t1)
    z1 = tconv(T, z[1]); z2 = tconv(T, z[2])
    tm_scale!(cp.zsz, cp.s11, z1 * z1, ws); tm_axpy!(cp.zsz, 2 * z1 * z2, cp.s12, ws); tm_axpy!(cp.zsz, z2 * z2, cp.s22, ws)
    a11up = tm_upper(A11)
    ups = (tr = tm_upper(cp.tr), det = tm_upper(cp.det), s11 = tm_upper(cp.s11), s22 = tm_upper(cp.s22), zsz = tm_upper(cp.zsz))
    # M ≻ 0: M11 > 0, M11 M22 − M12² > 0, det M > 0
    M11, M12, M13, M22, M23, M33 = M
    _fmul!(cp, cp.m2, M11, M22, ws); _fmul!(cp, cp.t1, M12, M12, ws); tm_sub!(cp.m2, cp.m2, cp.t1)
    # det M = M11 (M22 M33 − M23²) − M12 (M12 M33 − M23 M13) + M13 (M12 M23 − M22 M13)
    _fmul!(cp, cp.t1, M22, M33, ws); _fmul!(cp, cp.t2, M23, M23, ws); tm_sub!(cp.t1, cp.t1, cp.t2); _fmul!(cp, cp.m3, M11, cp.t1, ws)
    _fmul!(cp, cp.t1, M12, M33, ws); _fmul!(cp, cp.t2, M23, M13, ws); tm_sub!(cp.t1, cp.t1, cp.t2); _fmul!(cp, cp.t2, M12, cp.t1, ws); tm_sub!(cp.m3, cp.m3, cp.t2)
    _fmul!(cp, cp.t1, M12, M23, ws); _fmul!(cp, cp.t2, M22, M13, ws); tm_sub!(cp.t1, cp.t1, cp.t2); _fmul!(cp, cp.t2, M13, cp.t1, ws); tm_add!(cp.m3, cp.m3, cp.t2)
    mpos = tm_lower(M11) > 0 && tm_lower(cp.m2) > 0 && tm_lower(cp.m3) > 0
    which = :none
    for s in (:tr, :det, :zsz, :s11, :s22)
        if getfield(ups, s) < 0; which = s; break; end
    end
    ok = g1ok && mpos && a11up < 0 && which !== :none
    return (; ok, which, a11up, ups, g1lo, g2lo, mpos, g1ok)
end
