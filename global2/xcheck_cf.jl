# xcheck_cf.jl — independent exclusion checks by centred forms of the exact polynomials (no Taylor-model code).
#
# Vertices are affine in p: v_j(m + δ) = v_j(m) + D_j δ with D_j = DV[j,:,:] (exact dyadic numbers).  For affine
# vectors a, b the scalars cross(a,b), dot(a,b) are quadratics  q(δ) = c + g·δ + δᵀQδ  with c, g balls (Arb, centre
# m) and Q exact (Float64 products of dyadics).  On a box |δ_k| ≤ r_k:
#   q ≤ c + Σ|g_k| r_k + Σ_k max(Q_kk,0) r_k² + Σ_{k<l} |Q_kl+Q_lk| r_k r_l          (E2, E4: quadratic conditions)
# For the thin conditions F = cross(e_i, v_j − v_i)² − s²|e_i|²|v_k − v_l|² (quartic) the Taylor form with Lagrange
# remainder is used: F(δ) ≤ F(0) + Σ_k |∂_kF(0)| r_k + ½ Σ_kl max|∂²_klF| r_k r_l, the Hessian being enclosed over the
# box from enclosures of the quadratic factors and their gradients.  A condition is proved on a box by adaptive
# bisection (widest side, Float64 midpoint; the two halves cover the box) for ONE fixed condition index.
cfarb(x) = Arb(x; prec = PREC)
struct CFQ
    c::Arb
    g::Vector{Arb}
    Q::Matrix{Float64}
end
function cf_vert(j, m)
    vm = [let c = cfarb(V0[j, k]); for i = 1:4; c += cfarb(DV[j, k, i]) * cfarb(m[i]); end; c end for k = 1:2]
    return (vm, [DV[j, k, i] for k = 1:2, i = 1:4])
end
cf_sub(a, b) = (a[1] .- b[1], a[2] .- b[2])
function cf_cross(a, b)
    am, A = a; bm, B = b
    c = am[1] * bm[2] - am[2] * bm[1]
    g = [am[1] * B[2, k] - am[2] * B[1, k] + A[1, k] * bm[2] - A[2, k] * bm[1] for k = 1:4]
    Q = [A[1, k] * B[2, l] - A[2, k] * B[1, l] for k = 1:4, l = 1:4]
    CFQ(c, g, Q)
end
function cf_dot(a, b)
    am, A = a; bm, B = b
    c = am[1] * bm[1] + am[2] * bm[2]
    g = [am[1] * B[1, k] + am[2] * B[2, k] + A[1, k] * bm[1] + A[2, k] * bm[2] for k = 1:4]
    Q = [A[1, k] * B[1, l] + A[2, k] * B[2, l] for k = 1:4, l = 1:4]
    CFQ(c, g, Q)
end
"q1 − η q2 (η a ball); the Q part stays exact only if η q2.Q is exact — so it is bounded as a ball matrix"
function cf_comb_upper(q1::CFQ, q2::CFQ, eta::Arb, r)
    ra = cfarb.(r); s = q1.c - eta * q2.c
    for k = 1:4; s += abs(q1.g[k] - eta * q2.g[k]) * ra[k]; end
    for k = 1:4
        qkk = cfarb(q1.Q[k, k]) - eta * cfarb(q2.Q[k, k])       # δ_k² ∈ [0, r_k²]: contributes ≤ max(q_kk, 0) r_k²
        ub = Arblib.ubound(qkk)
        ub > 0 && (s += Arb(ub; prec = PREC) * ra[k]^2)
        for l = k+1:4
            s += abs(cfarb(q1.Q[k, l] + q1.Q[l, k]) - eta * cfarb(q2.Q[k, l] + q2.Q[l, k])) * ra[k] * ra[l]
        end
    end
    return Arblib.ubound(s)
end
cf_upper(q::CFQ, r) = cf_comb_upper(q, CFQ(cfarb(0), [cfarb(0) for _ = 1:4], zeros(4, 4)), cfarb(0), r)
"enclosures over the box: value ball, gradient balls, constant Hessian S = Q + Qᵀ"
function cf_box(q::CFQ, r)
    ra = cfarb.(r); w = cfarb(0)
    for k = 1:4; w += abs(q.g[k]) * ra[k]; end
    for k = 1:4, l = 1:4; w += abs(cfarb(q.Q[k, l])) * ra[k] * ra[l]; end
    val = q.c + Arb((-1, 1); prec = PREC) * w
    S = q.Q .+ q.Q'
    grad = [q.g[k] + Arb((-1, 1); prec = PREC) * sum(abs(S[k, l]) * ra[l] for l = 1:4) for k = 1:4]
    return val, grad, S
end
"upper bound over the box of F = C² − s2·E·Dd (C, E, Dd quadratics at the same centre)"
function cf_quartic_upper(C::CFQ, E::CFQ, Dd::CFQ, s2::Arb, r)
    ra = cfarb.(r)
    F0 = C.c^2 - s2 * E.c * Dd.c
    s = F0
    for k = 1:4
        s += abs(2 * C.c * C.g[k] - s2 * (E.c * Dd.g[k] + Dd.c * E.g[k])) * ra[k]
    end
    Cb, gC, SC = cf_box(C, r); Eb, gE, SE = cf_box(E, r); Db, gD, SD = cf_box(Dd, r)
    for k = 1:4, l = 1:4
        H = 2 * (Cb * SC[k, l] + gC[k] * gC[l]) - s2 * (Eb * SD[k, l] + Db * SE[k, l] + gE[k] * gD[l] + gD[k] * gE[l])
        s += abs(H) * ra[k] * ra[l] / 2
    end
    return Arblib.ubound(s)
end
"does condition `cond` (a closure (m, r) -> upper bound) hold on [lo, hi]?  adaptive bisection, at most `budget` boxes"
function cf_holds(cond, lo, hi, depth, budget::Ref{Int})
    budget[] -= 1
    budget[] < 0 && return false
    m, r = box_center_radius(lo, hi)
    cond(m, r) <= 0 && return true
    depth == 0 && return false
    d = argmax(hi .- lo); mid = (lo[d] + hi[d]) / 2
    (lo[d] < mid < hi[d]) || return false
    h1 = copy(hi); h1[d] = mid; l2 = copy(lo); l2[d] = mid
    return cf_holds(cond, lo, h1, depth - 1, budget) && cf_holds(cond, l2, hi, depth - 1, budget)
end
const CF_ETA = Arb(numerator(RAT_ETA); prec = PREC) / Arb(denominator(RAT_ETA); prec = PREC)
"independent check of exclusion `kind` on [lo, hi] by centred forms: the proved condition index, or nothing"
function check_E_cf_which(kind, lo, hi; depth = 10, budget = 2000)
    if kind == "E1"
        return all(abs(Rational{BigInt}(lo[i])) <= CUBE && abs(Rational{BigInt}(hi[i])) <= CUBE for i = 1:4) ? (0,) : nothing
    end
    cvq(m, i) = (v = [cf_vert(j, m) for j = 1:4]; cf_cross(cf_sub(v[mod1(i + 1, 4)], v[i]), cf_sub(v[mod1(i + 2, 4)], v[i])))
    if kind == "E4"
        for i = 1:4
            cf_holds((m, r) -> cf_upper(cvq(m, i), r), lo, hi, depth, Ref(budget)) && return (i,)
        end
    elseif kind == "E2"
        for i = 1:4
            cf_holds((m, r) -> cf_comb_upper(cvq(m, i), cvq(m, mod1(i + 2, 4)), CF_ETA, r), lo, hi, depth, Ref(budget)) && return (i,)
        end
    elseif kind == "E3" || kind == "E3T"
        srat = kind == "E3" ? RAT_S0 : RAT_S1
        s2 = Arb(numerator(srat^2); prec = PREC) / Arb(denominator(srat^2); prec = PREC)
        for i = 1:4, kk = 1:4, ll = kk+1:4
            cond = function (m, r)
                v = [cf_vert(j, m) for j = 1:4]
                e = cf_sub(v[mod1(i + 1, 4)], v[i]); E = cf_dot(e, e)
                dkl = cf_sub(v[kk], v[ll]); Dd = cf_dot(dkl, dkl)
                u = Arf(-Inf)
                for j in (mod1(i + 2, 4), mod1(i + 3, 4))
                    C = cf_cross(e, cf_sub(v[j], v[i]))
                    x = cf_quartic_upper(C, E, Dd, s2, r); x > u && (u = x)
                end
                u
            end
            cf_holds(cond, lo, hi, depth, Ref(budget)) && return (i, kk, ll)
        end
    end
    return nothing
end
check_E_cf(kind, lo, hi; kw...) = check_E_cf_which(kind, lo, hi; kw...) !== nothing
