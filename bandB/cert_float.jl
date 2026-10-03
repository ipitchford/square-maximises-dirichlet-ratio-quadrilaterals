# cert_float.jl — floating-point layer of the box certificate (NOT rigorous; it only chooses E, the mesh and the
# Ritz test functions f1, f2, and estimates losses for the bisection heuristic).  Author: Anonymous.
using LinearAlgebra

const C6F = 1 / 6 + 1 / (4 * pi^2)

# ---------------- meshes ----------------
"""
    graded_nodes(X, Y; delta, p, dxmax, hmax) -> Vector{Float64}

Nodes on [0,1] for a PWL profile with breakpoints X (floats, X[1]=0, X[end]=1) and values Y ≥ 0.  On each linear piece
the nodes equidistribute the weight w(h) = |d log h| · (h/hmax)^p (cell 'loss' ∝ relative variation of h, weighted
by how much the eigenfunction can live there), so that each cell carries weight ≤ delta; cell length ≤ dxmax.
Breakpoints are nodes.
"""
function graded_nodes(X::Vector{Float64}, Y::Vector{Float64}; delta = 0.004, p = 2.0, dxmax = 0.01, hmax = maximum(Y))
    xs = Float64[0.0]
    for j = 1:length(X)-1
        a, b = X[j], X[j+1]; ya, yb = Y[j], Y[j+1]
        b > a || continue
        # weight of the piece: |ya^p - yb^p| / (p hmax^p)
        wa = (ya / hmax)^p; wb = (yb / hmax)^p
        nw = ceil(Int, abs(wa - wb) / (p * delta))
        nl = ceil(Int, (b - a) / dxmax)
        n = max(nw, nl, 1)
        if nw >= nl && ya != yb
            # nodes at equal increments of h^p along the piece (h linear in x)
            for i = 1:n-1
                w = wa + (wb - wa) * i / n
                h = hmax * w^(1 / p)
                push!(xs, a + (b - a) * (h - ya) / (yb - ya))
            end
        else
            for i = 1:n-1
                push!(xs, a + (b - a) * i / n)
            end
        end
        push!(xs, b)
    end
    xs = sort(unique(xs))
    xs[1] = 0.0; xs[end] = 1.0
    return xs
end

# ---------------- step-potential minorant: float Sturm shooting ----------------
"float version of the rigorous Sturm test (same cell rules); returns (ok, v(1))"
function shoot_f(xs::Vector{Float64}, V::Vector{Float64}, E::Float64)
    v = 0.0; dv = 1.0
    for k in eachindex(V)
        L = xs[k+1] - xs[k]
        q = V[k] - E
        if q > 0
            κ = sqrt(q); c = cosh(κ * L); s = sinh(κ * L)
            vn = v * c + dv * s / κ; dvn = v * κ * s + dv * c
        elseif q < 0
            ω = sqrt(-q)
            ω * L < pi || return (false, NaN)
            c = cos(ω * L); s = sin(ω * L)
            vn = v * c + dv * s / ω; dvn = -v * ω * s + dv * c
        else
            vn = v + dv * L; dvn = dv
        end
        vn > 0 || return (false, vn)
        v, dv = vn, dvn
        m = max(abs(v), abs(dv))
        if m > 1e150 || m < 1e-150
            v /= m; dv /= m
        end
    end
    return (true, v)
end

"largest E with no zero (≈ μ1 of the step potential), by bisection to relative tolerance rtol"
function mu1_step_f(xs, V; rtol = 1e-13)
    lo = minimum(V) * (1 - 1e-12)
    shoot_f(xs, V, lo)[1] || return 0.0
    hi = max(lo * 2, 1.0)
    while shoot_f(xs, V, hi)[1]
        lo = hi; hi *= 2
        hi > 1e12 && return lo
    end
    while hi - lo > rtol * hi
        mid = (lo + hi) / 2
        shoot_f(xs, V, mid)[1] ? (lo = mid) : (hi = mid)
    end
    return lo
end

# ---------------- weighted 1-D Ritz form (upper bound) ----------------
# element integrals I00 = ∫(1−s)²/h, I01 = ∫s(1−s)/h, I11 = ∫s²/h over s∈[0,1], h = ha(1−s) + hb s, ha, hb > 0
const SER_C = [(n % 2 == 0 ? 0.5^n / (n + 1) : 0.0) / 4 + (iseven(n + 1) ? 0.5^(n + 1) / (n + 2) : 0.0) +
               (iseven(n + 2) ? 0.5^(n + 2) / (n + 3) : 0.0) for n = 0:80]
const SER_B = [(iseven(n) ? 0.5^n / (n + 1) : 0.0) / 4 - (iseven(n + 2) ? 0.5^(n + 2) / (n + 3) : 0.0) for n = 0:80]
function Iints_f(ha::Float64, hb::Float64)
    hm = (ha + hb) / 2; ε = (hb - ha) / hm; ρ = abs(ε) / 2
    if ρ <= 0.25
        s00 = 0.0; s11 = 0.0; s01 = 0.0; pw = 1.0
        for n = 0:60
            s00 += pw * SER_C[n+1]; s11 += (iseven(n) ? pw : -pw) * SER_C[n+1]; s01 += pw * SER_B[n+1]
            pw *= ε
            abs(pw) < 1e-18 && break
        end
        return s00 / hm, s01 / hm, s11 / hm
    else
        Δ = hb - ha; lg = log(hb / ha)
        I11 = (Δ * (hb - 3ha) / 2 + ha^2 * lg) / Δ^3
        I00 = (Δ * (ha - 3hb) / 2 + hb^2 * lg) / Δ^3
        I01 = (Δ - ha * lg) / Δ^2 - I11
        return I00, I01, I11
    end
end

"""
    ritz_matrices_f(xs, a, b) -> (Kd, Ke, Md, Me) tridiagonal (interior nodes 2..n-1)

P1 discretisation of N(f) = ∫ f'²h/2 + f f' h'/2 + f²(π²/h)(1/2 + a'b'/2 + h'²C6),  D(f) = ∫ f² h/2.
"""
function ritz_matrices_f(xs, a, b)
    n = length(xs)
    Kd = zeros(n); Ke = zeros(n - 1); Md = zeros(n); Me = zeros(n - 1)
    for e = 1:n-1
        L = xs[e+1] - xs[e]
        ha = b[e] - a[e]; hb = b[e+1] - a[e+1]
        ap = (a[e+1] - a[e]) / L; bp = (b[e+1] - b[e]) / L; hp = bp - ap
        ce = 0.5 + ap * bp / 2 + hp^2 * C6F
        k1 = (ha + hb) / (4L)
        if e == 1
            Kd[2] += k1 + hp / 4 + pi^2 * ce * L / (2hb); Md[2] += L * 3hb / 24
        elseif e == n - 1
            Kd[n-1] += k1 - hp / 4 + pi^2 * ce * L / (2ha); Md[n-1] += L * 3ha / 24
        else
            I00, I01, I11 = Iints_f(ha, hb)
            Kd[e] += k1 - hp / 4 + pi^2 * ce * L * I00
            Kd[e+1] += k1 + hp / 4 + pi^2 * ce * L * I11
            Ke[e] += -k1 + pi^2 * ce * L * I01
            Md[e] += L * (3ha + hb) / 24; Md[e+1] += L * (ha + 3hb) / 24; Me[e] += L * (ha + hb) / 24
        end
    end
    return Kd[2:n-1], Ke[2:n-2], Md[2:n-1], Me[2:n-2]
end

tri_mul(d, e, x) = [d[i] * x[i] + (i > 1 ? e[i-1] * x[i-1] : 0.0) + (i < length(d) ? e[i] * x[i+1] : 0.0) for i in eachindex(d)]

"two lowest Ritz pairs (lumped-mass eigenvectors + `iters` steps of block inverse iteration, then 2×2 Rayleigh-Ritz)"
function ritz2_f(xs, a, b; iters = 1)
    Kd, Ke, Md, Me = ritz_matrices_f(xs, a, b)
    m = length(Kd)
    ml = copy(Md); for i = 1:m-1; ml[i] += Me[i]; ml[i+1] += Me[i]; end
    sq = sqrt.(ml)
    A = SymTridiagonal(Kd ./ ml, Ke ./ (sq[1:end-1] .* sq[2:end]))
    F = eigen(A, 1:2)
    X = F.vectors ./ sq
    Kt = SymTridiagonal(Kd, Ke)
    for _ = 1:iters
        Y = hcat(tri_mul(Md, Me, X[:, 1]), tri_mul(Md, Me, X[:, 2]))
        X = Kt \ Y
        X ./= maximum(abs, X; dims = 1)
    end
    KX = hcat(tri_mul(Kd, Ke, X[:, 1]), tri_mul(Kd, Ke, X[:, 2])); MX = hcat(tri_mul(Md, Me, X[:, 1]), tri_mul(Md, Me, X[:, 2]))
    N2 = Symmetric(X' * KX); D2 = Symmetric(X' * MX)
    G = eigen(N2, D2)
    f = X * G.vectors
    f ./= maximum(abs, f; dims = 1)
    fz = vcat(zeros(1, 2), f, zeros(1, 2))
    return G.values, fz
end
