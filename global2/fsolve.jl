# fsolve.jl — Float64 (NON-rigorous) solver at a single parameter point p0.  Used only to CHOOSE the
# test functions (Rayleigh–Ritz eigenvectors), rho0, the stream-function polynomial psi, and as the pilot.
#
# Basis: b(x)·P_i(u)P_l(v), i+l ≤ n, u = (x−xc)/sx, v = (y−yc)/sy, (xc,yc,sx,sy) = bounding box of Q_{p0}.
# Flux for Lehmann–Goerisch (see README_global2.md):  w = ∇u1/ρ0 + τ + curl ψ,  τ = φ·e,  e = v3 − v1,
#   φ(x) = −(1/ρ0) ∫_0^{Z(x)} r(P0 + X(P1−P0) + Z'(P2−P0)) dZ'   on each triangle (P0,P1,P2) = (v1, v2|v4, v3),
#   r = Δu1 + ρ0 u1.  Then e·∇φ = −r/ρ0, so div w = Δu1/ρ0 − r/ρ0 = −u1 exactly; φ is continuous across the
#   diagonal (both sides integrate along the diagonal from v1), so τ ∈ H(div, Q).

using LinearAlgebra

struct Scaling
    xc::Float64; yc::Float64; sx::Float64; sy::Float64
end
function basis_scaling(V::AbstractMatrix)
    xlo, xhi = extrema(V[:, 1]); ylo, yhi = extrema(V[:, 2])
    Scaling((xlo + xhi) / 2, (ylo + yhi) / 2, (xhi - xlo) / 2, (yhi - ylo) / 2)
end

basis_index(n::Int) = [(i, l) for i = 0:n for l = 0:n-i]
psi_index(k::Int) = [(i, l) for i = 0:k for l = 0:k-i if (i, l) != (0, 0)]

# Legendre P_k, P_k', P_k'' (k = 0..K) at u (Float64)
function legendre012!(P, dP, d2P, u::Float64, K::Int)
    P[1] = 1.0; dP[1] = 0.0; d2P[1] = 0.0
    K == 0 && return
    P[2] = u; dP[2] = 1.0; d2P[2] = 0.0
    for k = 1:K-1
        a = 2k + 1
        P[k+2] = (a * u * P[k+1] - k * P[k]) / (k + 1)
        dP[k+2] = dP[k] + a * P[k+1]
        d2P[k+2] = d2P[k] + a * dP[k+1]
    end
end

# Float64 node set of the split quadrilateral: positions, weights, and bookkeeping
struct NodesF
    x::Vector{Float64}; y::Vector{Float64}; w::Vector{Float64}
    tri::Vector{Int}; i::Vector{Int}; j::Vector{Int}
end
function nodes_f(V::AbstractMatrix, R::RuleF)
    m = R.m; xs = Float64[]; ys = Float64[]; ws = Float64[]; tr = Int[]; ii = Int[]; jj = Int[]
    for (t, (i0, i1, i2, sg)) in enumerate(TRIS)
        P0 = V[i0, :]; P1 = V[i1, :]; P2 = V[i2, :]
        J = sg * cross2(P1[1]-P0[1], P1[2]-P0[2], P2[1]-P0[1], P2[2]-P0[2])
        for i = 1:m, j = 1:m
            X = R.X[i]; Y = R.X[j]
            push!(xs, P0[1] + X * (P1[1]-P0[1]) + Y * (1-X) * (P2[1]-P0[1]))
            push!(ys, P0[2] + X * (P1[2]-P0[2]) + Y * (1-X) * (P2[2]-P0[2]))
            push!(ws, R.W[i] * R.W[j] * (1 - X) * J)
            push!(tr, t); push!(ii, i); push!(jj, j)
        end
    end
    NodesF(xs, ys, ws, tr, ii, jj)
end

# bubble b, ∇b, Δb at (x,y) (Float64)
function bubble_f(V::AbstractMatrix, x, y)
    l = zeros(4); gx = zeros(4); gy = zeros(4)
    for k = 1:4
        k2 = mod1(k + 1, 4)
        ex = V[k2, 1] - V[k, 1]; ey = V[k2, 2] - V[k, 2]
        gx[k] = -ey; gy[k] = ex
        l[k] = gx[k] * (x - V[k, 1]) + gy[k] * (y - V[k, 2])
    end
    b = l[1] * l[2] * l[3] * l[4]
    bx = 0.0; by = 0.0; lap = 0.0
    for k = 1:4
        pr = 1.0; for j = 1:4; j != k && (pr *= l[j]); end
        bx += gx[k] * pr; by += gy[k] * pr
    end
    for i = 1:4, j = i+1:4
        pr = 1.0; for k = 1:4; (k != i && k != j) && (pr *= l[k]); end
        lap += 2 * (gx[i] * gx[j] + gy[i] * gy[j]) * pr
    end
    return b, bx, by, lap
end

"""
    fsolve(p, n; kpsi = n+6, nfun = 3) -> NamedTuple

Float64 Rayleigh–Ritz in the degree-n polynomial×bubble basis at Q_p, the first `nfun` eigenvectors
(M-orthonormal), ρ0 = θ1, and the LG flux (τ along the diagonal + least-squares curl ψ, deg ψ ≤ kpsi).
A1 = 1 by normalisation; A0 = θ1; A2 = ‖w‖².
"""
fsolve(p::AbstractVector{<:Real}, n::Int; kw...) = fsolve_V(verts(p), n; p = collect(Float64, p), kw...)
function fsolve_V(V::AbstractMatrix, n::Int; kpsi::Int = n + 6, nfun::Int = 3, p = Float64[],
                  sc::Union{Nothing,Scaling} = nothing, rho_fix::Union{Nothing,Float64} = nothing,
                  Cref::Union{Nothing,Matrix{Float64}} = nothing, Mref::Union{Nothing,Matrix{Float64}} = nothing)
    sc === nothing && (sc = basis_scaling(V))
    m = npoints(n); R = rule_f(m); nd = nodes_f(V, R); Q = length(nd.x)
    K = max(n, kpsi)
    idx = basis_index(n); N = length(idx)
    pidx = psi_index(kpsi); Np = length(pidx)
    Pu = zeros(K + 1); dPu = zeros(K + 1); d2Pu = zeros(K + 1)
    Pv = zeros(K + 1); dPv = zeros(K + 1); d2Pv = zeros(K + 1)
    Phi = zeros(N, Q); Gx = zeros(N, Q); Gy = zeros(N, Q); Lap = zeros(N, Q)
    Cx = zeros(Np, Q); Cy = zeros(Np, Q)
    isx = 1 / sc.sx; isy = 1 / sc.sy
    for q = 1:Q
        x = nd.x[q]; y = nd.y[q]
        legendre012!(Pu, dPu, d2Pu, (x - sc.xc) * isx, K)
        legendre012!(Pv, dPv, d2Pv, (y - sc.yc) * isy, K)
        b, bx, by, lapb = bubble_f(V, x, y)
        for (k, (i, l)) in enumerate(idx)
            P = Pu[i+1] * Pv[l+1]; Px = dPu[i+1] * Pv[l+1] * isx; Py = Pu[i+1] * dPv[l+1] * isy
            LP = d2Pu[i+1] * Pv[l+1] * isx^2 + Pu[i+1] * d2Pv[l+1] * isy^2
            Phi[k, q] = b * P; Gx[k, q] = bx * P + b * Px; Gy[k, q] = by * P + b * Py
            Lap[k, q] = lapb * P + 2 * (bx * Px + by * Py) + b * LP
        end
        for (k, (i, l)) in enumerate(pidx)
            Cx[k, q] = Pu[i+1] * dPv[l+1] * isy          # ∂ψ/∂y
            Cy[k, q] = -dPu[i+1] * Pv[l+1] * isx         # −∂ψ/∂x
        end
    end
    W = nd.w
    M = (Phi .* W') * Phi'; Kst = (Gx .* W') * Gx' + (Gy .* W') * Gy'
    M = (M + M') / 2; Kst = (Kst + Kst') / 2
    theta, U = robust_geneig(Kst, M)
    C = U[:, 1:nfun]
    if Cref !== nothing
        # align with the reference vectors (in the reference M-inner product): sign of c1, Procrustes for (c2,c3)
        G = Cref' * Mref * C
        G[1, 1] < 0 && (C[:, 1] .*= -1)
        F2 = svd(G[2:3, 2:3])              # maximise tr(Q' G23) over orthogonal Q
        Q = F2.V * F2.U'
        C[:, 2:3] = C[:, 2:3] * Q
    end
    rho0 = rho_fix === nothing ? theta[1] : rho_fix
    c1 = C[:, 1]
    u1 = Phi' * c1; ux = Gx' * c1; uy = Gy' * c1; lap1 = Lap' * c1
    r = lap1 .+ rho0 .* u1
    Ystar = ystar_f(V)
    phi = line_integrate_f(r, nd, R, rho0, Ystar)
    ex = V[3, 1] - V[1, 1]; ey = V[3, 2] - V[1, 2]
    gx = ux ./ rho0 .+ phi .* ex; gy = uy ./ rho0 .+ phi .* ey
    sw = sqrt.(abs.(W))
    A = vcat((Cx .* sw')', (Cy .* sw')')
    rhs = -vcat(sw .* gx, sw .* gy)
    # Tikhonov-regularised least squares (smooth dependence on p; μ tiny relative to the scale)
    mu = 1e-10 * opnorm(A, 1)
    dpsi = vcat(A, mu * Matrix(I, size(A, 2), size(A, 2))) \ vcat(rhs, zeros(size(A, 2)))
    wx = gx .+ Cx' * dpsi; wy = gy .+ Cy' * dpsi
    A0 = dot(W, ux .^ 2 .+ uy .^ 2); A1 = dot(W, u1 .^ 2); A2 = dot(W, wx .^ 2 .+ wy .^ 2)
    A2_nopsi = dot(W, gx .^ 2 .+ gy .^ 2)
    eps2 = dot(W, r .^ 2) / A1
    return (; p, n, kpsi, V, sc, theta, C, rho0, dpsi, A0, A1, A2, A2_nopsi, eps2, N, Np, condM = cond(M), M)
end

# Global polynomial φ with φ = 0 on the line v2v4 (the other diagonal).  In Duffy coordinates of both
# triangles that line is Y = Y* (Y* = fraction of v1v3 at which the diagonals cross), so on each Duffy line
#   φ_j = −(1/ρ0)(1−X_i) Σ_k (S[j,k] − I_k(Y*)) r_k,   I_k(Y) = ∫_0^Y L_k.
"Y* = parameter along v1→v3 of the intersection of the diagonals (Float64)."
function ystar_f(V)
    e = (V[3,1]-V[1,1], V[3,2]-V[1,2]); f = (V[4,1]-V[2,1], V[4,2]-V[2,2])
    return cross2(f[1], f[2], V[2,1]-V[1,1], V[2,2]-V[1,2]) / cross2(f[1], f[2], e[1], e[2])
end
"I_k(Y) = ∫_0^Y L_k (Float64), k = 1..m, by the m-point Gauss rule on [0,Y]."
function Ivec_f(R::RuleF, Y::Float64)
    m = R.m; I = zeros(m)
    for k = 1:m, g = 1:m
        y = Y * R.X[g]; L = 1.0
        for i = 1:m; i != k && (L *= (y - R.X[i]) / (R.X[k] - R.X[i])); end
        I[k] += R.W[g] * L
    end
    return Y .* I
end
function line_integrate_f(r::AbstractVector, nd::NodesF, R::RuleF, rho0::Float64, Ystar::Float64)
    m = R.m; phi = similar(r); I = Ivec_f(R, Ystar)
    for t = 1:2, i = 1:m
        base = (t - 1) * m * m + (i - 1) * m
        fac = -(1 - R.X[i]) / rho0
        pre = 0.0
        for k = 1:m; pre += I[k] * r[base + k]; end
        for j = 1:m
            acc = 0.0
            for k = 1:m
                acc += R.S[j, k] * r[base + k]
            end
            phi[base + j] = fac * (acc - pre)
        end
    end
    return phi
end

"generalized eigenproblem K c = θ M c with M possibly numerically singular: drop M-directions with
eigenvalue < 1e-13·max (the basis is nearly dependent for thin shapes at high degree)."
function robust_geneig(K::AbstractMatrix, M::AbstractMatrix)
    Fm = eigen(Symmetric(M))
    keep = Fm.values .> 1e-13 * maximum(Fm.values)
    B = Fm.vectors[:, keep] ./ sqrt.(Fm.values[keep])'
    F = eigen(Symmetric(B' * K * B))
    return F.values, B * F.vectors
end

# ---------------- Float64 a-priori lambda_2 lower bounds at a point (pilot only) ----------------
const J01 = 2.4048255576957728
beta_ks_f(p) = 2π * J01^2 / area(p)

"rectangle with sides along the unit direction (c,s): extents A (along) and B (across)"
function rect_extents_f(V, c, s)
    pr = [V[j, 1] * c + V[j, 2] * s for j = 1:4]
    qr_ = [-V[j, 1] * s + V[j, 2] * c for j = 1:4]
    return maximum(pr) - minimum(pr), maximum(qr_) - minimum(qr_)
end
lam2_rect(A, B) = π^2 * min(4 / A^2 + 1 / B^2, 1 / A^2 + 4 / B^2)
lam1_rect(A, B) = π^2 * (1 / A^2 + 1 / B^2)

"candidate directions at V: axis, the four edges, the two diagonals (angles in radians)"
function candidate_dirs(V)
    th = [0.0]
    for i = 1:4
        j = mod1(i + 1, 4); push!(th, atan(V[j, 2] - V[i, 2], V[j, 1] - V[i, 1]))
    end
    push!(th, atan(V[3, 2] - V[1, 2], V[3, 1] - V[1, 1])); push!(th, atan(V[4, 2] - V[2, 2], V[4, 1] - V[2, 1]))
    return th
end

function beta_point_f(p)
    V = verts(p)
    best = beta_ks_f(p)
    for th in candidate_dirs(V)
        A, B = rect_extents_f(V, cos(th), sin(th)); best = max(best, lam2_rect(A, B))
    end
    return best
end
function lmono_point_f(p)
    V = verts(p); best = 0.0
    for th in candidate_dirs(V)
        A, B = rect_extents_f(V, cos(th), sin(th)); best = max(best, lam1_rect(A, B))
    end
    return best
end

"LG lower bound for λ1 given β ≤ λ2 (requires βA1 > A0)."
lg_bound(A0, A1, A2, beta) = (beta * A1 - A0) / (beta * A2 - A1)

"""
    fsolve_d(p0, n; kpsi, h) -> NamedTuple

Float64 solve at p0 plus first-order p-derivatives (central differences, step h, same basis scaling,
fixed ρ0, eigenvectors aligned to those at p0) of the three test-polynomial coefficient vectors and of the
stream-function coefficients.  Used only to CHOOSE the p-dependent test functions
P_k(x;p) = P_k0(x) + Σ_i (p_i − p0_i) P_ki(x),  ψ(x;p) = ψ0 + Σ_i (p_i − p0_i) ψ_i.
"""
function fsolve_d(p0::AbstractVector{<:Real}, n::Int; kpsi::Int = n + 6, h::Float64 = 1e-4)
    p0 = collect(Float64, p0)
    s0 = fsolve(p0, n; kpsi)
    dC = zeros(size(s0.C, 1), size(s0.C, 2), 4); dpsi = zeros(length(s0.dpsi), 4)
    for i = 1:4
        e = zeros(4); e[i] = h
        sp = fsolve(p0 .+ e, n; kpsi, sc = s0.sc, rho_fix = s0.rho0, Cref = s0.C, Mref = s0.M)
        sm = fsolve(p0 .- e, n; kpsi, sc = s0.sc, rho_fix = s0.rho0, Cref = s0.C, Mref = s0.M)
        dC[:, :, i] = (sp.C .- sm.C) ./ (2h)
        dpsi[:, i] = (sp.dpsi .- sm.dpsi) ./ (2h)
    end
    return merge(s0, (; dC, ddpsi = dpsi, p0 = p0))
end
