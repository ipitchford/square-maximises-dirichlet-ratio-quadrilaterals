# transport.jl — the affine-transport certificate (TRANSPORT_SPEC.md).
#
# Reference quadrilateral Q0 = Q_{p0}, split along v1v3 into K1 = (v1,v2,v3) and K2 = (v1,v4,v3) (= TRIS).
# For p in the box, Ψ_p is affine on each K, mapping K(p0) → K(p):  x = v1(p) + A_K(p)(x̂ − v1(p0)),
# A_K(p) = M_K(p) M_K(p0)⁻¹,  M_K(p) = [v_{i1}(p) − v1(p), v3(p) − v1(p)]  (affine in p; A_K(p0) = I).
# The two maps agree on the diagonal v1v3, so Ψ_p is continuous and maps ∂Q0 onto ∂Q_p.
#
# Test functions  u_k(p) = û_k ∘ Ψ_p⁻¹,  û_k = b_{Q0}·P_k  (P_k = first three Rayleigh–Ritz eigenvectors at p0):
#   continuous, piecewise polynomial, zero on ∂Q_p  ⇒  u_k(p) ∈ H¹₀(Q_p).
#   K_ij(p) = Σ_K (1/d_K) Σ_{ab} (adj A_K adj A_Kᵀ)_{ab} I^K_{ab,ij},  I^K_{ab,ij} = ∫_{K(p0)} ∂_a û_i ∂_b û_j,
#   M_ij(p) = Σ_K d_K M^K_ij,  M^K_ij = ∫_{K(p0)} û_i û_j,   d_K = det A_K.
# Flux (Lehmann–Goerisch): reference ŵ = d_K τ̂ + curl ψ̂ on K(p0), physical w = A_K ŵ / d_K (contravariant Piola):
#   τ̂ = φ̂·e, e = v3(p0) − v1(p0),  φ̂(x̂) = −∫_{σ}^{0} û_1(x̂ + s e) ds with x̂ + σe on the line v2v4 (of Q0):
#   e·∇φ̂ = −û_1, so div τ̂ = −û_1;  τ̂ ∥ diagonal, so τ̂·n = 0 on v1v3;  ψ̂ a fixed polynomial.
#   div_x w = (1/d_K) div_x̂ ŵ = div τ̂ = −û_1 = −u_1(p) on each K; the normal flux across the diagonal is
#   curl ψ̂·n̂ dŝ from both sides (Piola preserves normal flux), so w ∈ H(div, Q_p) and div w + u_1(p) = 0 exactly.
#   A2(p) = ‖w‖² = Σ_K [ d_K ⟨G_K, TT^K⟩ + 2⟨G_K, TC^K⟩ + ⟨G_K, CC^K⟩/d_K ],  G_K = A_KᵀA_K,
#   TT^K_ab = ∫τ̂_aτ̂_b, TC^K_ab = ∫τ̂_a(curl ψ̂)_b, CC^K_ab = ∫(curl ψ̂)_a(curl ψ̂)_b over K(p0).
# All reference integrals: exact Duffy–Gauss quadrature on K(p0) (degrees ≤ 2n+10 with deg ψ̂ = n+6, m = n+6).
# Cleared denominators: the criterion is homogeneous in (K, M, A2), so we use D·K, D·M, D·A2 with
# D = d_1 d_2 > 0 (certified), which are polynomials in δ = p − p0:
#   D·K_ij = d_2 Σ(adj₁adj₁ᵀ)·I¹ + d_1 Σ(adj₂adj₂ᵀ)·I²,
#   D·M_ij = D (d_1 M¹ + d_2 M²),
#   D·A2  = D Σ_K (d_K⟨G_K,TT^K⟩ + 2⟨G_K,TC^K⟩) + d_2⟨G_1,CC¹⟩ + d_1⟨G_2,CC²⟩.

# ------------------------------------------------------------------ Float64 choice at p0
"""
    tsolve(p0, n; kpsi = n+6) -> NamedTuple

Float64 (NOT rigorous; it only chooses the test functions): Rayleigh–Ritz at Q_{p0} in the degree-n
polynomial×bubble basis (eigenvectors M-orthonormal), the reference flux τ̂ at the quadrature nodes, and the
Tikhonov-regularised least-squares ψ̂ (deg ≤ kpsi) minimising ‖τ̂ + curl ψ̂‖² on Q_{p0}.
"""
# Diagonal choice.  The split is along V'1V'3 of the relabelled quadrilateral V'_j = v_{perm[j]},
# perm = (1+s, 2+s, 3+s, 4+s) mod 4, s ∈ {0, 1}: s = 0 splits along v1v3, s = 1 along v2v4.  A vertex with an angle
# close to π must be an END of the diagonal (as the apex of a split triangle it makes that triangle flat and A_K
# badly conditioned over the box), so we choose s such that the vertex of largest angle lies on the diagonal.
vperm(s::Int) = [mod1(j + s, 4) for j = 1:4]
choose_shift(p0) = (a = angles_deg(verts(p0)); argmax(a) in (1, 3) ? 0 : 1)

function tsolve(p0::AbstractVector{<:Real}, n::Int; kpsi::Int = n + 6, shift::Int = choose_shift(p0))
    p0 = collect(Float64, p0)
    V = verts(p0)[vperm(shift), :]; sc = basis_scaling(V)
    m = npoints(n); R = rule_f(m); nd = nodes_f(V, R); Q = length(nd.x)
    Kd = max(n, kpsi)
    idx = basis_index(n); N = length(idx); pidx = psi_index(kpsi); Np = length(pidx)
    Pu = zeros(Kd + 1); dPu = zeros(Kd + 1); d2Pu = zeros(Kd + 1); Pv = zeros(Kd + 1); dPv = zeros(Kd + 1); d2Pv = zeros(Kd + 1)
    Phi = zeros(N, Q); Gx = zeros(N, Q); Gy = zeros(N, Q); Cx = zeros(Np, Q); Cy = zeros(Np, Q)
    isx = 1 / sc.sx; isy = 1 / sc.sy
    for q = 1:Q
        x = nd.x[q]; y = nd.y[q]
        legendre012!(Pu, dPu, d2Pu, (x - sc.xc) * isx, Kd); legendre012!(Pv, dPv, d2Pv, (y - sc.yc) * isy, Kd)
        b, bx, by, _ = bubble_f(V, x, y)
        for (k, (i, l)) in enumerate(idx)
            P = Pu[i+1] * Pv[l+1]
            Phi[k, q] = b * P; Gx[k, q] = bx * P + b * dPu[i+1] * Pv[l+1] * isx; Gy[k, q] = by * P + b * Pu[i+1] * dPv[l+1] * isy
        end
        for (k, (i, l)) in enumerate(pidx)
            Cx[k, q] = Pu[i+1] * dPv[l+1] * isy; Cy[k, q] = -dPu[i+1] * Pv[l+1] * isx
        end
    end
    W = nd.w
    M = Symmetric((Phi .* W') * Phi'); Kst = Symmetric((Gx .* W') * Gx' + (Gy .* W') * Gy')
    theta, U = robust_geneig(Matrix(Kst), Matrix(M))
    C = U[:, 1:3]
    u1 = Phi' * C[:, 1]
    Ys = ystar_f(V); Iv = Ivec_f(R, Ys)
    phi = zeros(Q)
    for t = 1:2, i = 1:m
        base = (t - 1) * m * m + (i - 1) * m
        pre = sum(Iv[k] * u1[base+k] for k = 1:m)
        for j = 1:m
            acc = sum(R.S[j, k] * u1[base+k] for k = 1:m)
            phi[base+j] = -(1 - R.X[i]) * (acc - pre)
        end
    end
    ex = V[3, 1] - V[1, 1]; ey = V[3, 2] - V[1, 2]
    tx = phi .* ex; ty = phi .* ey
    sw = sqrt.(abs.(W))
    A = vcat((Cx .* sw')', (Cy .* sw')'); rhs = -vcat(sw .* tx, sw .* ty)
    mu = 1e-10 * opnorm(A, 1)
    dpsi = vcat(A, mu * Matrix(I, Np, Np)) \ vcat(rhs, zeros(Np))
    wx = tx .+ Cx' * dpsi; wy = ty .+ Cy' * dpsi
    A2 = dot(W, wx .^ 2 .+ wy .^ 2)
    return (; p0, n, kpsi, shift, V, sc, theta, C, dpsi, A2pt = A2, A0pt = theta[1], A1pt = 1.0)
end

# ------------------------------------------------------------------ reference integrals (rigorous in Arb)
struct TRefs{T}
    I::Vector{Array{T,4}}        # I[K][a,b,i,j]
    M::Vector{Matrix{T}}         # M[K][i,j]
    TT::Vector{Matrix{T}}; TC::Vector{Matrix{T}}; CC::Vector{Matrix{T}}
    Minv0::Vector{Matrix{T}}     # M_K(p0)⁻¹
end

# P, ∂uP, ∂vP of Σ C[i+1,l+1] P_i(u)P_l(v) from derivative tables PU, PV (column 1 = value, 2 = first derivative)
function poly_vg(C::Matrix{T}, PU::Matrix{T}, PV::Matrix{T}, d::Int) where {T}
    p = zero_(T); pu = zero_(T); pv = zero_(T)
    for i = 0:d
        cv = zero_(T); dv = zero_(T)
        for l = 0:d-i
            cv += C[i+1, l+1] * PV[l+1, 1]; dv += C[i+1, l+1] * PV[l+1, 2]
        end
        p += PU[i+1, 1] * cv; pu += PU[i+1, 2] * cv; pv += PU[i+1, 1] * dv
    end
    return p, pu, pv
end
zero_(::Type{Float64}) = 0.0
zero_(::Type{Arb}) = Arb(0; prec = PREC)

"exact vertex coordinates at p0 in type T"
function verts_T(::Type{T}, p0) where {T}
    [begin
         c = tconv(T, V0[j, k])
         for i = 1:4; c += tconv(T, DV[j, k, i]) * tconv(T, p0[i]); end
         c
     end for j = 1:4, k = 1:2]
end

"""
    transport_refs(T, sol) -> TRefs{T}

Reference integrals over K1(p0), K2(p0) for the test functions and flux chosen in `sol` (exact quadrature; in Arb
every quantity is a ball containing the exact integral).
"""
function transport_refs(::Type{T}, sol) where {T}
    p0 = sol.p0; n = sol.n; kpsi = sol.kpsi; sc = sol.sc; Kd = max(n, kpsi)
    m = npoints(n); NCs = nodeconsts(T, m)
    V = verts_T(T, p0)[vperm(sol.shift), :]
    isx = tconv(T, 1.0) / tconv(T, sc.sx); isy = tconv(T, 1.0) / tconv(T, sc.sy)
    txc = tconv(T, sc.xc); tyc = tconv(T, sc.yc)
    # edge data at p0: A[k,j] = l_k(v_j), gradients g_k = (−e_ky, e_kx)
    E = [V[mod1(k + 1, 4), c] - V[k, c] for k = 1:4, c = 1:2]
    Al = [zero_(T) for _ = 1:4, _ = 1:4]
    for k = 1:4, j = 1:4
        (j == k || j == mod1(k + 1, 4)) && continue
        Al[k, j] = E[k, 1] * (V[j, 2] - V[k, 2]) - E[k, 2] * (V[j, 1] - V[k, 1])
    end
    gx = [-E[k, 2] for k = 1:4]; gy = [E[k, 1] for k = 1:4]
    ex = V[3, 1] - V[1, 1]; ey = V[3, 2] - V[1, 2]
    # Y* and I_k(Y*) = ∫_0^{Y*} L_k  (Gauss rule on [0, Y*], exact for the degree-(m−1) L_k)
    fx = V[4, 1] - V[2, 1]; fy = V[4, 2] - V[2, 2]
    Ys = (fx * (V[2, 2] - V[1, 2]) - fy * (V[2, 1] - V[1, 1])) / (fx * ey - fy * ex)
    Ra = rule_arb(m)
    Ivec = Vector{T}(undef, m)
    for k = 1:m
        acc = zero_(T)
        for g = 1:m
            y = Ys * tconv(T, Ra.X[g]); L = tconv(T, 1.0)
            for i = 1:m; i != k && (L *= (y - tconv(T, Ra.X[i])) / (tconv(T, Ra.X[k]) - tconv(T, Ra.X[i]))); end
            acc += tconv(T, Ra.W[g]) * L
        end
        Ivec[k] = Ys * acc
    end
    idx = basis_index(n); pidx = psi_index(kpsi)
    Cm = [coeffmat(T, sol.C[:, k], idx, n) for k = 1:3]
    Cpsi = coeffmat(T, sol.dpsi, pidx, kpsi)
    PU = [zero_(T) for _ = 1:Kd+1, _ = 1:Kd+1]; PV = [zero_(T) for _ = 1:Kd+1, _ = 1:Kd+1]
    Iout = [[zero_(T) for _ = 1:2, _ = 1:2, _ = 1:3, _ = 1:3] for _ = 1:2]
    Mout = [[zero_(T) for _ = 1:3, _ = 1:3] for _ = 1:2]
    TT = [[zero_(T) for _ = 1:2, _ = 1:2] for _ = 1:2]; TC = [[zero_(T) for _ = 1:2, _ = 1:2] for _ = 1:2]
    CC = [[zero_(T) for _ = 1:2, _ = 1:2] for _ = 1:2]
    Minv0 = Matrix{T}[]
    u1l = Vector{T}(undef, m); cxl = Vector{T}(undef, m); cyl = Vector{T}(undef, m); wl = Vector{T}(undef, m)
    for (t, (i0, i1, i2, sg)) in enumerate(TRIS)
        # M_K(p0) and its inverse
        Mk = [V[i1, 1] - V[i0, 1] V[i2, 1] - V[i0, 1]; V[i1, 2] - V[i0, 2] V[i2, 2] - V[i0, 2]]
        detM = Mk[1, 1] * Mk[2, 2] - Mk[1, 2] * Mk[2, 1]
        push!(Minv0, [Mk[2, 2] / detM -Mk[1, 2] / detM; -Mk[2, 1] / detM Mk[1, 1] / detM])
        Jabs = tconv(T, sg) * detM                      # |det| (sg = −1 for the cw triangle)
        for i = 1:m
            for j = 1:m
                lam0 = NCs.lam0[i, j]; xi = NCs.xi[i]; eta = NCs.eta[i, j]
                x = lam0 * V[i0, 1] + xi * V[i1, 1] + eta * V[i2, 1]
                y = lam0 * V[i0, 2] + xi * V[i1, 2] + eta * V[i2, 2]
                legendre_table!(PU, (x - txc) * isx, Kd); legendre_table!(PV, (y - tyc) * isy, Kd)
                l = [lam0 * Al[k, i0] + xi * Al[k, i1] + eta * Al[k, i2] for k = 1:4]
                L12 = l[1] * l[2]; L34 = l[3] * l[4]
                b = L12 * L34
                tt = (L34 * l[2], L34 * l[1], L12 * l[4], L12 * l[3])
                bx = gx[1] * tt[1] + gx[2] * tt[2] + gx[3] * tt[3] + gx[4] * tt[4]
                by = gy[1] * tt[1] + gy[2] * tt[2] + gy[3] * tt[3] + gy[4] * tt[4]
                u = Vector{T}(undef, 3); ux = Vector{T}(undef, 3); uy = Vector{T}(undef, 3)
                for k = 1:3
                    P, Pu_, Pv_ = poly_vg(Cm[k], PU, PV, n)
                    u[k] = b * P; ux[k] = bx * P + b * (Pu_ * isx); uy[k] = by * P + b * (Pv_ * isy)
                end
                _, psu, psv = poly_vg(Cpsi, PU, PV, kpsi)
                cx = psv * isy; cy = -(psu * isx)
                w = NCs.w[i, j] * Jabs
                g = ((ux[1], uy[1]), (ux[2], uy[2]), (ux[3], uy[3]))
                for a = 1:2, bb = 1:2, ii = 1:3, jj = ii:3
                    Iout[t][a, bb, ii, jj] += w * g[ii][a] * g[jj][bb]
                end
                for ii = 1:3, jj = ii:3; Mout[t][ii, jj] += w * u[ii] * u[jj]; end
                CC[t][1, 1] += w * cx * cx; CC[t][1, 2] += w * cx * cy; CC[t][2, 2] += w * cy * cy
                u1l[j] = u[1]; cxl[j] = cx; cyl[j] = cy; wl[j] = w
            end
            # line integration along the Duffy line i: φ̂_j = −(1−X_i)(Σ_k S_jk û1_k − Σ_k I_k û1_k)
            pre = zero_(T)
            for k = 1:m; pre += Ivec[k] * u1l[k]; end
            for j = 1:m
                acc = zero_(T)
                for k = 1:m; acc += NCs.S[j, k] * u1l[k]; end
                phi = -(NCs.omx[i] * (acc - pre))
                tx = phi * ex; ty = phi * ey; w = wl[j]
                TT[t][1, 1] += w * tx * tx; TT[t][1, 2] += w * tx * ty; TT[t][2, 2] += w * ty * ty
                TC[t][1, 1] += w * tx * cxl[j]; TC[t][1, 2] += w * tx * cyl[j]
                TC[t][2, 1] += w * ty * cxl[j]; TC[t][2, 2] += w * ty * cyl[j]
            end
        end
        # symmetric completion
        for a = 1:2, bb = 1:2, ii = 1:3, jj = 1:ii-1
            Iout[t][a, bb, ii, jj] = Iout[t][bb, a, jj, ii]
        end
        for ii = 1:3, jj = 1:ii-1; Mout[t][ii, jj] = Mout[t][jj, ii]; end
        TT[t][2, 1] = TT[t][1, 2]; CC[t][2, 1] = CC[t][1, 2]
    end
    return TRefs{T}(Iout, Mout, TT, TC, CC, Minv0)
end

"point values (A = I): K, M, A2 at p0 from the reference integrals"
function tpoint(R::TRefs{Float64})
    K = zeros(3, 3); M = zeros(3, 3); A2 = 0.0
    for t = 1:2
        K .+= R.I[t][1, 1, :, :] .+ R.I[t][2, 2, :, :]; M .+= R.M[t]
        A2 += R.TT[t][1, 1] + R.TT[t][2, 2] + 2 * (R.TC[t][1, 1] + R.TC[t][2, 2]) + R.CC[t][1, 1] + R.CC[t][2, 2]
    end
    return K, M, A2
end

# ------------------------------------------------------------------ the box Taylor models
struct TPool{T,R}
    A::Vector{Matrix{TM{T,R}}}; d::Vector{TM{T,R}}; Gadj::Vector{Matrix{TM{T,R}}}; GA::Vector{Matrix{TM{T,R}}}
    D::TM{T,R}; t1::TM{T,R}; t2::TM{T,R}; t3::TM{T,R}
    M::Vector{TM{T,R}}; K::Vector{TM{T,R}}; A2::TM{T,R}
    buf::Vector{T}
end
function tpool(::Type{T}) where {T}
    R = tmtype(T).parameters[2]
    TPool{T,R}([[newtm(T) for _ = 1:2, _ = 1:2] for _ = 1:2], [newtm(T) for _ = 1:2],
               [[newtm(T) for _ = 1:2, _ = 1:2] for _ = 1:2], [[newtm(T) for _ = 1:2, _ = 1:2] for _ = 1:2],
               newtm(T), newtm(T), newtm(T), newtm(T), [newtm(T) for _ = 1:6], [newtm(T) for _ = 1:6], newtm(T),
               [tconv(T, 0.0) for _ = 1:FTAB[].nc2])
end

"""
    transport_base!(tp, R, p0, rad, ws, shift) -> (; ok, dlo)

Taylor models over the box p0 ± rad of D·K, D·M, D·A2 (D = d_1 d_2) from the reference integrals R.
`dlo` = certified lower bounds of d_1, d_2 on the box (ok requires both > 0: no fold).
"""
function transport_base!(tp::TPool{T}, R::TRefs{T}, p0::Vector{Float64}, rad::Vector{Float64}, ws, shift::Int) where {T}
    G0 = box_geom(T, p0, rad, ws)
    G = (; V = G0.V[vperm(shift), :])
    fm!(h, f, g) = tm_mul_full!(h, f, g, ws, tp.buf)
    for (t, (i0, i1, i2, sg)) in enumerate(TRIS)
        # M_K(p) columns (exact linear TMs): [v_{i1} − v_{i0}, v_{i2} − v_{i0}]
        Mp = [newtm(T) for _ = 1:2, _ = 1:2]
        for c = 1:2
            tm_sub!(Mp[c, 1], G.V[i1, c], G.V[i0, c]); tm_sub!(Mp[c, 2], G.V[i2, c], G.V[i0, c])
        end
        A = tp.A[t]; Mi = R.Minv0[t]
        for a = 1:2, b = 1:2
            tm_zero!(A[a, b])
            tm_axpy!(A[a, b], Mi[1, b], Mp[a, 1], ws); tm_axpy!(A[a, b], Mi[2, b], Mp[a, 2], ws)
        end
        # d = det A, adj A, Gadj = adj adjᵀ, GA = AᵀA
        fm!(tp.d[t], A[1, 1], A[2, 2]); fm!(tp.t1, A[1, 2], A[2, 1]); tm_sub!(tp.d[t], tp.d[t], tp.t1)
        adj = (A[2, 2], A[1, 2], A[2, 1], A[1, 1])        # adj = [[A22, −A12], [−A21, A11]]
        # (adj adjᵀ)_11 = A22² + A12²,  _12 = −A22 A21 − A12 A11,  _22 = A21² + A11²
        Ga = tp.Gadj[t]
        fm!(Ga[1, 1], A[2, 2], A[2, 2]); fm!(tp.t1, A[1, 2], A[1, 2]); tm_add!(Ga[1, 1], Ga[1, 1], tp.t1)
        fm!(Ga[1, 2], A[2, 2], A[2, 1]); fm!(tp.t1, A[1, 2], A[1, 1]); tm_add!(Ga[1, 2], Ga[1, 2], tp.t1)
        tm_scale!(Ga[1, 2], Ga[1, 2], tconv(T, -1.0), ws); tm_set!(Ga[2, 1], Ga[1, 2])
        fm!(Ga[2, 2], A[2, 1], A[2, 1]); fm!(tp.t1, A[1, 1], A[1, 1]); tm_add!(Ga[2, 2], Ga[2, 2], tp.t1)
        GA = tp.GA[t]
        for a = 1:2, b = a:2
            fm!(GA[a, b], A[1, a], A[1, b]); fm!(tp.t1, A[2, a], A[2, b]); tm_add!(GA[a, b], GA[a, b], tp.t1)
        end
        tm_set!(GA[2, 1], GA[1, 2])
    end
    dlo = [tm_lower(tp.d[1]), tm_lower(tp.d[2])]
    ok = dlo[1] > 0 && dlo[2] > 0
    fm!(tp.D, tp.d[1], tp.d[2])
    # D·K_ij = d_2 (Gadj1 : I¹_ij) + d_1 (Gadj2 : I²_ij)
    for (q, (i, j)) in enumerate(KL)
        for (t, other) in ((1, 2), (2, 1))
            tm_zero!(tp.t2)
            for a = 1:2, b = 1:2
                tm_axpy!(tp.t2, R.I[t][a, b, i, j], tp.Gadj[t][a, b], ws)
            end
            if t == 1
                fm!(tp.K[q], tp.d[2], tp.t2)
            else
                fm!(tp.t3, tp.d[1], tp.t2); tm_add!(tp.K[q], tp.K[q], tp.t3)
            end
        end
        # D·M_ij = D (d_1 M¹ + d_2 M²)
        tm_scale!(tp.t2, tp.d[1], R.M[1][i, j], ws); tm_axpy!(tp.t2, R.M[2][i, j], tp.d[2], ws)
        fm!(tp.M[q], tp.D, tp.t2)
    end
    # D·A2 = D Σ_K (d_K⟨G_K,TT^K⟩ + 2⟨G_K,TC^K⟩) + d_2⟨G_1,CC¹⟩ + d_1⟨G_2,CC²⟩
    tm_zero!(tp.t2)                                   # Σ_K (d_K⟨G,TT⟩ + 2⟨G,TC⟩)
    two = tconv(T, 2.0)
    for t = 1:2
        tm_zero!(tp.t1)
        for a = 1:2, b = 1:2; tm_axpy!(tp.t1, R.TT[t][a, b], tp.GA[t][a, b], ws); end
        fm!(tp.t3, tp.d[t], tp.t1); tm_add!(tp.t2, tp.t2, tp.t3)
        for a = 1:2, b = 1:2; tm_axpy!(tp.t2, two * R.TC[t][a, b], tp.GA[t][a, b], ws); end
    end
    fm!(tp.A2, tp.D, tp.t2)
    for (t, other) in ((1, 2), (2, 1))
        tm_zero!(tp.t1)
        for a = 1:2, b = 1:2; tm_axpy!(tp.t1, R.CC[t][a, b], tp.GA[t][a, b], ws); end
        fm!(tp.t3, tp.d[other], tp.t1); tm_add!(tp.A2, tp.A2, tp.t3)
    end
    return (; ok, dlo)
end
