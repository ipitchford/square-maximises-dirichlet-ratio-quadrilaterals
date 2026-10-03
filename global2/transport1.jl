# transport1.jl — first-order corrected transport certificate.
#
# Reference test functions depend linearly on δ = p − p0:  û_k(x̂;δ) = b_{Q0}(x̂)·(P_k⁰(x̂) + Σ_i δ_i P_kⁱ(x̂)),
# reference flux ŵ = d_K τ̂(δ) + curl ψ̂(δ), τ̂(δ) = φ̂(δ)·e with e·∇φ̂(δ) = −û_1(·;δ) (φ̂ = 0 on line v2v4 of Q0),
# ψ̂(δ) = ψ̂⁰ + Σ_i δ_i ψ̂ⁱ.  For every p, u_k(p) = û_k(·;δ)∘Ψ_p⁻¹ ∈ H¹₀(Q_p) and the Piola flux w satisfies
# div w + u_1(p) = 0 exactly (same argument as transport.jl, since div τ̂(δ) = −û_1(·;δ) for each δ).
# Base quantities are polynomials in δ:  with f^μ = b P^μ (μ = 0 base, μ = i derivative), δ^0 := 1,
#   K_ij = Σ_K (1/d_K) Σ_ab (adj adjᵀ)_ab Σ_{μν} δ^μ δ^ν I^K_ab(i^μ, j^ν),   M_ij = Σ_K d_K Σ_{μν} δ^μ δ^ν M^K(i^μ, j^ν),
#   A2 = Σ_K [ d_K⟨G,TT(δ)⟩ + 2⟨G,TC(δ)⟩ + ⟨G,CC(δ)⟩/d_K ],  TT(δ) = Σ_{μν} δ^μδ^ν TT^K(μ,ν), etc.
# The P_kⁱ, ψ̂ⁱ are central differences (step h) of the optimal choices at p0 ± h e_i, obtained WITHOUT new
# quadrature from the full-basis reference integrals transported by A_K(p) (Float64; a choice only).

# ------------------------------------------------------------------ Float64: full-basis reference data at p0
"""
    tsolve1(p0, n; kpsi = n+6, shift, h = 1e-4) -> NamedTuple

Float64 choice of the first-order corrected test polynomials and stream function.
Fields: C (N×3) base eigenvectors, dC (N×3×4) derivatives, dpsi (Np), ddpsi (Np×4), plus the data of tsolve.
"""
function tsolve1(p0::AbstractVector{<:Real}, n::Int; kpsi::Int = n + 6, shift::Int = choose_shift(p0), h::Float64 = 1e-4)
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
    # τ̂ of every basis function (line integration), components along e
    Ys = ystar_f(V); Iv = Ivec_f(R, Ys)
    PhiT = zeros(N, Q)                      # φ̂ of basis function k at node q
    for t = 1:2, i = 1:m
        base = (t - 1) * m * m + (i - 1) * m
        for k = 1:N
            pre = sum(Iv[kk] * Phi[k, base+kk] for kk = 1:m)
            for j = 1:m
                acc = sum(R.S[j, kk] * Phi[k, base+kk] for kk = 1:m)
                PhiT[k, base+j] = -(1 - R.X[i]) * (acc - pre)
            end
        end
    end
    ex = V[3, 1] - V[1, 1]; ey = V[3, 2] - V[1, 2]
    # per-triangle full-basis reference matrices
    qs = [((t - 1) * m * m + 1):(t * m * m) for t = 1:2]
    IF = [[zeros(N, N) for _ = 1:2, _ = 1:2] for _ = 1:2]; MF = [zeros(N, N) for _ = 1:2]
    TTF = [[zeros(N, N) for _ = 1:2, _ = 1:2] for _ = 1:2]; TCF = [[zeros(N, Np) for _ = 1:2, _ = 1:2] for _ = 1:2]
    CCF = [[zeros(Np, Np) for _ = 1:2, _ = 1:2] for _ = 1:2]
    for t = 1:2
        r = qs[t]; w = W[r]
        G = (Gx[:, r], Gy[:, r]); Tf = (PhiT[:, r] .* ex, PhiT[:, r] .* ey); Cf = (Cx[:, r], Cy[:, r])
        for a = 1:2, b = 1:2
            IF[t][a, b] = (G[a] .* w') * G[b]'
            TTF[t][a, b] = (Tf[a] .* w') * Tf[b]'
            TCF[t][a, b] = (Tf[a] .* w') * Cf[b]'
            CCF[t][a, b] = (Cf[a] .* w') * Cf[b]'
        end
        MF[t] = (Phi[:, r] .* w') * Phi[:, r]'
    end
    MK0 = [let (i0, i1, i2, _) = TRIS[t]
               [V[i1, 1]-V[i0, 1] V[i2, 1]-V[i0, 1]; V[i1, 2]-V[i0, 2] V[i2, 2]-V[i0, 2]]
           end for t = 1:2]
    function pencil_at(p)
        Vp = verts(p)[vperm(shift), :]
        Kp = zeros(N, N); Mp = zeros(N, N); TT = zeros(N, N); TC = zeros(N, Np); CC = zeros(Np, Np)
        for (t, (i0, i1, i2, _)) in enumerate(TRIS)
            Mk = [Vp[i1, 1]-Vp[i0, 1] Vp[i2, 1]-Vp[i0, 1]; Vp[i1, 2]-Vp[i0, 2] Vp[i2, 2]-Vp[i0, 2]]
            A = Mk / MK0[t]; d = det(A); adj = [A[2, 2] -A[1, 2]; -A[2, 1] A[1, 1]]; Ga = adj * adj'; GA = A' * A
            for a = 1:2, b = 1:2
                Kp .+= (Ga[a, b] / d) .* IF[t][a, b]
                TT .+= (d * GA[a, b]) .* TTF[t][a, b]; TC .+= GA[a, b] .* TCF[t][a, b]; CC .+= (GA[a, b] / d) .* CCF[t][a, b]
            end
            Mp .+= d .* MF[t]
        end
        return (Kp + Kp') / 2, (Mp + Mp') / 2, TT, TC, (CC + CC') / 2
    end
    function choose_at(p, Cref, Mref)
        Kp, Mp, TT, TC, CC = pencil_at(p)
        theta, U = robust_geneig(Kp, Mp)
        C = U[:, 1:3]
        if Cref !== nothing
            Gm = Cref' * Mref * C
            Gm[1, 1] < 0 && (C[:, 1] .*= -1)
            F2 = svd(Gm[2:3, 2:3]); C[:, 2:3] = C[:, 2:3] * (F2.V * F2.U')
        end
        c1 = C[:, 1]
        # ψ minimising A2(ψ) = c1ᵀTT c1 + 2 c1ᵀTC ψ + ψᵀCC ψ  (Tikhonov-regularised)
        mu = 1e-10 * opnorm(CC, 1)
        psi = -((CC + mu * I) \ (TC' * c1))
        return C, psi, theta, Mp
    end
    C0, psi0, theta0, M0 = choose_at(p0, nothing, nothing)
    dC = zeros(N, 3, 4); dpsi = zeros(Np, 4)
    for i = 1:4
        e = zeros(4); e[i] = h
        Cp, psip, _, _ = choose_at(p0 .+ e, C0, M0)
        Cm, psim, _, _ = choose_at(p0 .- e, C0, M0)
        dC[:, :, i] = (Cp .- Cm) ./ (2h); dpsi[:, i] = (psip .- psim) ./ (2h)
    end
    return (; p0, n, kpsi, shift, V, sc, theta = theta0, C = C0, dpsi = psi0, dC, ddpsi = dpsi, pencil_at)
end

# ------------------------------------------------------------------ rigorous reference integrals (matrix form)
struct TRefs1{T}
    I::Vector{Matrix{Matrix{T}}}   # I[K][a,b] :: 15×15,  rows/cols (k,μ) ↦ 5(k−1)+μ+1
    M::Vector{Matrix{T}}           # 15×15
    TT::Vector{Matrix{T}}          # 5×5   ∫ φ̂^μ φ̂^ν   (τ̂^μ_a = φ̂^μ e_a)
    TC::Vector{Matrix{T}}          # 5×10  ∫ φ̂^μ c^ν_b   (column (b,ν) ↦ 5(b−1)+ν+1)
    CC::Vector{Matrix{T}}          # 10×10 ∫ c^μ_a c^ν_b
    e::Vector{T}
    Minv0::Vector{Matrix{T}}
end
_mat(::Type{Float64}, A) = A
_mat(::Type{Arb}, A) = ArbMatrix(A; prec = PREC)
_unmat(::Type{Float64}, C) = C
_unmat(::Type{Arb}, C) = [Arb(C[i, j]; prec = PREC) for i = 1:size(C, 1), j = 1:size(C, 2)]
mmT(::Type{T}, A, B) where {T} = _unmat(T, _mat(T, A) * _mat(T, B))

# P_k(u), P_k'(u), k = 0..K  (generic, three-term recurrences; exact integer constants)
function legendre01(::Type{T}, u, K::Int) where {T}
    P = Vector{T}(undef, K + 1); D = Vector{T}(undef, K + 1)
    P[1] = tconv(T, 1.0); D[1] = tconv(T, 0.0)
    K == 0 && return P, D
    P[2] = u; D[2] = tconv(T, 1.0)
    for k = 1:K-1
        P[k+2] = ((2k + 1) * u * P[k+1] - k * P[k]) / (k + 1)
        D[k+2] = D[k] + (2k + 1) * P[k+1]
    end
    return P, D
end

"""
    transport1_refs(T, sol) -> TRefs1{T}

Reference integrals over K1(p0), K2(p0) of the 15 test-polynomial parts, the 5 flux parts and the 5 curl parts
(exact Duffy–Gauss quadrature; in Arb every entry is a ball containing the exact integral).
"""
function transport1_refs(::Type{T}, sol) where {T}
    p0 = sol.p0; n = sol.n; kpsi = sol.kpsi; sc = sol.sc; Kd = max(n, kpsi)
    m = npoints(n); NCs = nodeconsts(T, m); Qk = m * m
    V = verts_T(T, p0)[vperm(sol.shift), :]
    isx = tconv(T, 1.0) / tconv(T, sc.sx); isy = tconv(T, 1.0) / tconv(T, sc.sy)
    txc = tconv(T, sc.xc); tyc = tconv(T, sc.yc)
    E = [V[mod1(k + 1, 4), c] - V[k, c] for k = 1:4, c = 1:2]
    Al = [tconv(T, 0.0) for _ = 1:4, _ = 1:4]
    for k = 1:4, j = 1:4
        (j == k || j == mod1(k + 1, 4)) && continue
        Al[k, j] = E[k, 1] * (V[j, 2] - V[k, 2]) - E[k, 2] * (V[j, 1] - V[k, 1])
    end
    gx = [-E[k, 2] for k = 1:4]; gy = [E[k, 1] for k = 1:4]
    ex = V[3, 1] - V[1, 1]; ey = V[3, 2] - V[1, 2]
    fx_ = V[4, 1] - V[2, 1]; fy_ = V[4, 2] - V[2, 2]
    Ys = (fx_ * (V[2, 2] - V[1, 2]) - fy_ * (V[2, 1] - V[1, 1])) / (fx_ * ey - fy_ * ex)
    Ra = rule_arb(m)
    Ivec = Vector{T}(undef, m)
    for k = 1:m
        acc = tconv(T, 0.0)
        for g = 1:m
            y = Ys * tconv(T, Ra.X[g]); L = tconv(T, 1.0)
            for i = 1:m; i != k && (L *= (y - tconv(T, Ra.X[i])) / (tconv(T, Ra.X[k]) - tconv(T, Ra.X[i]))); end
            acc += tconv(T, Ra.W[g]) * L
        end
        Ivec[k] = Ys * acc
    end
    idx = basis_index(n); N = length(idx); pidx = psi_index(kpsi); Np = length(pidx)
    # coefficient matrices: 15 × N (rows (k,μ)) and 5 × Np
    Cmat = Matrix{T}(undef, 15, N)
    for k = 1:3, μ = 0:4, c = 1:N
        Cmat[5(k-1)+μ+1, c] = tconv(T, μ == 0 ? sol.C[c, k] : sol.dC[c, k, μ])
    end
    Pmat = Matrix{T}(undef, 5, Np)
    for μ = 0:4, c = 1:Np
        Pmat[μ+1, c] = tconv(T, μ == 0 ? sol.dpsi[c] : sol.ddpsi[c, μ])
    end
    Iout = Matrix{Matrix{T}}[]; Mout = Matrix{T}[]; TTo = Matrix{T}[]; TCo = Matrix{T}[]; CCo = Matrix{T}[]; Minv0 = Matrix{T}[]
    for (t, (i0, i1, i2, sg)) in enumerate(TRIS)
        Mk = [V[i1, 1]-V[i0, 1] V[i2, 1]-V[i0, 1]; V[i1, 2]-V[i0, 2] V[i2, 2]-V[i0, 2]]
        detM = Mk[1, 1] * Mk[2, 2] - Mk[1, 2] * Mk[2, 1]
        push!(Minv0, [Mk[2, 2]/detM -Mk[1, 2]/detM; -Mk[2, 1]/detM Mk[1, 1]/detM])
        Jabs = tconv(T, sg) * detM
        B0 = Matrix{T}(undef, N, Qk); Bu = similar(B0); Bv = similar(B0)
        Su = Matrix{T}(undef, Np, Qk); Sv = similar(Su)
        bq = Vector{T}(undef, Qk); bxq = similar(bq); byq = similar(bq); wq = similar(bq)
        for i = 1:m, j = 1:m
            q = (i - 1) * m + j
            lam0 = NCs.lam0[i, j]; xi = NCs.xi[i]; eta = NCs.eta[i, j]
            x = lam0 * V[i0, 1] + xi * V[i1, 1] + eta * V[i2, 1]
            y = lam0 * V[i0, 2] + xi * V[i1, 2] + eta * V[i2, 2]
            PU, DU = legendre01(T, (x - txc) * isx, Kd); PV, DV_ = legendre01(T, (y - tyc) * isy, Kd)
            for (c, (a, bb)) in enumerate(idx)
                B0[c, q] = PU[a+1] * PV[bb+1]; Bu[c, q] = (DU[a+1] * PV[bb+1]) * isx; Bv[c, q] = (PU[a+1] * DV_[bb+1]) * isy
            end
            for (c, (a, bb)) in enumerate(pidx)
                Su[c, q] = (DU[a+1] * PV[bb+1]) * isx; Sv[c, q] = (PU[a+1] * DV_[bb+1]) * isy
            end
            l = [lam0 * Al[k, i0] + xi * Al[k, i1] + eta * Al[k, i2] for k = 1:4]
            L12 = l[1] * l[2]; L34 = l[3] * l[4]
            bq[q] = L12 * L34
            tt = (L34 * l[2], L34 * l[1], L12 * l[4], L12 * l[3])
            bxq[q] = gx[1] * tt[1] + gx[2] * tt[2] + gx[3] * tt[3] + gx[4] * tt[4]
            byq[q] = gy[1] * tt[1] + gy[2] * tt[2] + gy[3] * tt[3] + gy[4] * tt[4]
            wq[q] = NCs.w[i, j] * Jabs
        end
        Pv = mmT(T, Cmat, B0); Pu = mmT(T, Cmat, Bu); Pw = mmT(T, Cmat, Bv)        # 15 × Q
        F = similar(Pv); Fx = similar(Pv); Fy = similar(Pv)
        for q = 1:Qk, r = 1:15
            F[r, q] = bq[q] * Pv[r, q]
            Fx[r, q] = bxq[q] * Pv[r, q] + bq[q] * Pu[r, q]
            Fy[r, q] = byq[q] * Pv[r, q] + bq[q] * Pw[r, q]
        end
        # flux parts: φ̂^μ of rows (1, μ) by line integration on each Duffy line
        Phi = Matrix{T}(undef, 5, Qk)
        for μ = 0:4, i = 1:m
            r = μ + 1; base = (i - 1) * m
            pre = tconv(T, 0.0)
            for k = 1:m; pre += Ivec[k] * F[r, base+k]; end
            for j = 1:m
                acc = tconv(T, 0.0)
                for k = 1:m; acc += NCs.S[j, k] * F[r, base+k]; end
                Phi[r, base+j] = -(NCs.omx[i] * (acc - pre))
            end
        end
        Cv = mmT(T, Pmat, Sv); Cu = mmT(T, Pmat, Su)                                 # 5 × Q: ψ_y, ψ_x
        Cb = Matrix{T}(undef, 10, Qk)                                                  # rows (b,ν): c_x = ψ_y, c_y = −ψ_x
        for q = 1:Qk, ν = 1:5
            Cb[ν, q] = Cv[ν, q]; Cb[5+ν, q] = -Cu[ν, q]
        end
        # weighted copies and products
        Fw = similar(F); Fxw = similar(F); Fyw = similar(F); Phw = similar(Phi); Cbw = similar(Cb)
        for q = 1:Qk
            for r = 1:15; Fw[r, q] = wq[q] * F[r, q]; Fxw[r, q] = wq[q] * Fx[r, q]; Fyw[r, q] = wq[q] * Fy[r, q]; end
            for r = 1:5; Phw[r, q] = wq[q] * Phi[r, q]; end
            for r = 1:10; Cbw[r, q] = wq[q] * Cb[r, q]; end
        end
        G = (Fx, Fy); Gw = (Fxw, Fyw)
        push!(Iout, [mmT(T, Gw[a], permutedims(G[b])) for a = 1:2, b = 1:2])
        push!(Mout, mmT(T, Fw, permutedims(F)))
        push!(TTo, mmT(T, Phw, permutedims(Phi)))
        push!(TCo, mmT(T, Phw, permutedims(Cb)))
        push!(CCo, mmT(T, Cbw, permutedims(Cb)))
    end
    return TRefs1{T}(Iout, Mout, TTo, TCo, CCo, [ex, ey], Minv0)
end

# exact quadratic TM  Σ_{μ,ν} δ^μ δ^ν X[r(μ), c(ν)],  δ^0 = 1, δ_i = rad_i t_i
function quadform_tm!(h::TM{T}, X::Matrix{T}, rows, cols, rad::Vector{Float64}, scale) where {T}
    tab = TAB[]
    tm_zero!(h)
    rT = [tconv(T, x) for x in rad]
    _acc!(h, 1, scale * X[rows[1], cols[1]])
    for l = 1:4
        _acc!(h, tab.lin[l], scale * rT[l] * (X[rows[1], cols[l+1]] + X[rows[l+1], cols[1]]))
    end
    for k = 1:4, l = k:4
        v = k == l ? X[rows[k+1], cols[k+1]] : X[rows[k+1], cols[l+1]] + X[rows[l+1], cols[k+1]]
        _acc!(h, tab.quad[k, l], scale * rT[k] * rT[l] * v)
    end
    h
end
_acc!(h::TM{Float64}, k, x) = (h.c[k] += x)
_acc!(h::TM{Arb}, k, x) = Arblib.add!(h.c[k], h.c[k], x)

"""
    transport1_base!(tp, R, p0, rad, ws, shift) -> (; ok, dlo)

Taylor models over the box of D·K, D·M, D·A2 for the first-order corrected transport certificate.
"""
function transport1_base!(tp::TPool{T}, R::TRefs1{T}, p0::Vector{Float64}, rad::Vector{Float64}, ws, shift::Int) where {T}
    G0 = box_geom(T, p0, rad, ws)
    Gv = G0.V[vperm(shift), :]
    fm!(h, f, g) = tm_mul_full!(h, f, g, ws, tp.buf)
    for (t, (i0, i1, i2, sg)) in enumerate(TRIS)
        Mp = [newtm(T) for _ = 1:2, _ = 1:2]
        for c = 1:2
            tm_sub!(Mp[c, 1], Gv[i1, c], Gv[i0, c]); tm_sub!(Mp[c, 2], Gv[i2, c], Gv[i0, c])
        end
        A = tp.A[t]; Mi = R.Minv0[t]
        for a = 1:2, b = 1:2
            tm_zero!(A[a, b]); tm_axpy!(A[a, b], Mi[1, b], Mp[a, 1], ws); tm_axpy!(A[a, b], Mi[2, b], Mp[a, 2], ws)
        end
        fm!(tp.d[t], A[1, 1], A[2, 2]); fm!(tp.t1, A[1, 2], A[2, 1]); tm_sub!(tp.d[t], tp.d[t], tp.t1)
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
    one_ = tconv(T, 1.0)
    rowsk(k) = [5(k - 1) + μ + 1 for μ = 0:4]
    q5 = newtm(T); s4 = newtm(T)
    for (q, (i, j)) in enumerate(KL)
        for (t, other) in ((1, 2), (2, 1))
            tm_zero!(tp.t2)                                # Σ_ab Gadj_ab ⊗ Iδ_ab
            for a = 1:2, b = 1:2
                quadform_tm!(q5, R.I[t][a, b], rowsk(i), rowsk(j), rad, one_)
                fm!(s4, tp.Gadj[t][a, b], q5); tm_add!(tp.t2, tp.t2, s4)
            end
            if t == 1
                fm!(tp.K[q], tp.d[2], tp.t2)
            else
                fm!(tp.t3, tp.d[1], tp.t2); tm_add!(tp.K[q], tp.K[q], tp.t3)
            end
        end
        tm_zero!(tp.t2)                                    # d_1 Mδ¹ + d_2 Mδ²
        for t = 1:2
            quadform_tm!(q5, R.M[t], rowsk(i), rowsk(j), rad, one_)
            fm!(s4, tp.d[t], q5); tm_add!(tp.t2, tp.t2, s4)
        end
        fm!(tp.M[q], tp.D, tp.t2)
    end
    # A2:  D Σ_K (d_K⟨G,TTδ⟩ + 2⟨G,TCδ⟩) + d_2⟨G_1,CCδ¹⟩ + d_1⟨G_2,CCδ²⟩
    r5 = collect(1:5)
    tm_zero!(tp.t2)
    for t = 1:2
        tm_zero!(tp.t1)
        quadform_tm!(q5, R.TT[t], r5, r5, rad, one_)              # Σ δδ ∫φ̂φ̂ ; TT_ab = e_a e_b ·
        egg = newtm(T); tm_zero!(egg)
        for a = 1:2, b = 1:2; tm_axpy!(egg, R.e[a] * R.e[b], tp.GA[t][a, b], ws); end   # ⟨G, e eᵀ⟩
        fm!(s4, egg, q5); fm!(tp.t3, tp.d[t], s4); tm_add!(tp.t2, tp.t2, tp.t3)
        for a = 1:2, b = 1:2                                       # 2⟨G, TCδ⟩, TC_ab = e_a Σ δδ ∫φ̂ c_b
            quadform_tm!(q5, R.TC[t], r5, [5(b - 1) + ν for ν = 1:5], rad, tconv(T, 2.0) * R.e[a])
            fm!(s4, tp.GA[t][a, b], q5); tm_add!(tp.t2, tp.t2, s4)
        end
    end
    fm!(tp.A2, tp.D, tp.t2)
    for (t, other) in ((1, 2), (2, 1))
        tm_zero!(tp.t1)
        for a = 1:2, b = 1:2
            quadform_tm!(q5, R.CC[t], [5(a - 1) + μ for μ = 1:5], [5(b - 1) + ν for ν = 1:5], rad, one_)
            fm!(s4, tp.GA[t][a, b], q5); tm_add!(tp.t1, tp.t1, s4)
        end
        fm!(tp.t3, tp.d[other], tp.t1); tm_add!(tp.A2, tp.A2, tp.t3)
    end
    return (; ok, dlo)
end
