# assemble.jl — Taylor models (order 2 in t = δ/r ∈ [−1,1]^4, δ = p − p0) of the base quantities
#   M_kl = (u_k, u_l),  K_kl = (∇u_k, ∇u_l)  (k ≤ l ∈ {1,2,3}),  A2 = ‖w‖²,
# u_k = b(x;p)·P_k(x) with FIXED polynomials P_k (coefficients from the Float64 solve at p0) and the bubble b,
# w = ∇u1/ρ0 + φ e + curl ψ (e = v3 − v1, φ = 0 on the line v2v4, ψ fixed), by EXACT quadrature:
# two triangles (v1,v2,v3), (v1,v4,v3), Duffy-collapsed m×m Gauss–Legendre product rule, m = n+6
# (exact for total degree ≤ 2m−2 = 2n+10 ≥ deg |w|² = 2(n+5); see README_global2.md).
# Every node is an affine function of p (fixed barycentric weights), so the node displacement is an exact
# linear TM; the fixed polynomials are evaluated by their FINITE Taylor expansion at the node centre with the
# derivative table computed at the (thin) centre (no dependency blow-up), remainder = Σ_{a+b≥3}|D_ab|·H^a/a!·H^b/b!.

# ------------------------------------------------------------------ kernels: Float64
function legendre_table!(PD::Matrix{Float64}, u::Float64, K::Int)
    fill!(PD, 0.0)
    PD[1, 1] = 1.0
    K == 0 && return PD
    PD[2, 1] = u; PD[2, 2] = 1.0
    for k = 1:K-1
        a = 2k + 1
        PD[k+2, 1] = (a * u * PD[k+1, 1] - k * PD[k, 1]) / (k + 1)
        for d = 1:k+1
            PD[k+2, d+1] = PD[k, d+1] + a * PD[k+1, d]
        end
    end
    PD
end
# D[a+1,b+1] = ∂_u^a ∂_v^b P(u0,v0), P = Σ_{i+l≤d} C[i+1,l+1] P_i(u) P_l(v)
function dtable!(D::Matrix{Float64}, Cb::Matrix{Float64}, C::Matrix{Float64}, PU::Matrix{Float64}, PV::Matrix{Float64}, d::Int)
    fill!(D, 0.0)
    for i = 0:d, b = 0:d-i
        s = 0.0
        for l = b:d-i; s += C[i+1, l+1] * PV[l+1, b+1]; end
        Cb[i+1, b+1] = s
    end
    for a = 0:d, b = 0:d-a
        s = 0.0
        for i = a:d-b; s += PU[i+1, a+1] * Cb[i+1, b+1]; end
        D[a+1, b+1] = s
    end
    D
end
function absdtable!(AD::Matrix{Float64}, D::Matrix{Float64}, d::Int)
    for a = 0:d, b = 0:d-a; AD[a+1, b+1] = abs(D[a+1, b+1]); end
    AD
end

# ------------------------------------------------------------------ kernels: Arb (in place)
function legendre_table!(PD::Matrix{Arb}, u::Arb, K::Int)
    for x in PD; Arblib.zero!(x); end
    Arblib.one!(PD[1, 1])
    K == 0 && return PD
    Arblib.set!(PD[2, 1], u); Arblib.one!(PD[2, 2])
    for k = 1:K-1
        a = 2k + 1
        z = PD[k+2, 1]
        Arblib.mul!(z, u, PD[k+1, 1]); Arblib.mul!(z, z, a)
        Arblib.submul!(z, PD[k, 1], k)
        Arblib.div!(z, z, k + 1)
        for d = 1:k+1
            z = PD[k+2, d+1]
            Arblib.set!(z, PD[k, d+1]); Arblib.addmul!(z, PD[k+1, d], a)
        end
    end
    PD
end
function dtable!(D::Matrix{Arb}, Cb::Matrix{Arb}, C::Matrix{Arb}, PU::Matrix{Arb}, PV::Matrix{Arb}, d::Int)
    for x in D; Arblib.zero!(x); end
    for i = 0:d, b = 0:d-i
        z = Cb[i+1, b+1]; Arblib.zero!(z)
        for l = b:d-i; Arblib.addmul!(z, C[i+1, l+1], PV[l+1, b+1]); end
    end
    for a = 0:d, b = 0:d-a
        z = D[a+1, b+1]
        for i = a:d-b; Arblib.addmul!(z, PU[i+1, a+1], Cb[i+1, b+1]); end
    end
    D
end
function absdtable!(AD::Matrix{Mag}, D::Matrix{Arb}, d::Int)
    for a = 0:d, b = 0:d-a; Arblib.get!(AD[a+1, b+1], D[a+1, b+1]); end
    AD
end

# ------------------------------------------------------------------ Taylor TM of a fixed polynomial at a moving node
# node displacement in scaled coordinates: hu(t) = Σ hu_i t_i, hv(t) = Σ hv_i t_i (exact, linear).
# HUV[a+1,b+1] = hu^a hv^b /(a! b!)  (exact TMs of degree a+b ≤ ORDER, no remainder)
# P^{(su,sv)}(x0 + h) = Σ_{a+b ≤ ORDER} D[a+su,b+sv] HUV[a,b]  +  [−ρ, ρ],
#   ρ = Σ_{a+b > ORDER} |D[a+su,b+sv]| HU[a] HV[b],  HU[a+1] ≥ (Σ|hu_i|)^a/a!.
function build_huv!(HUV::Matrix{TM{T,R}}, hu::Vector{T}, hv::Vector{T}, ws, invfact::Vector{T}) where {T,R}
    K = TAB[].order; lin = TAB[].lin
    tm_zero!(HUV[1, 1]); _setc!(HUV[1, 1], 1, tconv(T, 1.0))
    if K >= 1
        tm_zero!(HUV[2, 1]); tm_zero!(HUV[1, 2])
        for i = 1:4; _setc!(HUV[2, 1], lin[i], hu[i]); _setc!(HUV[1, 2], lin[i], hv[i]); end
    end
    for d = 2:K, a = 0:d
        b = d - a
        if a >= 1
            tm_mul!(HUV[a+1, b+1], HUV[a, b+1], HUV[2, 1], ws)       # (hu^{a-1} hv^b/((a-1)! b!)) · hu  then /a
            tm_scale!(HUV[a+1, b+1], HUV[a+1, b+1], invfact[a+1], ws)
        else
            tm_mul!(HUV[a+1, b+1], HUV[a+1, b], HUV[1, 2], ws)
            tm_scale!(HUV[a+1, b+1], HUV[a+1, b+1], invfact[b+1], ws)
        end
    end
    HUV
end
# invfact[k+1] = 1/k  (so that repeated multiplication builds 1/a!)
function taylor_tm!(h::TM{Float64}, D::Matrix{Float64}, AD::Matrix{Float64}, su::Int, sv::Int, d::Int,
                    HUV::Matrix{TM{Float64,Float64}}, HU, HV, K::Int = TAB[].order)
    tab = TAB[]
    fill!(h.c, 0.0)
    for a = 0:K, b = 0:K-a
        Dab = D[a+su+1, b+sv+1]
        Dab == 0.0 && continue
        hc = HUV[a+1, b+1].c
        for k in tab.degrange[a+b+1]; h.c[k] += Dab * hc[k]; end
    end
    rem = 0.0
    dd = d - su - sv
    for a = 0:dd, b = max(0, K + 1 - a):dd-a
        rem += AD[a+su+1, b+sv+1] * HU[a+1] * HV[b+1]
    end
    h.r[1] = rem
    h
end
function taylor_tm!(h::TM{Arb}, D::Matrix{Arb}, AD::Matrix{Mag}, su::Int, sv::Int, d::Int,
                    HUV::Matrix{TM{Arb,Mag}}, HU, HV, tmpm::Mag, K::Int = TAB[].order)
    tab = TAB[]
    for x in h.c; Arblib.zero!(x); end
    for a = 0:K, b = 0:K-a
        Dab = D[a+su+1, b+sv+1]
        hc = HUV[a+1, b+1].c
        for k in tab.degrange[a+b+1]; Arblib.addmul!(h.c[k], Dab, hc[k]); end
    end
    r = h.r[1]; Arblib.zero!(r)
    dd = d - su - sv
    for a = 0:dd, b = max(0, K + 1 - a):dd-a
        Arblib.mul!(tmpm, AD[a+su+1, b+sv+1], HU[a+1]); Arblib.mul!(tmpm, tmpm, HV[b+1]); Arblib.add!(r, r, tmpm)
    end
    h
end

# Family TM:  P(x0+h; p) = P0(x0+h) + Σ_i δ_i P_i(x0+h),  δ_i = rad_i t_i.
# Base part: Taylor TM of order ORDER; derivative parts: Taylor TM of order ORDER−1 times t_i (exact degree shift,
# remainder unchanged since |t_i| ≤ 1) times rad_i.
function family_tm!(h::TM{Float64}, Ds, ADs, su, sv, d, HUV, HU, HV, rad::Vector{Float64}, nder::Int, tmp::TM{Float64})
    K = TAB[].order; sh = TAB[].shift
    taylor_tm!(h, Ds[1], ADs[1], su, sv, d, HUV, HU, HV, K)
    for i = 1:nder
        taylor_tm!(tmp, Ds[1+i], ADs[1+i], su, sv, d, HUV, HU, HV, K - 1)
        ri = rad[i]
        for k in eachindex(tmp.c)
            j = sh[k, i]; j == 0 && continue
            h.c[j] += ri * tmp.c[k]
        end
        h.r[1] += ri * tmp.r[1]
    end
    h
end
function family_tm!(h::TM{Arb}, Ds, ADs, su, sv, d, HUV, HU, HV, rad::Vector{Arb}, nder::Int, tmp::TM{Arb}, tmpm::Mag, radm::Vector{Mag})
    K = TAB[].order; sh = TAB[].shift
    taylor_tm!(h, Ds[1], ADs[1], su, sv, d, HUV, HU, HV, tmpm, K)
    for i = 1:nder
        taylor_tm!(tmp, Ds[1+i], ADs[1+i], su, sv, d, HUV, HU, HV, tmpm, K - 1)
        ri = rad[i]
        for k in eachindex(tmp.c)
            j = sh[k, i]; j == 0 && continue
            Arblib.addmul!(h.c[j], ri, tmp.c[k])
        end
        Arblib.mul!(tmpm, radm[i], tmp.r[1]); Arblib.add!(h.r[1], h.r[1], tmpm)
    end
    h
end

# ------------------------------------------------------------------ scalar conversion
tconv(::Type{Float64}, x::Real) = Float64(x)
tconv(::Type{Arb}, x::Real) = Arb(x; prec = PREC)
tconv(::Type{Float64}, x::Arb) = Float64(Arblib.midref(x))
tconv(::Type{Arb}, x::Arb) = x
rconv(::Type{Float64}) = 0.0
rconv(::Type{Arb}) = Mag(0)
_abs_into!(m::Mag, x::Arb) = (Arblib.get!(m, x); m)

# polynomial coefficient matrix from a coefficient vector in the given index list
function coeffmat(::Type{T}, c::AbstractVector, idx, d::Int) where {T}
    C = [tconv(T, 0.0) for _ = 1:d+1, _ = 1:d+1]
    for (k, (i, l)) in enumerate(idx)
        C[i+1, l+1] = tconv(T, Float64(c[k]))
    end
    C
end

# ------------------------------------------------------------------ rule-derived node constants (cached)
struct NodeConsts{T}
    m::Int
    lam0::Matrix{T}; xi::Vector{T}; eta::Matrix{T}; w::Matrix{T}; omx::Vector{T}
    S::Matrix{T}
end
const _NC_CACHE = Dict{Tuple{DataType,Int},Any}()
const _NC_LOCK = ReentrantLock()
function nodeconsts(::Type{T}, m::Int) where {T}
    lock(_NC_LOCK) do
        get!(_NC_CACHE, (T, m)) do
            R = rule_arb(m)
            lam0 = [tconv(T, (1 - R.X[i]) * (1 - R.X[j])) for i = 1:m, j = 1:m]
            xi = [tconv(T, R.X[i]) for i = 1:m]
            eta = [tconv(T, R.X[j] * (1 - R.X[i])) for i = 1:m, j = 1:m]
            w = [tconv(T, R.W[i] * R.W[j] * (1 - R.X[i])) for i = 1:m, j = 1:m]
            omx = [tconv(T, 1 - R.X[i]) for i = 1:m]
            S = [tconv(T, R.S[j, k]) for j = 1:m, k = 1:m]
            NodeConsts{T}(m, lam0, xi, eta, w, omx, S)
        end
    end
end

# ------------------------------------------------------------------ the workspace
mutable struct AsmWork{T,R}
    ws::TMWork{T,R}
    pool::Dict{Symbol,Any}
end
function tmpool(::Type{T}) where {T}
    P = Dict{Symbol,Any}()
    for s in (:b, :bx, :by, :lapb, :L12, :L34, :L13, :L24, :L14, :L23, :t1, :t2, :t3, :t4, :tmp, :tmp2, :tmp3,
              :p1, :p1u, :p1v, :p1uu, :p1vv, :p2, :p2u, :p2v, :p3, :p3u, :p3v, :psu, :psv,
              :u1, :u2, :u3, :ux1, :ux2, :ux3, :uy1, :uy2, :uy3, :lap1, :W, :Wu, :Wux, :Wuy,
              :e_x, :e_y, :ystar, :num, :den, :inv, :pre, :phi, :wx, :wy, :Ww, :hlin, :hsq)
        P[s] = newtm(T)
    end
    P[:V] = [newtm(T) for _ = 1:4, _ = 1:2]
    P[:E] = [newtm(T) for _ = 1:4, _ = 1:2]          # edge vectors e_k = v_{k+1} − v_k
    P[:A] = [newtm(T) for _ = 1:4, _ = 1:4]          # A[k,j] = l_k(v_j)
    P[:G] = [newtm(T) for _ = 1:6]                    # g12, g34, g13, g24, g14, g23
    P[:ell] = [newtm(T) for _ = 1:4]
    P[:M] = [newtm(T) for _ = 1:6]; P[:K] = [newtm(T) for _ = 1:6]; P[:A2] = newtm(T)
    TMT = tmtype(T)
    P[:I] = TMT[]                                      # I_k(Y*), filled per m
    P[:rline] = TMT[]; P[:gxl] = TMT[]; P[:gyl] = TMT[]; P[:Wl] = TMT[]
    P
end
newasm(::Type{T}) where {T} = AsmWork{T,T == Arb ? Mag : Float64}(newwork(T), tmpool(T))

function ensure_lines!(aw::AsmWork{T}, m::Int) where {T}
    P = aw.pool
    for s in (:I, :rline, :gxl, :gyl, :Wl)
        while length(P[s]) < m; push!(P[s], newtm(T)); end
    end
end

# product index pairs for M/K storage: 1:(1,1) 2:(1,2) 3:(1,3) 4:(2,2) 5:(2,3) 6:(3,3)
const KL = ((1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3))

# ------------------------------------------------------------------ generic small helpers on TMs
# linear TM from constant + 4 coefficients
function tm_linear!(h::TM{T}, c0, lin) where {T}
    tm_zero!(h)
    _setc!(h, 1, c0)
    L = TAB[].lin
    for i = 1:4; _setc!(h, L[i], lin[i]); end
    h
end
_setc!(h::TM{Float64}, k, x) = (h.c[k] = Float64(x))
_setc!(h::TM{Arb}, k, x) = Arblib.set!(h.c[k], x)

# TM reciprocal.  F = a point value (midpoint of the constant coefficient), u = (f − F)/F (a TM, |u| ≤ U < 1),
# 1/f = (1/F)·(Σ_{j=0}^{K} (−u)^j + ε),  |ε| = |u|^{K+1}/|1+u| ≤ U^{K+1}/(1−U).
function tm_inv!(h::TM{Float64}, f::TM{Float64}, aw)
    P = aw.pool; ws = aw.ws; u = P[:hlin]; up = P[:hsq]; tmp = P[:tmp3]
    K = TAB[].order
    F = f.c[1]
    tm_set!(u, f); u.c[1] -= F; tm_scale!(u, u, 1 / F)
    U = polybound(u) + u.r[1]
    U < 1 || return false
    tm_const!(h, 1.0); tm_const!(up, 1.0)
    for j = 1:K
        tm_mul!(tmp, up, u, ws); tm_set!(up, tmp)
        tm_axpy!(h, (-1.0)^j, up)
    end
    h.r[1] += U^(K + 1) / (1 - U)
    tm_scale!(h, h, 1 / F)
    true
end
function tm_inv!(h::TM{Arb}, f::TM{Arb}, aw)
    P = aw.pool; ws = aw.ws; u = P[:hlin]::TM{Arb,Mag}; up = P[:hsq]::TM{Arb,Mag}; tmp = P[:tmp3]::TM{Arb,Mag}
    K = TAB[].order
    F = Arb(Arblib.midref(f.c[1]); prec = PREC)          # exact point
    Arblib.iszero(F) && return false
    iF = 1 / F
    tm_set!(u, f); Arblib.sub!(u.c[1], u.c[1], F); tm_scale!(u, u, iF, ws)
    Ubm = Mag(0); polybound!(Ubm, u, ws); Arblib.add!(Ubm, Ubm, u.r[1])
    Ua = Arb(Ubm; prec = PREC)
    Arblib.ubound(Ua) < 1 || return false
    tm_const!(h, Arb(1; prec = PREC)); tm_const!(up, Arb(1; prec = PREC))
    for j = 1:K
        tm_mul!(tmp, up, u, ws); tm_set!(up, tmp)
        tm_axpy!(h, Arb(isodd(j) ? -1 : 1; prec = PREC), up, ws)
    end
    rem = Ua^(K + 1) / (1 - Ua)
    rm = Mag(0); Arblib.get!(rm, rem); Arblib.add!(h.r[1], h.r[1], rm)
    tm_scale!(h, h, iF, ws)
    true
end

# polynomial q(Y) = Σ_j a_j (Y − Y0)^j at a TM Y = Y0 + hY: terms j ≤ ORDER by TM powers (hpow[j] = hY^j,
# computed by TM products, each with its own remainder), terms j > ORDER into the remainder via |hY| ≤ H.
function tm_poly_at!(h::TM{Float64}, a::Vector{Float64}, hpow::Vector{TM{Float64,Float64}}, H::Float64)
    K = TAB[].order
    tm_const!(h, a[1])
    for j = 1:min(K, length(a) - 1); tm_axpy!(h, a[j+1], hpow[j]); end
    rem = 0.0; Hp = H^(K + 1)
    for j = K+1:length(a)-1; rem += abs(a[j+1]) * Hp; Hp *= H; end
    h.r[1] += rem
    h
end
function tm_poly_at!(h::TM{Arb}, a::Vector{Arb}, hpow::Vector{TM{Arb,Mag}}, H::Mag, ws)
    K = TAB[].order
    tm_const!(h, a[1])
    for j = 1:min(K, length(a) - 1); tm_axpy!(h, a[j+1], hpow[j], ws); end
    Hp = Mag(1); for _ = 1:K+1; Arblib.mul!(Hp, Hp, H); end
    t = Mag(0)
    for j = K+1:length(a)-1
        Arblib.get!(t, a[j+1]); Arblib.mul!(t, t, Hp); Arblib.add!(h.r[1], h.r[1], t); Arblib.mul!(Hp, Hp, H)
    end
    h
end

# Taylor coefficients (in powers of (Y − Y0)) of I_k(Y) = ∫_0^Y L_k, k = 1..m, at Y0, in type T.
# I_k(Y0 + s) = I_k(Y0) + Σ_{j≥0} l_j s^{j+1}/(j+1),  L_k(Y0 + s) = Σ_j l_j s^j (degree m−1).
function Ik_taylor(::Type{T}, m::Int, Y0::Float64) where {T}
    R = rule_arb(m)
    out = Vector{Vector{T}}(undef, m)
    Y0a = Arb(Y0; prec = PREC)
    for k = 1:m
        # power series of L_k(Y0 + s) = Π_{i≠k} (Y0 − X_i + s)/(X_k − X_i)
        ser = [Arb(1; prec = PREC)]
        for i = 1:m
            i == k && continue
            den = R.X[k] - R.X[i]; c0 = (Y0a - R.X[i]) / den; c1 = 1 / den
            new = [Arb(0; prec = PREC) for _ = 1:length(ser)+1]
            for j = 1:length(ser)
                new[j] += ser[j] * c0; new[j+1] += ser[j] * c1
            end
            ser = new
        end
        # I_k(Y0) by the m-point Gauss rule on [0, Y0]
        acc = Arb(0; prec = PREC)
        for g = 1:m
            y = Y0a * R.X[g]; L = Arb(1; prec = PREC)
            for i = 1:m; i != k && (L *= (y - R.X[i]) / (R.X[k] - R.X[i])); end
            acc += R.W[g] * L
        end
        coeffs = [Y0a * acc]
        for j = 1:length(ser); push!(coeffs, ser[j] / j); end
        out[k] = [tconv(T, x) for x in coeffs]
    end
    out
end

# ------------------------------------------------------------------ main assembly
"""
    assemble_box!(aw, p0, rad, sol) -> (; M, K, A2, ok)

TMs over the box p0 + Π[−rad_i, rad_i] of the base quantities, for the test functions/flux in `sol`
(fields: n, kpsi, sc::Scaling, C (N×3), dpsi, rho0).  Mode is the element type of `aw`.
M[q], K[q] for q ↔ (k,l) in KL = ((1,1),(1,2),(1,3),(2,2),(2,3),(3,3)).
"""
function assemble_box!(aw::AsmWork{T,R}, p0::Vector{Float64}, rad::Vector{Float64}, sol) where {T,R}
    TMT = TM{T,R}
    P = aw.pool; ws = aw.ws
    n = sol.n; kpsi = sol.kpsi; sc = sol.sc; Kd = max(n, kpsi)
    m = npoints(n); ensure_lines!(aw, m)
    NCst = nodeconsts(T, m)::NodeConsts{T}
    isArb = T == Arb
    V = P[:V]::Matrix{TMT}; E = P[:E]::Matrix{TMT}; A = P[:A]::Matrix{TMT}; G = P[:G]::Vector{TMT}; ell = P[:ell]::Vector{TMT}
    tmp = P[:tmp]::TMT; tmp2 = P[:tmp2]::TMT; tmp3 = P[:tmp3]::TMT
    # ---- vertices (exact linear TMs)
    for j = 1:4, k = 1:2
        if isArb
            c0 = Arb(V0[j, k]; prec = PREC)
            for i = 1:4; Arblib.addmul!(c0, Arb(DV[j, k, i]; prec = PREC), Arb(p0[i]; prec = PREC)); end
            tm_linear!(V[j, k], c0, [Arb(DV[j, k, i]; prec = PREC) * Arb(rad[i]; prec = PREC) for i = 1:4])
        else
            tm_linear!(V[j, k], V0[j, k] + sum(DV[j, k, i] * p0[i] for i = 1:4), [DV[j, k, i] * rad[i] for i = 1:4])
        end
    end
    for k = 1:4, c = 1:2; tm_sub!(E[k, c], V[mod1(k + 1, 4), c], V[k, c]); end
    # A[k,j] = l_k(v_j) = cross(e_k, v_j − v_k) = e_kx (v_jy − v_ky) − e_ky (v_jx − v_kx)   (exact quadratics)
    for k = 1:4, j = 1:4
        tm_zero!(A[k, j])
        (j == k || j == mod1(k + 1, 4)) && continue
        tm_sub!(tmp, V[j, 1], V[k, 1]); tm_sub!(tmp2, V[j, 2], V[k, 2])
        tm_mul!(A[k, j], E[k, 1], tmp2, ws)
        tm_mul!(tmp3, E[k, 2], tmp, ws); tm_sub!(A[k, j], A[k, j], tmp3)
    end
    # gradients of l_k: g_k = (−e_ky, e_kx);  g_i·g_j = e_i·e_j
    pairs6 = ((1, 2), (3, 4), (1, 3), (2, 4), (1, 4), (2, 3))
    for (q, (i, j)) in enumerate(pairs6)
        tm_mul!(G[q], E[i, 1], E[j, 1], ws); tm_addmul!(G[q], E[i, 2], E[j, 2], ws)
    end
    ex = P[:e_x]::TMT; ey = P[:e_y]::TMT
    tm_sub!(ex, V[3, 1], V[1, 1]); tm_sub!(ey, V[3, 2], V[1, 2])
    # ---- Y* = cross(v4−v2, v2−v1) / cross(v4−v2, v3−v1)
    num = P[:num]::TMT; den = P[:den]::TMT
    tm_sub!(tmp, V[4, 1], V[2, 1]); tm_sub!(tmp2, V[4, 2], V[2, 2])
    tm_mul!(num, tmp, E[1, 2], ws); tm_mul!(tmp3, tmp2, E[1, 1], ws); tm_sub!(num, num, tmp3)
    tm_mul!(den, tmp, ey, ws); tm_mul!(tmp3, tmp2, ex, ws); tm_sub!(den, den, tmp3)
    tm_inv!(P[:inv]::TMT, den, aw) || return (; ok = false, why = "Y* denominator", M = P[:M], K = P[:K], A2 = P[:A2])
    ystar = P[:ystar]::TMT
    tm_mul!(ystar, num, P[:inv]::TMT, ws)
    Y0 = isArb ? Float64(Arblib.midref(ystar.c[1])) : Float64(ystar.c[1])
    hY = P[:hlin]::TMT; hY2 = P[:hsq]::TMT
    tm_set!(hY, ystar)
    if isArb
        Arblib.sub!(hY.c[1], hY.c[1], Arb(Y0; prec = PREC))
    else
        hY.c[1] -= Y0
    end
    K = TAB[].order
    hpow = TMT[newtm(T) for _ = 1:K]
    tm_set!(hpow[1], hY)
    for j = 2:K; tm_mul!(hpow[j], hpow[j-1], hY, ws); end
    Ik = Ik_taylor(T, m, Y0)
    Iv = P[:I]::Vector{TMT}
    if isArb
        Hm = Mag(0); polybound!(Hm, hY, ws); Arblib.add!(Hm, Hm, hY.r[1])
        for k = 1:m; tm_poly_at!(Iv[k], Ik[k], hpow, Hm, ws); end
    else
        H = polybound(hY) + hY.r[1]
        for k = 1:m; tm_poly_at!(Iv[k], Ik[k], hpow, H); end
    end
    # ---- polynomial coefficient matrices (families: base + first-order p-derivative polynomials)
    idx = basis_index(n); pidx = psi_index(kpsi)
    hasd = hasproperty(sol, :dC)
    dc = hasd ? p0 .- sol.p0 : zeros(4)            # box centre relative to the solve point (0 in the certified run)
    Cfam = Vector{Vector{Matrix{T}}}(undef, 4)
    for k = 1:3
        base = copy(sol.C[:, k])
        hasd && (for i = 1:4; base .+= dc[i] .* sol.dC[:, k, i]; end)
        Cfam[k] = [coeffmat(T, base, idx, n)]
        hasd && append!(Cfam[k], [coeffmat(T, sol.dC[:, k, i], idx, n) for i = 1:4])
    end
    basep = copy(sol.dpsi)
    hasd && (for i = 1:4; basep .+= dc[i] .* sol.ddpsi[:, i]; end)
    Cfam[4] = [coeffmat(T, basep, pidx, kpsi)]
    hasd && append!(Cfam[4], [coeffmat(T, sol.ddpsi[:, i], pidx, kpsi) for i = 1:4])
    nder = hasd ? 4 : 0
    Ml = P[:M]::Vector{TMT}; Kl = P[:K]::Vector{TMT}; A2 = P[:A2]::TMT
    for x in Ml; tm_zero!(x); end
    for x in Kl; tm_zero!(x); end
    tm_zero!(A2)
    _assemble_nodes!(aw, V, E, A, G, ell, ex, ey, Iv, Cfam, nder, rad, Ml, Kl, A2, n, kpsi, Kd, m, NCst, sc, sol.rho0)
    # A2 currently holds T = ‖τ + curl ψ‖².  Exact identity (u1 ∈ H¹0, div τ = −r/ρ0, div curl ψ = 0):
    #   ‖∇u1/ρ0 + τ + curl ψ‖² = (2ρ0 A1 − A0)/ρ0² + T.
    if isArb
        r0 = Arb(sol.rho0; prec = PREC); ir2 = 1 / (r0 * r0)
        tm_axpy!(A2, 2 * r0 * ir2, Ml[1], ws); tm_axpy!(A2, -ir2, Kl[1], ws)
    else
        r0 = sol.rho0
        tm_axpy!(A2, 2 / r0, Ml[1]); tm_axpy!(A2, -1 / r0^2, Kl[1])
    end
    return (; ok = true, why = "", M = Ml, K = Kl, A2 = A2)
end

function _assemble_nodes!(aw::AsmWork{T,R}, V::Matrix{TM{T,R}}, E::Matrix{TM{T,R}}, A::Matrix{TM{T,R}}, G::Vector{TM{T,R}},
                          ell::Vector{TM{T,R}}, ex::TM{T,R}, ey::TM{T,R}, Iv::Vector{TM{T,R}}, Cfam::Vector{Vector{Matrix{T}}}, nder::Int, rad::Vector{Float64},
                          Ml::Vector{TM{T,R}}, Kl::Vector{TM{T,R}}, A2::TM{T,R}, n::Int, kpsi::Int, Kd::Int, m::Int,
                          NCst::NodeConsts{T}, sc::Scaling, rho0f::Float64) where {T,R}
    TMT = TM{T,R}
    P = aw.pool; ws = aw.ws
    isArb = T == Arb
    # ---- scalars
    isx = isArb ? 1 / Arb(sc.sx; prec = PREC) : 1 / sc.sx
    isy = isArb ? 1 / Arb(sc.sy; prec = PREC) : 1 / sc.sy
    misx = -isx
    isx2 = isx * isx; isy2 = isy * isy
    txc = tconv(T, sc.xc); tyc = tconv(T, sc.yc)
    rho0 = tconv(T, rho0f); irho0 = isArb ? 1 / rho0 : 1 / rho0f
    two = tconv(T, 2.0); mone = tconv(T, -1.0); half = tconv(T, 0.5)
    # ---- work arrays
    PU = Matrix{T}(undef, Kd + 1, Kd + 1); PV = similar(PU); Cb = similar(PU)
    for x in eachindex(PU); PU[x] = tconv(T, 0.0); PV[x] = tconv(T, 0.0); Cb[x] = tconv(T, 0.0); end
    Dt = [[Matrix{T}(undef, Kd + 1, Kd + 1) for _ = 1:1+nder] for _ = 1:4]
    for DD in Dt, D in DD, x in eachindex(D); D[x] = tconv(T, 0.0); end
    AD = [[Matrix{R}(undef, Kd + 1, Kd + 1) for _ = 1:1+nder] for _ = 1:4]
    for DD in AD, D in DD, x in eachindex(D); D[x] = isArb ? Mag(0) : 0.0; end
    radT = T[tconv(T, rad[i]) for i = 1:4]; radM = [Mag(rad[i]) for i = 1:4]
    ftmp = newtm(T)
    hu = [tconv(T, 0.0) for _ = 1:4]; hv = [tconv(T, 0.0) for _ = 1:4]
    Kord = TAB[].order
    HUV = [newtm(T) for _ = 1:Kord+1, _ = 1:Kord+1]
    invfact = [k == 0 ? tconv(T, 1.0) : (isArb ? 1 / Arb(k; prec = PREC) : 1.0 / k) for k = 0:Kord]
    HU = R[isArb ? Mag(0) : 0.0 for _ = 1:Kd+1]; HV = R[isArb ? Mag(0) : 0.0 for _ = 1:Kd+1]
    tmpm = Mag(0); HuM = Mag(0); HvM = Mag(0)
    x0 = tconv(T, 0.0); y0 = tconv(T, 0.0); u0 = tconv(T, 0.0); v0 = tconv(T, 0.0)
    b = P[:b]::TMT; bx = P[:bx]::TMT; by = P[:by]::TMT; lapb = P[:lapb]::TMT
    L12 = P[:L12]::TMT; L34 = P[:L34]::TMT; L13 = P[:L13]::TMT; L24 = P[:L24]::TMT; L14 = P[:L14]::TMT; L23 = P[:L23]::TMT
    t1 = P[:t1]::TMT; t2 = P[:t2]::TMT; t3 = P[:t3]::TMT; t4 = P[:t4]::TMT; tmp = P[:tmp]::TMT
    p = (P[:p1]::TMT, P[:p2]::TMT, P[:p3]::TMT); pu = (P[:p1u]::TMT, P[:p2u]::TMT, P[:p3u]::TMT); pv = (P[:p1v]::TMT, P[:p2v]::TMT, P[:p3v]::TMT)
    u = (P[:u1]::TMT, P[:u2]::TMT, P[:u3]::TMT); ux = (P[:ux1]::TMT, P[:ux2]::TMT, P[:ux3]::TMT); uy = (P[:uy1]::TMT, P[:uy2]::TMT, P[:uy3]::TMT)
    p1uu = P[:p1uu]::TMT; p1vv = P[:p1vv]::TMT; psu = P[:psu]::TMT; psv = P[:psv]::TMT; lap1 = P[:lap1]::TMT
    Wu = P[:Wu]::TMT; Wux = P[:Wux]::TMT; Wuy = P[:Wuy]::TMT; pre = P[:pre]::TMT; phi = P[:phi]::TMT
    wxT = P[:wx]::TMT; wyT = P[:wy]::TMT; Ww = P[:Ww]::TMT
    rline = P[:rline]::Vector{TMT}; gxline = P[:gxl]::Vector{TMT}; gyline = P[:gyl]::Vector{TMT}; Wline = P[:Wl]::Vector{TMT}
    for (t, (i0, i1, i2, sg)) in enumerate(TRIS)
        Jfac = t == 1 ? A[1, 3] : A[3, 1]        # c_1 (ccw triangle) and c_3 = −cross(v4−v1, v3−v1)
        for i = 1:m
            for j = 1:m
                lam0 = NCst.lam0[i, j]; xi = NCst.xi[i]; eta = NCst.eta[i, j]
                # node centre and displacement
                if isArb
                    Arblib.mul!(x0, lam0, V[i0, 1].c[1]); Arblib.addmul!(x0, xi, V[i1, 1].c[1]); Arblib.addmul!(x0, eta, V[i2, 1].c[1])
                    Arblib.mul!(y0, lam0, V[i0, 2].c[1]); Arblib.addmul!(y0, xi, V[i1, 2].c[1]); Arblib.addmul!(y0, eta, V[i2, 2].c[1])
                    Arblib.sub!(u0, x0, txc); Arblib.mul!(u0, u0, isx)
                    Arblib.sub!(v0, y0, tyc); Arblib.mul!(v0, v0, isy)
                    for q = 1:4
                        Arblib.mul!(hu[q], lam0, V[i0, 1].c[1+q]); Arblib.addmul!(hu[q], xi, V[i1, 1].c[1+q]); Arblib.addmul!(hu[q], eta, V[i2, 1].c[1+q])
                        Arblib.mul!(hu[q], hu[q], isx)
                        Arblib.mul!(hv[q], lam0, V[i0, 2].c[1+q]); Arblib.addmul!(hv[q], xi, V[i1, 2].c[1+q]); Arblib.addmul!(hv[q], eta, V[i2, 2].c[1+q])
                        Arblib.mul!(hv[q], hv[q], isy)
                    end
                    Arblib.zero!(HuM); Arblib.zero!(HvM)
                    for q = 1:4
                        Arblib.get!(tmpm, hu[q]); Arblib.add!(HuM, HuM, tmpm)
                        Arblib.get!(tmpm, hv[q]); Arblib.add!(HvM, HvM, tmpm)
                    end
                    Arblib.one!(HU[1]); Arblib.one!(HV[1])
                    for a = 1:Kd
                        Arblib.mul!(HU[a+1], HU[a], HuM); Arblib.div!(HU[a+1], HU[a+1], UInt(a))
                        Arblib.mul!(HV[a+1], HV[a], HvM); Arblib.div!(HV[a+1], HV[a+1], UInt(a))
                    end
                else
                    x0 = lam0 * V[i0, 1].c[1] + xi * V[i1, 1].c[1] + eta * V[i2, 1].c[1]
                    y0 = lam0 * V[i0, 2].c[1] + xi * V[i1, 2].c[1] + eta * V[i2, 2].c[1]
                    u0 = (x0 - txc) * isx; v0 = (y0 - tyc) * isy
                    for q = 1:4
                        hu[q] = (lam0 * V[i0, 1].c[1+q] + xi * V[i1, 1].c[1+q] + eta * V[i2, 1].c[1+q]) * isx
                        hv[q] = (lam0 * V[i0, 2].c[1+q] + xi * V[i1, 2].c[1+q] + eta * V[i2, 2].c[1+q]) * isy
                    end
                    Hu = sum(abs, hu); Hv = sum(abs, hv)
                    HU[1] = 1.0; HV[1] = 1.0
                    for a = 1:Kd; HU[a+1] = HU[a] * Hu / a; HV[a+1] = HV[a] * Hv / a; end
                end
                build_huv!(HUV, hu, hv, ws, invfact)
                legendre_table!(PU, u0, Kd); legendre_table!(PV, v0, Kd)
                for f = 1:4
                    df = f <= 3 ? n : kpsi
                    for q = 1:1+nder
                        dtable!(Dt[f][q], Cb, Cfam[f][q], PU, PV, df); absdtable!(AD[f][q], Dt[f][q], df)
                    end
                end
                if isArb
                    for k = 1:3
                        family_tm!(p[k], Dt[k], AD[k], 0, 0, n, HUV, HU, HV, radT, nder, ftmp, tmpm, radM)
                        family_tm!(pu[k], Dt[k], AD[k], 1, 0, n, HUV, HU, HV, radT, nder, ftmp, tmpm, radM)
                        family_tm!(pv[k], Dt[k], AD[k], 0, 1, n, HUV, HU, HV, radT, nder, ftmp, tmpm, radM)
                    end
                    family_tm!(p1uu, Dt[1], AD[1], 2, 0, n, HUV, HU, HV, radT, nder, ftmp, tmpm, radM)
                    family_tm!(p1vv, Dt[1], AD[1], 0, 2, n, HUV, HU, HV, radT, nder, ftmp, tmpm, radM)
                    family_tm!(psu, Dt[4], AD[4], 1, 0, kpsi, HUV, HU, HV, radT, nder, ftmp, tmpm, radM)
                    family_tm!(psv, Dt[4], AD[4], 0, 1, kpsi, HUV, HU, HV, radT, nder, ftmp, tmpm, radM)
                else
                    for k = 1:3
                        family_tm!(p[k], Dt[k], AD[k], 0, 0, n, HUV, HU, HV, radT, nder, ftmp)
                        family_tm!(pu[k], Dt[k], AD[k], 1, 0, n, HUV, HU, HV, radT, nder, ftmp)
                        family_tm!(pv[k], Dt[k], AD[k], 0, 1, n, HUV, HU, HV, radT, nder, ftmp)
                    end
                    family_tm!(p1uu, Dt[1], AD[1], 2, 0, n, HUV, HU, HV, radT, nder, ftmp)
                    family_tm!(p1vv, Dt[1], AD[1], 0, 2, n, HUV, HU, HV, radT, nder, ftmp)
                    family_tm!(psu, Dt[4], AD[4], 1, 0, kpsi, HUV, HU, HV, radT, nder, ftmp)
                    family_tm!(psv, Dt[4], AD[4], 0, 1, kpsi, HUV, HU, HV, radT, nder, ftmp)
                end
                # scale derivatives to x,y: p_x = p_u/sx, p_y = p_v/sy
                for k = 1:3
                    tm_scale!(pu[k], pu[k], isx, ws); tm_scale!(pv[k], pv[k], isy, ws)
                end
                # ΔP1 = p1uu/sx² + p1vv/sy²  (into p1uu)
                tm_scale!(p1uu, p1uu, isx2, ws); tm_axpy!(p1uu, isy2, p1vv, ws)
                # ---- edge functions at the node (exact quadratics) and the bubble
                for k = 1:4
                    tm_zero!(ell[k])
                    tm_axpy!(ell[k], lam0, A[k, i0], ws); tm_axpy!(ell[k], xi, A[k, i1], ws); tm_axpy!(ell[k], eta, A[k, i2], ws)
                end
                tm_mul!(L12, ell[1], ell[2], ws); tm_mul!(L34, ell[3], ell[4], ws)
                tm_mul!(L13, ell[1], ell[3], ws); tm_mul!(L24, ell[2], ell[4], ws)
                tm_mul!(L14, ell[1], ell[4], ws); tm_mul!(L23, ell[2], ell[3], ws)
                tm_mul!(b, L12, L34, ws)
                tm_mul!(t1, L34, ell[2], ws); tm_mul!(t2, L34, ell[1], ws)
                tm_mul!(t3, L12, ell[4], ws); tm_mul!(t4, L12, ell[3], ws)
                # bx = Σ g_kx t_k = Σ (−e_ky) t_k ;  by = Σ e_kx t_k
                tm_mul!(bx, E[1, 2], t1, ws); tm_addmul!(bx, E[2, 2], t2, ws); tm_addmul!(bx, E[3, 2], t3, ws); tm_addmul!(bx, E[4, 2], t4, ws)
                tm_scale!(bx, bx, mone, ws)
                tm_mul!(by, E[1, 1], t1, ws); tm_addmul!(by, E[2, 1], t2, ws); tm_addmul!(by, E[3, 1], t3, ws); tm_addmul!(by, E[4, 1], t4, ws)
                # Δb = 2 Σ_{i<j} (g_i·g_j) Π_{k≠i,j} l_k
                tm_mul!(lapb, G[1], L34, ws); tm_addmul!(lapb, G[2], L12, ws); tm_addmul!(lapb, G[3], L24, ws)
                tm_addmul!(lapb, G[4], L13, ws); tm_addmul!(lapb, G[5], L23, ws); tm_addmul!(lapb, G[6], L14, ws)
                tm_scale!(lapb, lapb, two, ws)
                # ---- test functions
                for k = 1:3
                    tm_mul!(u[k], b, p[k], ws)
                    tm_mul!(ux[k], bx, p[k], ws); tm_addmul!(ux[k], b, pu[k], ws)
                    tm_mul!(uy[k], by, p[k], ws); tm_addmul!(uy[k], b, pv[k], ws)
                end
                tm_mul!(lap1, lapb, p[1], ws)
                tm_mul!(tmp, bx, pu[1], ws); tm_addmul!(tmp, by, pv[1], ws); tm_axpy!(lap1, two, tmp, ws)
                tm_addmul!(lap1, b, p1uu, ws)
                # r = Δu1 + ρ0 u1  (stored per line)
                rl = rline[j]; tm_set!(rl, lap1); tm_axpy!(rl, rho0, u[1], ws)
                # correction flux part: curl ψ = (ψ_v/sy, −ψ_u/sx)   (the ∇u1/ρ0 part is handled by the identity
                # ‖w‖² = (2ρ0 A1 − A0)/ρ0² + ‖τ + curl ψ‖², see below)
                gxl = gxline[j]; gyl = gyline[j]
                tm_scale!(gxl, psv, isy, ws)
                tm_scale!(gyl, psu, misx, ws)
                # weight W = w_ij · Jfac
                Wl = Wline[j]; tm_scale!(Wl, Jfac, NCst.w[i, j], ws)
                # ---- accumulate M, K
                for k = 1:3
                    tm_mul!(Wu, Wl, u[k], ws); tm_mul!(Wux, Wl, ux[k], ws); tm_mul!(Wuy, Wl, uy[k], ws)
                    for (q, (a, c)) in enumerate(KL)
                        a == k || continue
                        tm_addmul!(Ml[q], Wu, u[c], ws)
                        tm_addmul!(Kl[q], Wux, ux[c], ws); tm_addmul!(Kl[q], Wuy, uy[c], ws)
                    end
                end
            end
            # ---- line integration along the Duffy line i: φ_j = −(1/ρ0)(1−X_i)(Σ_k S_jk r_k − Σ_k I_k r_k)
            tm_zero!(pre)
            for k = 1:m; tm_addmul!(pre, Iv[k], rline[k], ws); end
            fac = -(NCst.omx[i] * irho0)
            for j = 1:m
                tm_scale!(phi, pre, mone, ws)
                for k = 1:m; tm_axpy!(phi, NCst.S[j, k], rline[k], ws); end
                tm_scale!(phi, phi, fac, ws)
                tm_set!(wxT, gxline[j]); tm_addmul!(wxT, phi, ex, ws)
                tm_set!(wyT, gyline[j]); tm_addmul!(wyT, phi, ey, ws)
                tm_mul!(Ww, Wline[j], wxT, ws); tm_addmul!(A2, Ww, wxT, ws)
                tm_mul!(Ww, Wline[j], wyT, ws); tm_addmul!(A2, Ww, wyT, ws)
            end
        end
    end
    nothing
end
