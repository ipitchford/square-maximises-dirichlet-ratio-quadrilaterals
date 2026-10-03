# bounds_tm.jl — p-DEPENDENT rigorous eigenvalue bounds as Taylor models over a box (pointwise valid):
#   λ1(Q_p) ≥ L0(p)  and  λ2(Q_p) ≥ β(p)   for every convex Q_p, p in the box.
# A box-constant bound has to take the worst corner while the Rayleigh–Ritz value is largest at another one; for
# thin shapes λ1 ≈ π²/w² varies by tens of per cent over a box, so constants destroy the criterion.
# Candidates (each a valid pointwise bound; the criterion may use any one of them):
#   strip_i:  Q_p lies in the strip of width w_i(p) = max_j cross(e_i, v_j − v_i)/|e_i| ⇒ λ1 ≥ π²|e_i|²/c_i(p)²,
#             c_i(p) ≥ max_j cross(e_i, v_j − v_i) on the box (dominant vertex + certified slack);
#   rect_u:   Q_p lies in the rectangle with sides Â(p), B̂(p) along a FIXED direction u (Â, B̂ affine upper bounds of the
#             extents: dominant extreme vertices + certified slack) ⇒ λ1 ≥ π²(1/Â² + 1/B̂²), λ2 ≥ π²·min(4/Â²+1/B̂², 1/Â²+4/B̂²);
#   KS:       λ2 ≥ 2π j01²/|Q_p| = 2π j01²/(1 − a² − d²);
#   AC:       λ2 ≥ L0(p) + 3π²/D̂(p)², D̂² ≥ max_{j<k}|v_j − v_k|² (dominant pair + slack) (Andrews–Clutterbuck, convex).
# Reciprocals: 1/f = (1/F)(Σ_{j≤K}(−u)^j) + (1/F)·U^{K+1}/(1−U)·[−1,1],  u = (f−F)/F, |u| ≤ U < 1 on the box.

"TM reciprocal (generic); returns false if the range of f is not bounded away from 0 relative to its centre value."
function tm_recip!(h::TM{T}, f::TM{T}, ws) where {T}
    K = TAB[].order
    u = newtm(T); up = newtm(T); tmp = newtm(T)
    if T == Arb
        F = Arb(Arblib.midref(f.c[1]); prec = PREC)
        Arblib.iszero(F) && return false
        iF = 1 / F
        tm_set!(u, f); Arblib.sub!(u.c[1], u.c[1], F); tm_scale!(u, u, iF, ws)
        Ubm = Mag(0); polybound!(Ubm, u, ws); Arblib.add!(Ubm, Ubm, u.r[1])
        Ua = Arb(Ubm; prec = PREC)
        Arblib.ubound(Ua) < 1 || return false
        tm_const!(h, Arb(1; prec = PREC)); tm_const!(up, Arb(1; prec = PREC))
        for j = 1:K
            tm_mul_full!(tmp, up, u, ws); tm_set!(up, tmp)
            tm_axpy!(h, Arb(isodd(j) ? -1 : 1; prec = PREC), up, ws)
        end
        rem = Ua^(K + 1) / (1 - Ua)
        rm = Mag(0); Arblib.get!(rm, rem); Arblib.add!(h.r[1], h.r[1], rm)
        tm_scale!(h, h, iF, ws)
    else
        F = f.c[1]; F == 0 && return false
        tm_set!(u, f); u.c[1] -= F; tm_scale!(u, u, 1 / F)
        U = polybound(u) + u.r[1]
        U < 1 || return false
        tm_const!(h, 1.0); tm_const!(up, 1.0)
        for j = 1:K
            tm_mul_full!(tmp, up, u, ws); tm_set!(up, tmp); tm_axpy!(h, (-1.0)^j, up)
        end
        h.r[1] += U^(K + 1) / (1 - U)
        tm_scale!(h, h, 1 / F)
    end
    true
end

# upper bound TM of max_j f_j: the f_j dominant at the box centre plus a certified constant slack
function tm_maxbound!(h::TM{T}, fs::Vector, ws) where {T}
    c0 = [T == Arb ? Float64(Arblib.midref(f.c[1])) : f.c[1] for f in fs]
    jstar = argmax(c0)
    tm_set!(h, fs[jstar])
    slack = 0.0
    d = newtm(T)
    for (j, f) in enumerate(fs)
        j == jstar && continue
        tm_sub!(d, f, fs[jstar]); slack = max(slack, tm_upper(d))
    end
    if slack > 0
        T == Arb ? Arblib.add!(h.c[1], h.c[1], Arb(slack; prec = PREC)) : (h.c[1] += slack)
    end
    h
end
tcp(::Type{T}, x) where {T} = tconv(T, Float64(x))

"""
    bound_tms(T, p0, rad, ws) -> (; L0s, betas)

Lists of (name, TM) candidates over the box p0 ± rad: λ1 lower bounds L0s and λ2 lower bounds betas.
Every TM is a pointwise lower bound valid for all convex Q_p in the box (when its construction succeeded).
"""
function bound_tms(::Type{T}, p0::Vector{Float64}, rad::Vector{Float64}, ws) where {T}
    G = box_geom(T, p0, rad, ws)
    pi2 = T == Arb ? PI_A^2 : π^2
    L0s = Tuple{String,TM}[]; betas = Tuple{String,TM}[]
    tmp = newtm(T); tmp2 = newtm(T)
    # ---- strips: λ1 ≥ π² |e_i|² / c_i²
    for i = 1:4
        e2 = newtm(T); tm_sq2!(e2, G.E[i, 1], G.E[i, 2], ws)
        cr = newtm(T); tm_maxbound!(cr, [G.A[i, mod1(i + 2, 4)], G.A[i, mod1(i + 3, 4)]], ws)
        tm_lower(cr) > 0 || continue
        c2 = newtm(T); tm_mul_full!(c2, cr, cr, ws)
        ic2 = newtm(T); tm_recip!(ic2, c2, ws) || continue
        L = newtm(T); tm_mul_full!(L, e2, ic2, ws); tm_scale!(L, L, pi2, ws)
        push!(L0s, ("strip$i", L))
    end
    # ---- rectangles along fixed directions
    Vf = verts(p0)
    for (kdir, th) in enumerate(candidate_dirs(Vf))
        cs = cos(th); sn = sin(th)
        nrm = T == Arb ? sqrt(Arb(cs)^2 + Arb(sn)^2) : hypot(cs, sn)
        pu = [newtm(T) for _ = 1:4]; pv = [newtm(T) for _ = 1:4]       # projections v_j·u, v_j·u⊥ (exact affine)
        for j = 1:4
            tm_scale!(pu[j], G.V[j, 1], tcp(T, cs), ws); tm_axpy!(pu[j], tcp(T, sn), G.V[j, 2], ws)
            tm_scale!(pv[j], G.V[j, 2], tcp(T, cs), ws); tm_axpy!(pv[j], tcp(T, -sn), G.V[j, 1], ws)
        end
        sides = TM[]
        for proj in (pu, pv)
            mx = newtm(T); tm_maxbound!(mx, proj, ws)
            neg = [(x = newtm(T); tm_scale!(x, q, tcp(T, -1.0), ws); x) for q in proj]
            mn = newtm(T); tm_maxbound!(mn, neg, ws)                   # ≥ −min_j
            ext = newtm(T); tm_add!(ext, mx, mn)                       # ≥ max − min
            tm_scale!(ext, ext, T == Arb ? 1 / nrm : 1 / nrm, ws)
            push!(sides, ext)
        end
        (tm_lower(sides[1]) > 0 && tm_lower(sides[2]) > 0) || continue
        inv2 = TM[]
        okr = true
        for sd in sides
            s2 = newtm(T); tm_mul_full!(s2, sd, sd, ws)
            is2 = newtm(T); okr &= tm_recip!(is2, s2, ws); push!(inv2, is2)
        end
        okr || continue
        l1 = newtm(T); tm_add!(l1, inv2[1], inv2[2]); tm_scale!(l1, l1, pi2, ws)
        push!(L0s, ("rect$kdir", l1))
        b1 = newtm(T); tm_scale!(b1, inv2[1], tcp(T, 4.0), ws); tm_add!(b1, b1, inv2[2])
        b2 = newtm(T); tm_scale!(b2, inv2[2], tcp(T, 4.0), ws); tm_add!(b2, b2, inv2[1])
        # λ2(rect) = π² min(b1, b2): use the branch that is certified ≤ the other on the box
        tm_sub!(tmp, b1, b2)
        if tm_upper(tmp) <= 0
            b = b1
        elseif tm_lower(tmp) >= 0
            b = b2
        else
            continue
        end
        bb = newtm(T); tm_scale!(bb, b, pi2, ws)
        push!(betas, ("rect$kdir", bb))
    end
    # ---- Krahn–Szegő
    at = newtm(T); dt = newtm(T)
    tm_linear!(at, tcp(T, p0[1]), [i == 1 ? tcp(T, rad[1]) : tcp(T, 0.0) for i = 1:4])
    tm_linear!(dt, tcp(T, p0[4]), [i == 4 ? tcp(T, rad[4]) : tcp(T, 0.0) for i = 1:4])
    ar = newtm(T); tm_mul_full!(ar, at, at, ws); tm_mul_full!(tmp, dt, dt, ws); tm_add!(ar, ar, tmp)
    tm_scale!(ar, ar, tcp(T, -1.0), ws); T == Arb ? Arblib.add!(ar.c[1], ar.c[1], Arb(1)) : (ar.c[1] += 1.0)
    if tm_lower(ar) > 0
        iar = newtm(T)
        if tm_recip!(iar, ar, ws)
            ks = newtm(T); tm_scale!(ks, iar, T == Arb ? 2 * PI_A * J01_A^2 : 2π * 2.4048255576957728^2, ws)
            push!(betas, ("KS", ks))
        end
    end
    # ---- Andrews–Clutterbuck: λ2 ≥ L0(p) + 3π²/D̂² (best L0 candidate)
    d2s = TM[]
    for j = 1:4, k = j+1:4
        dx = newtm(T); dy = newtm(T); tm_sub!(dx, G.V[j, 1], G.V[k, 1]); tm_sub!(dy, G.V[j, 2], G.V[k, 2])
        d2 = newtm(T); tm_sq2!(d2, dx, dy, ws); push!(d2s, d2)
    end
    D2 = newtm(T); tm_maxbound!(D2, d2s, ws)
    iD2 = newtm(T)
    gap = nothing
    if !isempty(L0s) && tm_lower(D2) > 0 && tm_recip!(iD2, D2, ws)
        gap = newtm(T); tm_scale!(gap, iD2, T == Arb ? 3 * PI_A^2 : 3π^2, ws)
        for (nm, L) in L0s            # one Andrews–Clutterbuck candidate per λ1 bound (names identify them exactly)
            ac = newtm(T); tm_add!(ac, L, gap)
            push!(betas, ("AC+" * nm, ac))
        end
    end
    if gap === nothing && tm_lower(D2) > 0 && tm_recip!(iD2, D2, ws)
        gap = newtm(T); tm_scale!(gap, iD2, T == Arb ? 3 * PI_A^2 : 3π^2, ws)
    end
    return (; L0s, betas, gap)
end

"""
    boot_tm(cp, base, beta1, gap, ws) -> Union{Nothing,TM}

Pointwise Andrews–Clutterbuck bootstrap: with a valid β1(p) ≤ λ2(Q_p) (constant or Taylor model) and the LG
quantities of the box (g2 = β1A1 − A0 > 0 and g1 = β1A2 − A1 > 0 certified), L1(p) = g2/g1 ≤ λ1(Q_p) pointwise
(Lehmann–Goerisch), hence β2(p) = L1(p) + 3π²/D̂(p)² ≤ λ1(Q_p) + 3π²/diam(Q_p)² ≤ λ2(Q_p) (convex Q_p).
The quotient is homogeneous of degree 0 in the (normalised) base quantities.  Returns β2 as a Taylor model.
"""
function boot_tm(cp, base, beta1, gap, ws)
    gap === nothing && return nothing
    T = eltype(base.M[1].c)
    r = crit_eval!(cp, base.M, base.K, base.A2, ws; mode = :LG, beta = beta1, L0 = 0.0)
    r.g1ok || return nothing
    ig1 = newtm(T); tm_recip!(ig1, cp.g1, ws) || return nothing
    L1 = newtm(T); tm_mul_full!(L1, cp.g2, ig1, ws)
    b2 = newtm(T); tm_add!(b2, L1, gap)
    return b2
end
