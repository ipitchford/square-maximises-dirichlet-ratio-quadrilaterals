# bound1d.jl — rigorous 1-D reduction lower bound for λ1 of all Q_p in a box (box constant).
#
# Lemma (separation of variables).  Let u be a unit vector, x = ⟨·,u⟩, y = ⟨·,u⊥⟩, Q convex, h(x) = length of the
# slice Q ∩ {x = const}, I ⊇ proj_u(Q).  For φ ∈ H¹₀(Q): ∫|∂_yφ|²dy ≥ π²/h(x)² ∫φ²dy on each slice (1-D Poincaré),
# and v(x) = (∫φ(x,y)²dy)^{1/2} ∈ H¹₀(I) with |v'(x)| ≤ (∫|∂_xφ|²dy)^{1/2} (φ = 0 at the slice ends).  Hence
#     ∫_Q |∇φ|² ≥ ∫_I (v'² + V v²) dx ≥ μ1(V, I) ∫_Q φ²   for every V with 0 ≤ V ≤ π²/h² on proj_u(Q),
# i.e. λ1(Q) ≥ μ1(−d²/dx² + V, Dirichlet on I).  We take I = [x0, xK] ⊇ ∪_p proj(Q_p), K slabs J_k, and the step
# potential V_k = π²/H_k² with H_k ≥ sup_{p ∈ box} sup_{x ∈ J_k} h_p(x) (slab part of Q_p ⊂ hull of the vertices in
# the slab and the edge/slab-boundary intersections; all in ball arithmetic).
# μ1 > μ is certified by Sturm oscillation: the solution of −v'' + (V − μ)v = 0, v(x0) = 0, v'(x0) = 1, has no zero in
# (x0, xK] (exact transfer matrices on each slab; hyperbolic slabs: at most one zero, so positivity at both ends
# suffices; trigonometric slabs: phase condition).

# slab upper bounds of the slice length, for direction (cs, sn) (Float64 numbers; normalised in Arb)
function slab_heights(p0, rad, cs, sn, K)
    nrm = sqrt(Arb(cs)^2 + Arb(sn)^2)
    ux = Arb(cs) / nrm; uy = Arb(sn) / nrm
    # vertex coordinates over the box (exact per-vertex ranges)
    V = [let c = Arb(V0[j, k]; prec = PREC)
             for i = 1:4; c += Arb(DV[j, k, i]; prec = PREC) * Arb((p0[i] - rad[i], p0[i] + rad[i]); prec = PREC); end
             c
         end for j = 1:4, k = 1:2]
    X = [V[j, 1] * ux + V[j, 2] * uy for j = 1:4]
    Y = [-V[j, 1] * uy + V[j, 2] * ux for j = 1:4]
    x0 = minimum(Float64(Arblib.lbound(x), RoundDown) for x in X); xK = maximum(Float64(Arblib.ubound(x), RoundUp) for x in X)
    xi = [x0 + (xK - x0) * k / K for k = 0:K]; xi[end] = xK
    H = zeros(K)
    hullhi(a, b) = max(Arblib.ubound(a), Arblib.ubound(b)); hulllo(a, b) = min(Arblib.lbound(a), Arblib.lbound(b))
    for k = 1:K
        a = xi[k]; b = xi[k+1]
        ymax = Arf(-Inf); ymin = Arf(Inf)
        for j = 1:4                                      # vertices that may lie in the slab
            if Arblib.ubound(X[j]) >= a && Arblib.lbound(X[j]) <= b
                ymax = max(ymax, Arblib.ubound(Y[j])); ymin = min(ymin, Arblib.lbound(Y[j]))
            end
        end
        for j = 1:4, ξ in (a, b)                         # edge ∩ slab boundary
            jn = mod1(j + 1, 4)
            xlo = min(Arblib.lbound(X[j]), Arblib.lbound(X[jn])); xhi = max(Arblib.ubound(X[j]), Arblib.ubound(X[jn]))
            (xlo <= ξ <= xhi) || continue
            dx = X[jn] - X[j]
            yhi = hullhi(Y[j], Y[jn]); ylo = hulllo(Y[j], Y[jn])
            if !Arblib.contains_zero(dx)
                t = (Arb(ξ) - X[j]) / dx
                t = Arb((max(Arblib.lbound(t), Arf(0)), min(Arblib.ubound(t), Arf(1))))
                yv = Y[j] + t * (Y[jn] - Y[j])
                yhi = min(yhi, Arblib.ubound(yv)); ylo = max(ylo, Arblib.lbound(yv))
            end
            ymax = max(ymax, yhi); ymin = min(ymin, ylo)
        end
        H[k] = ymax > ymin ? Float64(Arblib.ubound(Arb(ymax) - Arb(ymin)), RoundUp) : 0.0
    end
    return xi, H
end

# Sturm check: does the solution with parameter μ stay positive on (x0, xK]?
function no_zero(xi, H, μ::Arb)
    v = Arb(0); dv = Arb(1); started = false
    pi2 = PI_A^2
    for k = 1:length(H)
        L = Arb(xi[k+1]) - Arb(xi[k])
        Vk = H[k] > 0 ? pi2 / Arb(H[k])^2 : Arb(1e30)
        q = Vk - μ
        if Arblib.ispositive(q)                          # hyperbolic: v = v cosh(κs) + dv/κ sinh(κs)
            κ = sqrt(q); c = cosh(κ * L); s = sinh(κ * L)
            vn = v * c + dv * s / κ; dvn = v * κ * s + dv * c
            # at most one zero on the slab: positive at both ends (or zero at x0 with dv > 0) ⇒ positive inside
            Arblib.ispositive(vn) || return false
            (started ? Arblib.ispositive(v) : Arblib.ispositive(dv)) || return false
            v = vn; dv = dvn
        elseif Arblib.isnegative(q)                      # trigonometric: v = A cos(ω s − φ)
            ω = sqrt(-q); c = cos(ω * L); s = sin(ω * L)
            vn = v * c + dv * s / ω; dvn = -v * ω * s + dv * c
            # no zero in (0, L]: phase φ = atan2(dv/ω, v) ∈ (−π/2, π] (v ≥ 0 at the start); first zero at s* = (π/2 + φ)/ω
            if !started
                # v(0) = 0, dv > 0: v = (dv/ω) sin(ωs) > 0 on (0, L] iff ωL < π
                Arblib.ispositive(dv) || return false
                Arblib.ispositive(PI_A - ω * L) || return false
            else
                Arblib.ispositive(v) || return false
                φ = atan(dv / ω, v)
                Arblib.ispositive(PI_A / 2 + φ - ω * L) || return false
            end
            Arblib.ispositive(vn) || return false
            v = vn; dv = dvn
        else
            return false                                 # q straddles 0: refine μ instead
        end
        started = true
    end
    return true
end

"""
    lambda1_1d(p0, rad; K = 48, dirs) -> Float64

Rigorous lower bound of λ1(Q_p) for all convex Q_p in the box p0 ± rad (1-D reduction, best of the directions).
"""
function lambda1_1d(p0, rad; K::Int = 48, dirs = nothing)
    Vf = verts(p0)
    if dirs === nothing
        dirs = Float64[]
        for i = 1:4, j = i+1:4                         # directions of all vertex pairs (diameter, diagonals, edges)
            push!(dirs, atan(Vf[j, 2] - Vf[i, 2], Vf[j, 1] - Vf[i, 1]))
        end
    end
    best = 0.0
    for th in dirs
        xi, H = slab_heights(p0, rad, cos(th), sin(th), K)
        maximum(H) > 0 || continue
        # bisection on μ in [0, μhi]; μhi from the largest slab potential bound
        lo = 0.0; hi = maximum(π^2 / h^2 for h in H if h > 0) + π^2 / (xi[end] - xi[1])^2
        no_zero(xi, H, Arb(lo)) || continue
        for _ = 1:30
            mid = (lo + hi) / 2
            no_zero(xi, H, Arb(mid)) ? (lo = mid) : (hi = mid)
        end
        best = max(best, lo)
    end
    return best
end
