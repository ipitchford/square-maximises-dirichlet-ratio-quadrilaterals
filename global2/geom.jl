# geom.jl — Endo–Osting parametrisation, D4 action, Float64 geometry helpers.
#
# p = (a,b,c,d).  Vertices (ccw; p = 0 is the unit square [0,1]^2):
#   v1 = (0,0) + ½(a+d−b, −a+d+c)      v2 = (1,0) + ½(−a+d+b, −a−d−c)
#   v3 = (1,1) + ½(−a−d−b,  a−d+c)     v4 = (0,1) + ½( a−d+b,  a+d−c)
# c_i(p) = (v_{i+1}−v_i) × (v_{i+2}−v_i)  (indices mod 4), convex non-degenerate iff all c_i > 0.
# |Q_p| = 1 − a² − d².

const V0 = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
# DV[j, k, i] = ∂ (v_j)_k / ∂ p_i     (j vertex, k = 1 (x) / 2 (y), i parameter)
const DV = let D = zeros(4, 2, 4)
    D[1, 1, :] = [0.5, -0.5, 0.0, 0.5];   D[1, 2, :] = [-0.5, 0.0, 0.5, 0.5]
    D[2, 1, :] = [-0.5, 0.5, 0.0, 0.5];   D[2, 2, :] = [-0.5, 0.0, -0.5, -0.5]
    D[3, 1, :] = [-0.5, -0.5, 0.0, -0.5]; D[3, 2, :] = [0.5, 0.0, 0.5, -0.5]
    D[4, 1, :] = [0.5, 0.5, 0.0, -0.5];   D[4, 2, :] = [0.5, 0.0, -0.5, 0.5]
    D
end

"Float64 vertices, 4×2 matrix (rows = v1..v4)."
function verts(p::AbstractVector{<:Real})
    V = copy(V0)
    for j = 1:4, k = 1:2, i = 1:4
        V[j, k] += DV[j, k, i] * p[i]
    end
    return V
end

cross2(ax, ay, bx, by) = ax * by - ay * bx

"c_i(p), i = 1..4 (Float64)."
function cvals(p)
    V = verts(p)
    [cross2(V[mod1(i+1,4),1]-V[i,1], V[mod1(i+1,4),2]-V[i,2], V[mod1(i+2,4),1]-V[i,1], V[mod1(i+2,4),2]-V[i,2]) for i = 1:4]
end

area(p) = 1 - p[1]^2 - p[4]^2

# D4 generators (xi is invariant)
sigma_h(p) = [p[1], p[2], -p[3], -p[4]]
sigma_v(p) = [p[1], -p[2], p[3], -p[4]]
rho_rot(p) = [-p[1], -p[3], -p[2], p[4]]

"map any p into the fundamental domain F = {a,b,c ≥ 0}"
function to_fundamental(p)
    q = collect(Float64, p)
    if q[1] < 0; q = rho_rot(q); end
    if q[2] < 0; q = sigma_v(q); end
    if q[3] < 0; q = sigma_h(q); end
    return q
end

"width (min over edges of the max vertex distance from the edge line) and diameter (Float64)."
function width_diam(V)
    ws = Float64[]
    for i = 1:4
        ex = V[mod1(i+1,4),1] - V[i,1]; ey = V[mod1(i+1,4),2] - V[i,2]; L = hypot(ex, ey)
        push!(ws, maximum(abs(cross2(ex, ey, V[j,1]-V[i,1], V[j,2]-V[i,2])) / L for j = 1:4))
    end
    D = maximum(hypot(V[i,1]-V[j,1], V[i,2]-V[j,2]) for i = 1:4 for j = i+1:4)
    return minimum(ws), D
end

function angles_deg(V)
    A = Float64[]
    for i = 1:4
        ux = V[mod1(i-1,4),1]-V[i,1]; uy = V[mod1(i-1,4),2]-V[i,2]
        wx = V[mod1(i+1,4),1]-V[i,1]; wy = V[mod1(i+1,4),2]-V[i,2]
        push!(A, rad2deg(acos(clamp((ux*wx+uy*wy)/hypot(ux,uy)/hypot(wx,wy), -1, 1))))
    end
    return A
end

"Float64 classification of a point: :cube, :nonconvex, :neartri, :thin, or :R."
function classify_point(p; cube, eta, s0)
    maximum(abs, p) <= cube && return :cube
    c = cvals(p)
    minimum(c) <= 0 && return :nonconvex
    any(c[i] <= eta * c[mod1(i+2,4)] for i = 1:4) && return :neartri
    w, D = width_diam(verts(p))
    w / D <= s0 && return :thin
    return :R
end
