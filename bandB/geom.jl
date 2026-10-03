# geom.jl — exact rational geometry for the Lemma B′ band cover (BAND_SPEC.md §1–§3).  Author: Anonymous.
#
# Chart I  (chart = 1), coordinates p = (y3, t, x1, x3):  P0 = (0,0), P1 = (x1, −t·y3), P2 = (1,0), P3 = (x3, y3).
# Chart II (chart = 2), coordinates p = (y2, t, x2, x3):  P0 = (0,0), P1 = (1,0), P2 = (x2, y2), P3 = (x3, t·y2).
# Every vertex coordinate and every c_i = (P_{i+1} − P_i) × (P_{i+2} − P_i) is MULTILINEAR in p (degree ≤ 1 in each
# coordinate separately; checked symbolically in sympy_check.py), so its exact range over a box is [min, max] of its
# 16 corner values.  All exclusion tests below are decided exactly in rationals with this fact.

const QQ = Rational{BigInt}

struct Box
    chart::Int
    lo::NTuple{4,QQ}
    hi::NTuple{4,QQ}
end

const ROOT1 = Box(1, (QQ(1, 20), QQ(0), QQ(0), QQ(0)), (QQ(2, 5), QQ(1), QQ(1), QQ(1, 2)))
const ROOT2 = Box(2, (QQ(1, 10), QQ(0), QQ(0), QQ(0)), (QQ(2, 5), QQ(1), QQ(1), QQ(1)))
root_box(c::Int) = c == 1 ? ROOT1 : ROOT2
const COORD_NAMES = Dict(1 => ("y3", "t", "x1", "x3"), 2 => ("y2", "t", "x2", "x3"))

volume(B::Box) = prod(B.hi[k] - B.lo[k] for k = 1:4)

"bisect B in coordinate d (1-based) at the exact midpoint"
function bisect(B::Box, d::Int)
    m = (B.lo[d] + B.hi[d]) / 2
    lower = Box(B.chart, B.lo, ntuple(k -> k == d ? m : B.hi[k], 4))
    upper = Box(B.chart, ntuple(k -> k == d ? m : B.lo[k], 4), B.hi)
    return lower, upper
end

# Path encoding: one character per bisection, 'a'+(d−1) = lower half of coordinate d, 'A'+(d−1) = upper half.
pathchar(d::Int, upper::Bool) = upper ? Char('A' + d - 1) : Char('a' + d - 1)
function box_from_path(chart::Int, path::AbstractString)
    B = root_box(chart)
    for ch in path
        if 'a' <= ch <= 'd'
            B = bisect(B, Int(ch - 'a') + 1)[1]
        elseif 'A' <= ch <= 'D'
            B = bisect(B, Int(ch - 'A') + 1)[2]
        else
            error("bad path character $ch")
        end
    end
    return B
end

corners(B::Box) = [ntuple(k -> ((m >> (k - 1)) & 1 == 1 ? B.hi[k] : B.lo[k]), 4) for m = 0:15]

"vertices P0..P3 (counter-clockwise for valid chart points) at the exact point p"
function verts(chart::Int, p)
    z = zero(QQ); o = one(QQ)
    if chart == 1
        y3, t, x1, x3 = p
        return ((z, z), (x1, -t * y3), (o, z), (x3, y3))
    else
        y2, t, x2, x3 = p
        return ((z, z), (o, z), (x2, y2), (x3, t * y2))
    end
end
cross2(a, b) = a[1] * b[2] - a[2] * b[1]
sub2(a, b) = (a[1] - b[1], a[2] - b[2])
"c_i = (P_{i+1} − P_i) × (P_{i+2} − P_i), returned for paper indices i = 0..3 at Julia positions 1..4"
cvals(V) = ntuple(i -> cross2(sub2(V[mod1(i + 1, 4)], V[i]), sub2(V[mod1(i + 2, 4)], V[i])), 4)

const ETA = QQ(7, 200)      # Lemma A threshold (paper §2.2: ξ ≤ (7/3)(1+η)² = 2.4995… < 5/2)
const S_LO = QQ(1, 10)      # Lemma B threshold
const S_HI = QQ(2, 5)       # top of the band

sqmin(lo, hi) = (lo <= 0 <= hi) ? zero(lo) : min(lo^2, hi^2)   # min of g² for g ∈ [lo, hi]
sqmax(lo, hi) = max(lo^2, hi^2)                                 # max of g² for g ∈ [lo, hi]

"""
    exclusion(B) -> (status::Symbol, reason::String) or nothing

Exact rational decision of the exclusion rules of BAND_SPEC §2 on the WHOLE closed box B (single rule per box).
Ranges of multilinear functions are exact (corner values); |P_j − P_k|² is bounded componentwise.
"""
function exclusion(B::Box)
    C = corners(B)
    V = [verts(B.chart, p) for p in C]
    cv = [cvals(v) for v in V]
    cmin(i) = minimum(c[i] for c in cv)
    cmax(i) = maximum(c[i] for c in cv)
    function d2range(j, k)      # (lower, upper) bounds of |P_j − P_k|² over B
        dx = [v[j][1] - v[k][1] for v in V]; dy = [v[j][2] - v[k][2] for v in V]
        ax, bx = extrema(dx); ay, by = extrema(dy)
        return sqmin(ax, bx) + sqmin(ay, by), sqmax(ax, bx) + sqmax(ay, by)
    end
    # OUTSIDE_CHART: no point of B is a valid chart point
    for j = 1:4, k = j+1:4
        d2range(j, k)[1] > 1 && return (:OUTSIDE_CHART, "|P$(j-1)P$(k-1)| > 1")
    end
    if B.chart == 2
        for i = 1:4
            cmax(i) <= 0 && return (:OUTSIDE_CHART, "c$(i-1) <= 0")
        end
        B.lo[4] >= B.hi[3] && return (:OUTSIDE_CHART, "x3 >= x2")
    end
    # LEMMA_A: c_i <= (7/200) c_{i+2} at every point (c_i − η c_{i+2} is multilinear)
    for i = 1:4
        j = mod1(i + 2, 4)
        maximum(c[i] - ETA * c[j] for c in cv) <= 0 && return (:LEMMA_A, "c$(i-1) <= 7/200 c$(j-1)")
    end
    # LEMMA_B: an upper bound of s is <= 1/10
    if B.chart == 1
        B.hi[1] * (1 + B.hi[2]) <= S_LO && return (:LEMMA_B, "y1+y3 <= 1/10")
        for i = 1:4                      # s <= max(c_i, c_{i+3}) / |e_i| (width normal to edge i)
            mc = max(cmax(i), cmax(mod1(i + 3, 4)))
            mc < 0 && continue
            mc^2 <= S_LO^2 * d2range(i, mod1(i + 1, 4))[1] && return (:LEMMA_B, "width normal to edge $(i-1) <= 1/10")
        end
    else
        B.hi[1] <= S_LO && return (:LEMMA_B, "y2 <= 1/10")
    end
    # OUT_OF_BAND: a lower bound of s is > 2/5
    if B.chart == 1
        B.lo[1] > S_HI && return (:OUT_OF_BAND, "y3 > 2/5")
        allok = true
        for i = 1:4                      # s = min_i max(c_i, c_{i+3}) / |e_i|
            m = max(cmin(i), cmin(mod1(i + 3, 4)))
            if !(m > 0 && m^2 > S_HI^2 * d2range(i, mod1(i + 1, 4))[2])
                allok = false; break
            end
        end
        allok && return (:OUT_OF_BAND, "width lower bound > 2/5")
    else
        B.lo[1] > S_HI && return (:OUT_OF_BAND, "y2 > 2/5")
    end
    return nothing
end

# ---------------------------------------------------------------------------------------------------------------
# Piecewise-linear profiles with exact rational breakpoints on [0, 1].
struct PWL
    X::Vector{QQ}     # strictly increasing, X[1] = 0, X[end] = 1
    Y::Vector{QQ}
end
function pwl_eval(p::PWL, x::QQ)
    j = clamp(searchsortedlast(p.X, x), 1, length(p.X) - 1)
    return p.Y[j] + (p.Y[j+1] - p.Y[j]) * (x - p.X[j]) / (p.X[j+1] - p.X[j])
end
"PWL through the given abscissae (deduplicated, sorted) of a function that is linear between them"
function pwl_from(f, xs)
    X = sort(unique(QQ.(xs)))
    @assert X[1] == 0 && X[end] == 1
    return PWL(X, [f(x) for x in X])
end

Tplus(x, l, u) = x < l ? x / l : (x <= u ? one(QQ) : (1 - x) / (1 - u))     # max_{a∈[l,u]} tent(x; a)
Tminus(x, l, u) = min(x / u, (1 - x) / (1 - l))                               # min_{a∈[l,u]} tent(x; a)
tminus_kink(l, u) = u / (1 + u - l)

orient(o, a, b) = (a[1] - o[1]) * (b[2] - o[2]) - (a[2] - o[2]) * (b[1] - o[1])
"upper concave envelope over [0,1] of a point set containing (0,0), (1,0), all x in [0,1]; as a PWL"
function upper_envelope(pts)
    P = sort(unique([(QQ(p[1]), QQ(p[2])) for p in pts]))
    H = Tuple{QQ,QQ}[]
    for p in P
        while length(H) >= 2 && orient(H[end-1], H[end], p) >= 0
            pop!(H)
        end
        push!(H, p)
    end
    X = QQ[]; Y = QQ[]
    for (x, y) in H                       # vertical runs (only possible at x = 0, 1): keep the top
        if !isempty(X) && X[end] == x
            Y[end] = max(Y[end], y)
        else
            push!(X, x); push!(Y, y)
        end
    end
    @assert X[1] == 0 && X[end] == 1
    return PWL(X, Y)
end

"""
    profiles(B) -> (hout, ain, bin)

hout ≥ h_Q (slice length) for every Q in B (concave PWL);  ain ≤ ... : Q_in = {0<x<1, ain < y < bin} ⊂ Q for every
convex Q in B, with ain convex (≤ 0) and bin concave.  Proofs: BAND_README.md §3.
"""
function profiles(B::Box)
    if B.chart == 1
        y3l, tl, x1l, x3l = B.lo; y3u, tu, x1u, x3u = B.hi
        y1l = tl * y3l; y1u = tu * y3u
        hout = pwl_from(x -> y1u * Tplus(x, x1l, x1u) + y3u * Tplus(x, x3l, x3u), [0, x1l, x1u, x3l, x3u, 1])
        k1 = tminus_kink(x1l, x1u); k3 = tminus_kink(x3l, x3u)
        ain = pwl_from(x -> -y1l * Tminus(x, x1l, x1u), [0, k1, 1])
        bin = pwl_from(x -> y3l * Tminus(x, x3l, x3u), [0, k3, 1])
        return hout, ain, bin
    else
        y2l, tl, x2l, x3l = B.lo; y2u, tu, x2u, x3u = B.hi
        y3l = tl * y2l; y3u = tu * y2u
        pts = [(QQ(0), QQ(0)), (QQ(1), QQ(0))]
        for x in (x2l, x2u), y in (y2l, y2u); push!(pts, (x, y)); end
        for x in (x3l, x3u), y in (y3l, y3u); push!(pts, (x, y)); end
        hout = upper_envelope(pts)
        apex(l, u, yl) = (u / (1 + u - l), yl / (1 + u - l))
        bin = upper_envelope([(QQ(0), QQ(0)), (QQ(1), QQ(0)), apex(x2l, x2u, y2l), apex(x3l, x3u, y3l)])
        ain = PWL([QQ(0), QQ(1)], [QQ(0), QQ(0)])
        return hout, ain, bin
    end
end

# ---------------------------------------------------------------------------------------------------------------
# Exact slice of the quadrilateral at a point (for tests): Q = conv(P0..P3); returns (lower, upper) y of Q ∩ {x}.
function slice_exact(V, x::QQ)
    ys = QQ[]
    for i = 1:4
        a = V[i]; b = V[mod1(i + 1, 4)]
        if a[1] == b[1]
            a[1] == x && (push!(ys, a[2]); push!(ys, b[2]))
        elseif min(a[1], b[1]) <= x <= max(a[1], b[1])
            push!(ys, a[2] + (b[2] - a[2]) * (x - a[1]) / (b[1] - a[1]))
        end
    end
    isempty(ys) && return nothing
    return minimum(ys), maximum(ys)
end
"exact convexity (all c_i > 0) and diameter (all distances ≤ 1) test of a chart point"
function valid_point(chart, p)
    V = verts(chart, p)
    all(c -> c > 0, cvals(V)) || return false
    for j = 1:4, k = j+1:4
        d = sub2(V[j], V[k]); d[1]^2 + d[2]^2 <= 1 || return false
    end
    chart == 2 && !(V[4][1] < V[3][1]) && return false
    return true
end
"exact minimal width of a convex quadrilateral, squared: min_i max(c_i, c_{i+3})² / |e_i|²"
function width2_exact(V)
    c = cvals(V)
    minimum(max(c[i], c[mod1(i + 3, 4)])^2 / (sub2(V[mod1(i + 1, 4)], V[i])[1]^2 + sub2(V[mod1(i + 1, 4)], V[i])[2]^2) for i = 1:4)
end
