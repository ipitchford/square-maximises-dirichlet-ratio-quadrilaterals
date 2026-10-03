# tm.jl — Taylor models of order ORDER in 4 variables over the unit box t ∈ [−1,1]^4.
#
# A TM f is (c, R): coefficients c_γ of the monomials t^γ, |γ| ≤ ORDER (graded order), and a remainder bound
# R ≥ 0.  Semantics (pointwise): for every t in the unit box, f(t) ∈ Σ_γ c_γ t^γ + [−R, R].
# Modes:
#   * Arb mode  (T = Arb, R = Mag): coefficients are Arb balls, remainders Mag (upward rounded).  RIGOROUS.
#   * Float64 mode (T = Float64, R = Float64): the same formulas in floating point.  NOT rigorous (pilot only).
# Product h = f·g: all coefficients with |γ| ≤ ORDER exactly (ball arithmetic); every dropped term
# (|α| + |β| > ORDER) is bounded with |t^γ| ≤ 1:  dropped ≤ Σ_{d1+d2 > ORDER} N_{d1}(f) N_{d2}(g),
# N_d(f) = Σ_{|γ|=d} |c_γ|;  R_h = |f|_B R_g + |g|_B R_f + R_f R_g + dropped,  |f|_B = Σ_d N_d(f).
# Truncation to order 2 (fold!): every coefficient with |γ| ≥ 3 is added to R in absolute value.
# With ball coefficients the product is valid pointwise (ball arithmetic encloses all products of points of the
# input balls; |·| uses upper bounds of absolute values of balls).

using Arblib

const ORDER = Ref(3)
const NV = 4

struct TMTables
    order::Int
    nc::Int
    mono::Vector{NTuple{4,Int}}
    deg::Vector{Int}
    degrange::Vector{UnitRange{Int}}       # indices of degree d monomials: degrange[d+1]
    idx::Dict{NTuple{4,Int},Int}
    pa::Vector{Vector{Int}}                # for output k: input pairs (pa[k][q], pb[k][q])
    pb::Vector{Vector{Int}}
    lin::Vector{Int}                       # index of t_i
    quad::Matrix{Int}                      # index of t_i t_j
    shift::Matrix{Int}                     # shift[k,i] = index of t^γ·t_i (0 if degree > order)
end
function make_tables(K::Int)
    mono = NTuple{4,Int}[]
    for d = 0:K
        for a = d:-1:0, b = d-a:-1:0, c = d-a-b:-1:0
            push!(mono, (a, b, c, d - a - b - c))
        end
    end
    nc = length(mono); deg = [sum(m) for m in mono]
    idx = Dict(m => k for (k, m) in enumerate(mono))
    degrange = [findfirst(==(d), deg):findlast(==(d), deg) for d = 0:K]
    pa = [Int[] for _ = 1:nc]; pb = [Int[] for _ = 1:nc]
    for a = 1:nc, b = 1:nc
        deg[a] + deg[b] <= K || continue
        k = idx[mono[a] .+ mono[b]]
        push!(pa[k], a); push!(pb[k], b)
    end
    lin = [idx[ntuple(k -> k == i ? 1 : 0, 4)] for i = 1:4]
    quad = zeros(Int, 4, 4)
    if K >= 2
        for i = 1:4, j = 1:4
            quad[i, j] = idx[ntuple(k -> (k == i) + (k == j), 4)]
        end
    end
    shift = zeros(Int, nc, 4)
    for k = 1:nc, i = 1:4
        m = ntuple(j -> mono[k][j] + (j == i), 4)
        shift[k, i] = sum(m) <= K ? idx[m] : 0
    end
    TMTables(K, nc, mono, deg, degrange, idx, pa, pb, lin, quad, shift)
end
const TAB = Ref(make_tables(3))

struct TM{T,R}
    c::Vector{T}
    r::Vector{R}      # length 1 (mutable cell)
end
newtm(::Type{Float64}) = TM{Float64,Float64}(zeros(TAB[].nc), [0.0])
newtm(::Type{Arb}) = TM{Arb,Mag}([Arb(0; prec = PREC) for _ = 1:TAB[].nc], [Mag(0)])
tmtype(::Type{Float64}) = TM{Float64,Float64}
tmtype(::Type{Arb}) = TM{Arb,Mag}

struct TMWork{T,R}
    m::Vector{R}       # scratch magnitudes
    t::Vector{T}       # scratch scalars
    nf::Vector{R}      # per-degree norms of f
    ng::Vector{R}      # per-degree norms of g
end
newwork(::Type{Float64}) = TMWork{Float64,Float64}(zeros(16), zeros(8), zeros(8), zeros(8))
newwork(::Type{Arb}) = TMWork{Arb,Mag}([Mag(0) for _ = 1:16], [Arb(0; prec = PREC) for _ = 1:8], [Mag(0) for _ = 1:8], [Mag(0) for _ = 1:8])

# ================================================================ Float64 mode
function tm_zero!(h::TM{Float64}); fill!(h.c, 0.0); h.r[1] = 0.0; h; end
function tm_set!(h::TM{Float64}, f::TM{Float64}); copyto!(h.c, f.c); h.r[1] = f.r[1]; h; end
function tm_const!(h::TM{Float64}, x::Real); fill!(h.c, 0.0); h.c[1] = x; h.r[1] = 0.0; h; end
function tm_add!(h::TM{Float64}, f::TM{Float64}, g::TM{Float64})
    @inbounds for k in eachindex(h.c); h.c[k] = f.c[k] + g.c[k]; end
    h.r[1] = f.r[1] + g.r[1]; h
end
function tm_sub!(h::TM{Float64}, f::TM{Float64}, g::TM{Float64})
    @inbounds for k in eachindex(h.c); h.c[k] = f.c[k] - g.c[k]; end
    h.r[1] = f.r[1] + g.r[1]; h
end
function tm_axpy!(h::TM{Float64}, a::Real, f::TM{Float64}, ws = nothing)        # h += a f
    @inbounds for k in eachindex(h.c); h.c[k] += a * f.c[k]; end
    h.r[1] += abs(a) * f.r[1]; h
end
function tm_scale!(h::TM{Float64}, f::TM{Float64}, a::Real, ws = nothing)        # h = a f
    @inbounds for k in eachindex(h.c); h.c[k] = a * f.c[k]; end
    h.r[1] = abs(a) * f.r[1]; h
end
@inline function _degnorms!(N::Vector{Float64}, f::TM{Float64})
    tab = TAB[]
    @inbounds for d = 0:tab.order
        s = 0.0
        for k in tab.degrange[d+1]; s += abs(f.c[k]); end
        N[d+1] = s
    end
    N
end
polybound(f::TM{Float64}) = sum(abs, f.c)
function _mulrem(f::TM{Float64}, g::TM{Float64}, ws::TMWork{Float64})
    K = TAB[].order
    nf = _degnorms!(ws.nf, f); ng = _degnorms!(ws.ng, g)
    Pf = 0.0; Pg = 0.0
    @inbounds for d = 0:K; Pf += nf[d+1]; Pg += ng[d+1]; end
    drop = 0.0
    @inbounds for d1 = 1:K, d2 = K+1-d1:K
        drop += nf[d1+1] * ng[d2+1]
    end
    return Pf * g.r[1] + Pg * f.r[1] + f.r[1] * g.r[1] + drop
end
function tm_mul!(h::TM{Float64}, f::TM{Float64}, g::TM{Float64}, ws::TMWork{Float64})
    tab = TAB[]
    rem = _mulrem(f, g, ws)
    @inbounds for k = 1:tab.nc
        A = tab.pa[k]; B = tab.pb[k]; s = 0.0
        for q in eachindex(A); s += f.c[A[q]] * g.c[B[q]]; end
        h.c[k] = s
    end
    h.r[1] = rem; h
end
function tm_addmul!(h::TM{Float64}, f::TM{Float64}, g::TM{Float64}, ws::TMWork{Float64})   # h += f g
    tab = TAB[]
    rem = _mulrem(f, g, ws)
    @inbounds for k = 1:tab.nc
        A = tab.pa[k]; B = tab.pb[k]; s = 0.0
        for q in eachindex(A); s += f.c[A[q]] * g.c[B[q]]; end
        h.c[k] += s
    end
    h.r[1] += rem; h
end
tm_addrem!(h::TM{Float64}, x::Real) = (h.r[1] += x; h)
remval(h::TM{Float64}) = h.r[1]
"fold all coefficients of degree ≥ dmax+1 into the remainder (|t^γ| ≤ 1 on the unit box)"
function tm_fold!(h::TM{Float64}, dmax::Int)
    tab = TAB[]
    for d = dmax+1:tab.order, k in tab.degrange[d+1]
        h.r[1] += abs(h.c[k]); h.c[k] = 0.0
    end
    h
end

# ================================================================ Arb mode
function tm_zero!(h::TM{Arb}); for x in h.c; Arblib.zero!(x); end; Arblib.zero!(h.r[1]); h; end
function tm_set!(h::TM{Arb}, f::TM{Arb}); for k in eachindex(h.c); Arblib.set!(h.c[k], f.c[k]); end; Arblib.set!(h.r[1], f.r[1]); h; end
function tm_const!(h::TM{Arb}, x); tm_zero!(h); Arblib.set!(h.c[1], x); h; end
function tm_add!(h::TM{Arb}, f::TM{Arb}, g::TM{Arb})
    for k in eachindex(h.c); Arblib.add!(h.c[k], f.c[k], g.c[k]); end
    Arblib.add!(h.r[1], f.r[1], g.r[1]); h
end
function tm_sub!(h::TM{Arb}, f::TM{Arb}, g::TM{Arb})
    for k in eachindex(h.c); Arblib.sub!(h.c[k], f.c[k], g.c[k]); end
    Arblib.add!(h.r[1], f.r[1], g.r[1]); h
end
function tm_axpy!(h::TM{Arb}, a::Arb, f::TM{Arb}, ws::TMWork{Arb})        # h += a f
    for k in eachindex(h.c); Arblib.addmul!(h.c[k], a, f.c[k]); end
    m = ws.m[1]; Arblib.get!(m, a); Arblib.mul!(m, m, f.r[1]); Arblib.add!(h.r[1], h.r[1], m); h
end
function tm_scale!(h::TM{Arb}, f::TM{Arb}, a::Arb, ws::TMWork{Arb})        # h = a f   (h may alias f)
    for k in eachindex(h.c); Arblib.mul!(h.c[k], f.c[k], a); end
    m = ws.m[1]; Arblib.get!(m, a); Arblib.mul!(h.r[1], f.r[1], m); h
end
function _degnorms!(N::Vector{Mag}, f::TM{Arb}, tmp::Mag)
    tab = TAB[]
    for d = 0:tab.order
        s = N[d+1]; Arblib.zero!(s)
        for k in tab.degrange[d+1]; Arblib.get!(tmp, f.c[k]); Arblib.add!(s, s, tmp); end
    end
    N
end
function polybound!(out::Mag, f::TM{Arb}, ws::TMWork{Arb})
    tmp = ws.m[13]
    Arblib.zero!(out)
    for k in eachindex(f.c); Arblib.get!(tmp, f.c[k]); Arblib.add!(out, out, tmp); end
    out
end
function _mulrem!(out::Mag, f::TM{Arb}, g::TM{Arb}, ws::TMWork{Arb})
    K = TAB[].order
    m = ws.m; tmp = m[8]; Pf = m[9]; Pg = m[14]
    nf = _degnorms!(ws.nf, f, tmp); ng = _degnorms!(ws.ng, g, tmp)
    Arblib.zero!(Pf); Arblib.zero!(Pg)
    for d = 0:K; Arblib.add!(Pf, Pf, nf[d+1]); Arblib.add!(Pg, Pg, ng[d+1]); end
    Arblib.mul!(out, Pf, g.r[1])
    Arblib.mul!(tmp, Pg, f.r[1]); Arblib.add!(out, out, tmp)
    Arblib.mul!(tmp, f.r[1], g.r[1]); Arblib.add!(out, out, tmp)
    for d1 = 1:K, d2 = K+1-d1:K
        Arblib.mul!(tmp, nf[d1+1], ng[d2+1]); Arblib.add!(out, out, tmp)
    end
    out
end
function tm_mul!(h::TM{Arb}, f::TM{Arb}, g::TM{Arb}, ws::TMWork{Arb})     # h must not alias f, g
    tab = TAB[]
    _mulrem!(h.r[1], f, g, ws)
    for k = 1:tab.nc
        A = tab.pa[k]; B = tab.pb[k]; hk = h.c[k]
        Arblib.mul!(hk, f.c[A[1]], g.c[B[1]])
        for q = 2:length(A); Arblib.addmul!(hk, f.c[A[q]], g.c[B[q]]); end
    end
    h
end
function tm_addmul!(h::TM{Arb}, f::TM{Arb}, g::TM{Arb}, ws::TMWork{Arb})  # h += f g  (h must not alias f, g)
    tab = TAB[]
    rm = ws.m[15]
    _mulrem!(rm, f, g, ws)
    Arblib.add!(h.r[1], h.r[1], rm)
    for k = 1:tab.nc
        A = tab.pa[k]; B = tab.pb[k]; hk = h.c[k]
        for q in eachindex(A); Arblib.addmul!(hk, f.c[A[q]], g.c[B[q]]); end
    end
    h
end
tm_addrem!(h::TM{Arb}, x::Mag) = (Arblib.add!(h.r[1], h.r[1], x); h)
remval(h::TM{Arb}) = Float64(Arblib.ubound(Arb(h.r[1]; prec = PREC)), RoundUp)
function tm_fold!(h::TM{Arb}, dmax::Int)
    tab = TAB[]; tmp = Mag(0)
    for d = dmax+1:tab.order, k in tab.degrange[d+1]
        Arblib.get!(tmp, h.c[k]); Arblib.add!(h.r[1], h.r[1], tmp); Arblib.zero!(h.c[k])
    end
    h
end

# ================================================================ range bounds over the unit box
# Quadratic part: exact 1-D range of c_i t_i + c_ii t_i^2 per variable plus |c_ij| for cross terms;
# terms of degree ≥ 3: |c_γ| each; then the remainder.
function _range1d(ci::Float64, cii::Float64)        # [min, max] of ci t + cii t^2 on [-1,1]
    lo = min(ci + cii, -ci + cii); hi = max(ci + cii, -ci + cii)
    if cii != 0
        ts = -ci / (2cii)
        if -1 < ts < 1
            v = ci * ts + cii * ts^2; lo = min(lo, v); hi = max(hi, v)
        end
    end
    return lo, hi
end
function tm_bounds(f::TM{Float64})
    tab = TAB[]
    lo = f.c[1]; hi = f.c[1]
    for i = 1:4
        l, h = _range1d(f.c[tab.lin[i]], tab.order >= 2 ? f.c[tab.quad[i, i]] : 0.0); lo += l; hi += h
    end
    if tab.order >= 2
        for i = 1:4, j = i+1:4
            a = abs(f.c[tab.quad[i, j]]); lo -= a; hi += a
        end
    end
    for d = 3:tab.order, k in tab.degrange[d+1]
        a = abs(f.c[k]); lo -= a; hi += a
    end
    return lo - f.r[1], hi + f.r[1]
end
tm_upper(f::TM{Float64}) = tm_bounds(f)[2]
tm_lower(f::TM{Float64}) = tm_bounds(f)[1]

const UNIT = Arb((-1, 1); prec = PREC)
const UNIT01 = Arb((0, 1); prec = PREC)
"Arb enclosure of the range of the TM over the unit box (rigorous)."
function tm_enclose(f::TM{Arb})
    tab = TAB[]
    acc = Arb(f.c[1]; prec = PREC)
    for i = 1:4
        ci = f.c[tab.lin[i]]
        cii = tab.order >= 2 ? f.c[tab.quad[i, i]] : Arb(0)
        if Arblib.contains_zero(cii)
            acc += ci * UNIT + cii * UNIT01
        else
            e1 = ci + cii; e2 = cii - ci
            lo = min(Arblib.lbound(e1), Arblib.lbound(e2)); hi = max(Arblib.ubound(e1), Arblib.ubound(e2))
            ts = -ci / (2 * cii)
            if Arblib.ubound(ts) > -1 && Arblib.lbound(ts) < 1
                v = -(ci * ci) / (4 * cii)
                lo = min(lo, Arblib.lbound(v)); hi = max(hi, Arblib.ubound(v))
            end
            acc += Arb((lo, hi); prec = PREC)
        end
    end
    if tab.order >= 2
        for i = 1:4, j = i+1:4
            acc += f.c[tab.quad[i, j]] * UNIT
        end
    end
    for d = 3:tab.order, k in tab.degrange[d+1]
        acc += f.c[k] * UNIT
    end
    acc += Arb(f.r[1]; prec = PREC) * UNIT
    return acc
end
tm_upper(f::TM{Arb}) = Float64(Arblib.ubound(tm_enclose(f)), RoundUp)
tm_lower(f::TM{Arb}) = Float64(Arblib.lbound(tm_enclose(f)), RoundDown)

"evaluate the polynomial part at a point t (Float64; testing only)"
function tm_eval(f::TM{Float64}, t)
    tab = TAB[]; s = 0.0
    for k = 1:tab.nc
        m = tab.mono[k]; s += f.c[k] * t[1]^m[1] * t[2]^m[2] * t[3]^m[3] * t[4]^m[4]
    end
    s
end

# ================================================================ full products (per-box use)
# h = f·g with ALL product coefficients computed exactly (degree ≤ 2·ORDER); the coefficients of degree > ORDER
# are then folded into the remainder with their exact (ball) values: |Σ_{deg>K} c_γ t^γ| ≤ Σ |c_γ|.
# This keeps the cancellation inside each product (never worse than the norm bound used by tm_mul!).
struct FullTables
    nc2::Int
    deg2::Vector{Int}
    trip::Vector{NTuple{3,Int}}     # (a, b, k2)
end
function make_full(K::Int)
    t1 = make_tables(K); t2 = make_tables(2K)
    trip = NTuple{3,Int}[]
    for a = 1:t1.nc, b = 1:t1.nc
        push!(trip, (a, b, t2.idx[t1.mono[a] .+ t1.mono[b]]))
    end
    # the first nc monomials of t2 coincide with t1 (same generation order)
    @assert t2.mono[1:t1.nc] == t1.mono
    FullTables(t2.nc, t2.deg, trip)
end
const FTAB = Ref(make_full(3))
set_order!(K::Int) = (ORDER[] = K; TAB[] = make_tables(K); FTAB[] = make_full(K); nothing)

function tm_mul_full!(h::TM{Float64}, f::TM{Float64}, g::TM{Float64}, ws, buf::Vector{Float64} = zeros(FTAB[].nc2))
    ft = FTAB[]; nc = TAB[].nc
    fill!(buf, 0.0)
    @inbounds for (a, b, k) in ft.trip; buf[k] += f.c[a] * g.c[b]; end
    Pf = polybound(f); Pg = polybound(g)
    rem = Pf * g.r[1] + Pg * f.r[1] + f.r[1] * g.r[1]
    @inbounds for k = nc+1:ft.nc2; rem += abs(buf[k]); end
    copyto!(h.c, 1, buf, 1, nc); h.r[1] = rem
    h
end
function tm_mul_full!(h::TM{Arb}, f::TM{Arb}, g::TM{Arb}, ws, buf::Vector{Arb} = [Arb(0; prec = PREC) for _ = 1:FTAB[].nc2])
    ft = FTAB[]; nc = TAB[].nc
    for x in buf; Arblib.zero!(x); end
    for (a, b, k) in ft.trip; Arblib.addmul!(buf[k], f.c[a], g.c[b]); end
    Pf = Mag(0); Pg = Mag(0); polybound!(Pf, f, ws); polybound!(Pg, g, ws)
    rem = Mag(0); tmp = Mag(0)
    Arblib.mul!(rem, Pf, g.r[1])
    Arblib.mul!(tmp, Pg, f.r[1]); Arblib.add!(rem, rem, tmp)
    Arblib.mul!(tmp, f.r[1], g.r[1]); Arblib.add!(rem, rem, tmp)
    for k = nc+1:ft.nc2; Arblib.get!(tmp, buf[k]); Arblib.add!(rem, rem, tmp); end
    for k = 1:nc; Arblib.set!(h.c[k], buf[k]); end
    Arblib.set!(h.r[1], rem)
    h
end
