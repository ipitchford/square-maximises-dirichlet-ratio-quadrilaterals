# rules.jl — rigorous Gauss–Legendre rules on [0,1] (Arb balls) and the line-integration matrix S.
#
# Gauss–Legendre: nodes/weights from Arb's arb_hypgeom_legendre_p_ui_root (rigorous balls).
# Line-integration matrix: S[j,k] = ∫_0^{X_j} L_k(y) dy, with L_k the Lagrange basis polynomials on the
# m Gauss nodes X_1..X_m.  For every polynomial f of degree ≤ m−1:  ∫_0^{X_j} f = Σ_k S[j,k] f(X_k)  exactly.
# S is computed with the same m-point Gauss rule on [0, X_j] (exact: L_k has degree m−1 ≤ 2m−1), in ball
# arithmetic at 256 bits, then rounded to PREC-bit balls.  All quantities are continuous functions of the
# exact Gauss nodes, which lie in the input balls, so every ball encloses the exact value.

using Arblib
const PREC = 128
setprecision(Arb, PREC)

struct RuleArb
    m::Int
    X::Vector{Arb}
    W::Vector{Arb}
    S::Matrix{Arb}
end

struct RuleF
    m::Int
    X::Vector{Float64}
    W::Vector{Float64}
    S::Matrix{Float64}
end

function gauss01(m::Int, prec::Int)
    xs = Vector{Arb}(undef, m); ws = Vector{Arb}(undef, m)
    for k = 0:m-1
        r = Arb(prec = prec); w = Arb(prec = prec)
        Arblib.hypgeom_legendre_p_ui_root!(r, w, UInt(m), UInt(k))
        xs[k+1] = (r + 1) / 2
        ws[k+1] = w / 2
    end
    # sort increasing (roots are returned in decreasing order)
    perm = sortperm([Float64(Arblib.midref(x)) for x in xs])
    return xs[perm], ws[perm]
end

const _RULES_ARB = Dict{Int,RuleArb}()
const _RULES_F = Dict{Int,RuleF}()
const _RULE_LOCK = ReentrantLock()

function rule_arb(m::Int)
    lock(_RULE_LOCK) do
        get!(_RULES_ARB, m) do
            hp = 256
            X, W = gauss01(m, hp)
            S = Matrix{Arb}(undef, m, m)
            for j = 1:m, k = 1:m
                acc = Arb(0; prec = hp)
                for g = 1:m
                    y = X[j] * X[g]
                    Lk = Arb(1; prec = hp)
                    for i = 1:m
                        i == k && continue
                        Lk = Lk * (y - X[i]) / (X[k] - X[i])
                    end
                    acc += W[g] * Lk
                end
                S[j, k] = Arb(X[j] * acc; prec = PREC)
            end
            RuleArb(m, [Arb(x; prec = PREC) for x in X], [Arb(w; prec = PREC) for w in W], S)
        end
    end
end

function rule_f(m::Int)
    lock(_RULE_LOCK) do
        get!(_RULES_F, m) do
            R = rule_arb(m)
            mf(x) = Float64(Arblib.midref(x))
            RuleF(m, mf.(R.X), mf.(R.W), mf.(R.S))
        end
    end
end

# The two triangles of the split along the diagonal v1v3, both with P0 = v1, P2 = v3 (so the Duffy lines
# X = const are parallel to the diagonal e = v3 − v1 in both).  Triangle 1 (v1,v2,v3) is ccw (J = c_1 > 0
# for convex Q); triangle 2 (v1,v4,v3) is cw, J = −c_3, so its weight uses −J.
const TRIS = ((1, 2, 3, 1.0), (1, 4, 3, -1.0))
npoints(n::Int) = n + 6
