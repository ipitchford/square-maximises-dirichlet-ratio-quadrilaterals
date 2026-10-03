# Darts: analytic lemmas and the rigour of the box certificate (2026-10-01)

Everything here is about the Dirichlet Laplacian. ξ(Ω) = λ₂(Ω)/λ₁(Ω). The target is ξ(Q) < 5/2 for every dart.

## 0. Setting, charts, symmetries

A **dart** is a simple quadrilateral with one reflex vertex. Label the reflex vertex D and the opposite vertex B, and write
the remaining two vertices A, C so that the boundary is A→B→C→D→A. The diagonal BD lies inside Q, and

  Q = int( ABD ∪ CBD ∪ BD ),   T := ABC (convex hull),   N := ACD (notch),   Q = int(T) ∖ N.

**BD chart.** By a similarity, B = (0,0), D = (1,0), A = (a_x, a_y) with a_y > 0 and C = (c_x, c_y) with c_y < 0. Every
quadrilateral of the form int(ABD ∪ CBD ∪ BD) with A, C on opposite sides of BD is simple; it is a dart with reflex vertex D iff

  f_D := (1 − c_x) a_y + c_y (a_x − 1) < 0       (cross product (D−C)×(A−D) < 0).

A dart whose reflex vertex is B is the same shape as one with reflex vertex D after relabelling, so the darts are exactly the
points (A, C) of the chart with f_D < 0. The reflection in the line BD gives the symmetry (a_x,a_y,c_x,c_y) ↦ (c_x,−c_y,a_x,−a_y).

Angles: b_A := ∠ABD, b_C := ∠CBD (so θ_B = ∠ABC = b_A + b_C), β := dist(D, AC)/dist(B, AC) (notch depth).
The shape space of darts is the open 4-simplex {(∠BAD, b_A, b_C, ∠BCD) > 0, sum < π}; its faces are: the shallow face
(reflex angle → π, D → AC, β → 0), the faces b_A → 0 or b_C → 0 (one triangle becomes a sliver along BD: collapse A → D or a
thin spike), and the faces ∠BAD → 0 or ∠BCD → 0 (an arm of infinite length).

## 1. Tools (all standard; proofs given or referenced)

**1.1 Domain monotonicity.** Ω ⊂ Ω' (open) ⇒ λ_k(Ω) ≥ λ_k(Ω').

**1.2 Mixed problems.** For a bounded Lipschitz domain Ω with ∂Ω = Γ_D ∪ Γ_N (Γ_N relatively open, a finite union of segments)
let V(Ω) be the H¹-closure of {u ∈ C¹(Ω̄): u = 0 near Γ_D} and μ_k the min–max values of ∫|∇u|²/∫u² on V.
*Neumann bracketing:* if Ω ∖ Σ = Ω₁ ∪ Ω₂ (disjoint open, Σ a finite union of segments) then λ₁(Ω) ≥ min μ₁(Ω_i), where on Ω_i the
cut Σ ∩ ∂Ω_i is Neumann and ∂Ω ∩ ∂Ω_i Dirichlet (restriction maps H¹₀(Ω) into V(Ω₁) ⊕ V(Ω₂)).
*Monotonicity for mixed problems:* if Ω ⊂ Ω', ∂Ω ∩ Ω' ⊂ Γ̄_D and Γ_N ⊂ Γ'_N, then zero extension maps V(Ω) isometrically into
V(Ω'), so μ₁(Ω) ≥ μ₁(Ω').

**1.3 Chord bound.** Let e be a unit vector. If every connected component of Ω ∩ (x + ℝe) is a segment of length ≤ ℓ whose two
endpoints lie in Γ̄_D, then μ₁(Ω) ≥ π²/ℓ². *Proof.* For u ∈ C¹ vanishing near Γ_D, on each component I the 1-D Dirichlet
Poincaré inequality gives ∫_I |∂_e u|² ≥ (π/|I|)² ∫_I u² ≥ (π/ℓ)² ∫_I u²; integrate over the transverse variable (Fubini) and use
density. ∎ (Pure Dirichlet case: Γ_N = ∅.)

**1.4 Parallelogram.** For P = {(x,y): 0 < x < L, sx < y < sx + h},  λ₂(P) ≤ π²(1+s²)/h² + 4π²/L².
*Proof.* u_k = sin(kπx/L) sin(π(y − sx)/h), k = 1, 2. With ξ = x/L, η = (y−sx)/h: ∂_x u = (kπ/L)cos(kπξ)sin(πη) − s(π/h)sin(kπξ)cos(πη),
∂_y u = (π/h) sin(kπξ)cos(πη). Every mixed product contains ∫₀¹ sin(πη)cos(πη)dη = 0 or ∫₀¹ cos(πξ)cos(2πξ)dξ = 0 or
∫₀¹ sin(πξ)sin(2πξ)dξ = 0, so u₁, u₂ are orthogonal in L² and in energy, with Rayleigh quotients π²(1+s²)/h² + k²π²/L². ∎

**1.5 Parallelogram in a convex region with an apex.** Let K = {(x,y): 0 < x < X, β(x) < y < τ(x)} with τ concave, β convex,
τ(0) = β(0) = 0, and let c(x) = τ(x) − β(x) attain its maximum ℓ at x* ∈ (0, X]. Then for every θ ∈ (0,1), K contains the
parallelogram {(1−θ)x* < x < x*, |y − sx| < (1−θ)ℓ/2} with s = (τ(x*) + β(x*))/(2x*). *Proof.* For 0 ≤ x ≤ x*, concavity and
τ(0) = 0 give τ(x) ≥ (x/x*)τ(x*); convexity and β(0) = 0 give β(x) ≤ (x/x*)β(x*). The vertical section of K at x therefore
contains the interval centred at sx of length (x/x*)ℓ ≥ (1−θ)ℓ. ∎

**1.6 Temple's inequality for a matrix pencil.** K symmetric, M symmetric positive definite, eigenvalues λ₁ ≤ λ₂ ≤ … of Kv = λMv.
For x ≠ 0 let ρ = xᵀKx/xᵀMx and ε² = (Kx−ρMx)ᵀM⁻¹(Kx−ρMx)/xᵀMx. If σ ≤ λ₂ and ρ < σ then λ₁ ≥ ρ − ε²/(σ − ρ).
*Proof.* Expand x = Σ c_k φ_k (M-orthonormal eigenvectors), b = xᵀMx. Then Σ(λ_k−ρ)(λ_k−σ)c_k² = ε²b (the cross term vanishes since
Σ(λ_k−ρ)c_k² = 0), and Σ(λ_k−σ)(λ_k−λ₁)c_k² ≥ 0 termwise (λ_k ≥ σ for k ≥ 2). Subtracting, ε²b + (ρ−λ₁)(ρ−σ)b ≥ 0. ∎

**1.7 Liu's lower bound for Crouzeix–Raviart eigenvalues** (X. Liu, Appl. Math. Comput. 267 (2015); interpolation constant from
X. Liu, F. Kikuchi, J. Math. Sci. Univ. Tokyo 17 (2010)). Let a mixed problem (1.2) on a polygon be discretised by
Crouzeix–Raviart elements on a conforming triangulation that resolves Γ_D (DOFs on Γ_D-edges set to 0, DOFs on Γ_N-edges free),
max element diameter h, discrete eigenvalues λ_{k,h} of (K, M). Then λ_k ≥ λ_{k,h}/(1 + C_h²λ_{k,h}), C_h = 0.1893 h.
*Proof (reproduced, so that the only external input is the constant).* Let Π be the CR interpolation (edge means). For u ∈ V,
Πu is admissible (edge means vanish on Γ_D), and on each element K, ∫_K ∇(u−Πu)·∇v_h = ∇v_h·∫_{∂K}(u−Πu)n = 0 for every
CR function v_h; hence a_h(u − Πu, v_h) = 0 (broken energy) and a(u,u) = a_h(Πu,Πu) + a_h(u−Πu,u−Πu). On every triangle,
u − Πu has zero mean on the three edges, so ‖u − Πu‖_{L²(K)} ≤ 0.1893 h_K ‖∇(u − Πu)‖_{L²(K)} (Liu–Kikuchi), and summing,
‖u − Πu‖ ≤ C_h ‖u − Πu‖_{a_h}. Let E_k be the span of the first k exact eigenfunctions. If C_h²λ_k ≥ 1 the claim is trivial
(x/(1+C²x) < 1/C²). Otherwise Π is injective on E_k (Πu = 0 would give ‖u‖ ≤ C_h‖u‖_a ≤ C_h√λ_k‖u‖ < ‖u‖), and for u ∈ E_k
with ‖u‖ = 1, A := ‖u‖_a² ≤ λ_k, t := ‖u−Πu‖_{a_h}: RQ(Πu) ≤ (A − t²)/(1 − C_h t)² ≤ A/(1 − C_h²A) (maximise over t; the maximum is
at t = C_hA). Min–max on Π(E_k): λ_{k,h} ≤ λ_k/(1 − C_h²λ_k), i.e. λ_k ≥ λ_{k,h}/(1 + C_h²λ_{k,h}). ∎
*Check of the constant (numerical, not part of the proof):* the sharp constant sup ‖w‖/(h‖∇w‖) over w ∈ H¹(K) with zero means on
all three edges, computed by P3 FEM with the three linear constraints (two refinement levels agree to 5 digits), is 0.18918 for the
equilateral triangle, 0.1885 for a needle, 0.1681 for the right isosceles, 0.162–0.175 for other shapes — consistent with 0.1893 as
the supremum over triangles.

**1.8 Rayleigh–Ritz.** For linearly independent u₁, u₂ ∈ H¹₀(Ω), λ₂(Ω) ≤ largest eigenvalue of the 2×2 pencil
((∇u_i,∇u_j), (u_i,u_j)). Continuous P2 finite element functions vanishing on ∂Ω are in H¹₀(Ω).

**1.9 Verified inertia.** If the float LDLᵀ factorisation (no pivoting, band b) of a symmetric float matrix A_fl has exactly one
negative pivot and no zero pivot, then L̂D̂L̂ᵀ = A_fl + ΔA with |ΔA| ≤ γ_{b+1}|L̂||D̂||L̂ᵀ| (Higham, ASNA, Thm 9.3 and the remark on
banded matrices; symmetric elimination is Doolittle with U = DLᵀ), and by Sylvester A_fl + ΔA has exactly one negative and no zero
eigenvalue. With A_ex = K − σ'M (exact) and ‖A_fl − A_ex‖₂ ≤ e_A, Weyl gives λ₂(A_ex) > −(e_A + ‖ΔA‖₂) (strict, since λ₂(A_fl + ΔA) > 0); since M is diagonal
positive (CR), λ₂(K − σM) ≥ λ₂(A_ex) + (σ'−σ)min_i M_ii > 0 for σ := σ' − (e_A + ‖ΔA‖₂)/min M_ii. Hence the pencil has at most one
eigenvalue ≤ σ, i.e. λ_{2,h} > σ. (‖·‖₂ ≤ max row sum for the symmetric bound matrices.) Implemented in `ldl.py`, `sigma.py`.

**1.10 Strain transfer between parameter points.** In the BD chart each triangle X ∈ {A, C} is the image of one fixed reference mesh
under F_X = [[1, x],[0, y]] (X = (x, y)); for CR (or P2) vectors v, a_p(v) = Σ_X ∫∇̂v̂ᵀ G_X(p)∇̂v̂, m_p(v) = Σ_X J_X(p) m̂_X(v), with
G_X = |det F_X| F_X⁻¹F_X⁻ᵀ (det G_X = 1) and J_X = |y|. Since vᵀG_X(p)v/vᵀG_X(q)v ≥ (η/η_q) λ_min(T Tᵀ), T = F_X(q)F_X(p)⁻¹
= [[1, (x_q − x)/η],[0, η_q/η]] (η = |y|), min–max gives
λ_{k,h}(p) ≥ λ_{k,h}(q) · min_X[(η_X/η_{X,q}) λ_min(T_XT_Xᵀ)] / max_X(J_X(p)/J_X(q)),  λ_min = 1 + (a − √(a² + 4b²))/2,
a = b² + c² − 1, b = (x_q − x)/η, c = η_q/η, evaluated in interval arithmetic over a box (monotone in a and b²).

## 2. Prop 5.1 (shallow notch) — from paper §5, restated
If β ≤ 0.0339 then ξ(Q) ≤ (7/3)(1 − β)⁻² < 5/2. (Homothety at B: T_s := B + (1−β)(T − B) ⊂ Q ⊂ T; triangle theorem.)
In the chart, β = |n·(D − A)|/|n·(B − A)| with n = (a_y − c_y, c_x − a_x) normal to AC; on a box this is enclosed by interval
arithmetic (`cover.beta_max`), and boxes with β ≤ 0.0339 throughout are discarded.

## 3. Lemma 3.1 (collapse collar: the reflex vertex close to a neighbouring vertex)
Let δ := dist(A, segment BD) and r := inradius of the triangle CBD. Then ξ(Q) ≤ (7/3)(1 + δ/r)², which is < 5/2 when δ ≤ 0.035 r.
(Symmetrically with A and C exchanged.)

*Proof.* dist(·, segment BD) is convex, so its maximum on the triangle ABD is attained at a vertex, i.e. it is δ; hence
ABD ⊂ N̄_δ(BD) ⊂ N̄_δ(CBD). With I the incentre of CBD, each side line of CBD is at distance r from I, so
N̄_δ(CBD) ⊂ I + (1 + δ/r)(CBD − I) =: H. Thus Q ⊂ int H and λ₁(Q) ≥ λ₁(CBD)/(1+δ/r)²; Q ⊃ CBD gives λ₂(Q) ≤ λ₂(CBD)
≤ (7/3)λ₁(CBD) (triangle theorem: Siudeja 2010; Arbon–Mannan–Psenka–Ragavan 2022). Finally (7/3)·1.035² = 2.4995… ∎

*Coverage.* In the BD chart, A = D + ρ(cos φ, sin φ): δ ≤ ρ. For C equilateral (|BD| = 1, r = 0.2887) this is ρ ≤ 0.0101.
This collar is a neighbourhood of the face b_A = 0 *only where A is close to D*; where the sliver ABD is long (A far along the ray
beyond D: a thin spike) the lemma does not apply (points of ABD near A are far from CBD).

## 4. Lemma 4.1 (BD-chord bound) — valid for every dart
λ₁(Q) ≥ π²/|BD|².

*Proof.* Write points as P = B + a e_A + c e_C (e_A, e_C unit vectors along BA, BC; a, c ≥ 0 on T) and D = B + d_A e_A + d_C e_C,
d_A, d_C > 0. (i) If P ∈ T with a ≥ d_A and c ≥ d_C then P ∈ N: indeed P − D ∈ cone(e_A, e_C), and e_A ∈ cone(A−D, C−D) because
e_A = μ(A−D) + ν(C−D) with ν = μ d_C/(|BC| − d_C) ≥ 0 and μ·|BA||BC|(1 − d_A/|BA| − d_C/|BC|)/(|BC| − d_C) = 1, μ > 0 (D is
interior to T); similarly e_C. So P ∈ (D + cone(A−D, C−D)) ∩ T = triangle ACD = N (the cone meets T exactly in ACD because D
and B lie on the same side of AC). Hence Q ⊂ W := {a > 0, c > 0, (a < d_A or c < d_C)}.
(ii) A line parallel to BD, parametrised as (a, c) = (a₀ + τd_A, c₀ + τd_C), meets W exactly for
τ_min < τ < τ_min + 1, τ_min := max(−a₀/d_A, −c₀/d_C): the conditions a, c > 0 give τ > τ_min and "a < d_A or c < d_C" gives
τ < 1 + τ_min. So every chord of Q parallel to BD lies in a segment of length |BD|; apply 1.3. ∎

## 5. Lemma 5.1 (thin wedge at B) — a complete analytic treatment of the corner θ_B → 0
If **both angles are acute**, b_A, b_C ∈ (0, π/2) (equivalently a_x > 0 and c_x > 0), and tan b_A + tan b_C ≤ 1/10, then
ξ(Q) ≤ 2.4223 < 5/2. In particular this holds if θ_B = b_A + b_C ≤ 0.0996 rad (5.7°): then both angles are acute and, tan being
superadditive on [0, π/2), tan b_A + tan b_C ≤ tan θ_B ≤ tan 0.0996 = 0.09993… < 1/10.
(The acute-angle hypothesis is essential: tan b_A + tan b_C ≤ 1/10 alone admits obtuse angles with negative tangents, e.g.
b_A = 10°, b_C = 160°, ρ_A = ρ_C = 10, a dart with c_x < 0 to which the proof below does not apply. The rigorous predicate
`regions.P_lemma51` checks a_x > 0, c_x > 0, a_y ≥ 0 ≥ c_y and the tangent sum on the whole box; the family runs use Lemma 5.3,
whose hypothesis θ_B ≤ 0.1772 implies acute angles.)
(The proof does not use the reflex condition; it holds for every quadrilateral int(ABD ∪ CBD ∪ BD) of the chart with these angles.)

*Proof.* BD chart; b_A, b_C < π/2, so A and C have positive abscissae and Q ⊂ {x > 0}. Put τ_B := tan b_A + tan b_C.
(i) *Stem.* For 0 < x < 1 the vertical section of ABD is [0, h_A(x)], where h_A is linear (if a_x ≥ 1) or a tent with apex at a_x
(if a_x < 1); in both cases h_A is concave on [0,1], h_A(0) = 0 and h_A(x) ≤ x tan b_A (ABD lies in the angle ∠DBA). Likewise CBD
gives [−h_C(x), 0]. Hence for 0 < x < 1 the section of Q is the single interval (−h_C(x), h_A(x)), of length c(x) := h_A + h_C ≤ xτ_B,
and the stem S := {0 < x < 1, −h_C(x) < y < h_A(x)} ⊂ Q is convex.
(ii) *Arms.* For x > 1 the section of ABD is empty if a_x ≤ 1, and otherwise is [y_DA(x), y_BA(x)], of length decreasing linearly
from tan b_A (= h_A(1)) at x = 1 to 0 at x = a_x; it lies above the x-axis. Likewise CBD below. The two pieces are disjoint (the open
ray beyond D is not in Q), so every vertical chord component for x ≥ 1 has length ≤ max(h_A(1), h_C(1)) ≤ c(1⁻).
(iii) Hence every vertical chord of Q has length ≤ ℓ := max_{(0,1]} c = c(x*), and λ₁(Q) ≥ π²/ℓ² by 1.3.
(iv) By 1.5 applied to S (τ = h_A, β = −h_C), with θ = 1/4, S contains a parallelogram with vertical sides (3/4)ℓ, length x*/4 and
shear s, |s| ≤ max(tan b_A, tan b_C)/2 ≤ τ_B/2 (because h_A(x*) ≤ x* tan b_A and h_C(x*) ≤ x* tan b_C). By 1.4,
λ₂(Q) ≤ π²(1 + s²)(4/(3ℓ))² + 4π²(4/x*)².
(v) With ℓ ≤ x*τ_B: ξ(Q) ≤ (16/9)(1 + τ_B²/4) + 64τ_B² ≤ (16/9)(1.0025) + 0.64 = 2.42222… ∎

## 5a. Lemma 5.3 (thin wedge at B, sector version) — extends Lemma 5.1 to θ_B ≤ 0.1772 rad (10.15°)
If θ_B = b_A + b_C ≤ 0.1772 (radians), then ξ(Q) < 5/2.

*Proof.* Use the notation of the proof of Lemma 5.1 (BD chart, stem S, maximiser x* of the vertical chord length c, ℓ = c(x*)).
Let φ_hi := arctan(h_A(x*)/x*) ≤ b_A and φ_lo := arctan(h_C(x*)/x*) ≤ b_C, α := φ_hi + φ_lo ≤ θ_B, T := ℓ/x* = tan φ_hi + tan φ_lo
≤ tan α (tan is superadditive on [0, π/2)). The triangle W* = {0 < x < x*, −x tan φ_lo < y < x tan φ_hi} lies in S ⊂ Q (proof of 5.1)
and contains the circular sector Σ(α, x*) = {0 < r < x*, −φ_lo < θ < φ_hi} (a point of the sector has abscissa r cos θ < x*).
(a) If α ≤ 0.0997 then T ≤ tan α ≤ tan 0.0997 = 0.1000316611… (note: slightly above 1/10) and the parallelogram argument of
Lemma 5.1 applies verbatim with τ_B replaced by T (it used only ℓ ≤ x*T and |s| ≤ T/2), giving
ξ ≤ (16/9)(1 + T²/4) + 64T² ≤ 2.4226304 < 5/2.
(b) If α ∈ [α_k, α_{k+1}] with α_k = 0.0997 + 0.0025k (exact decimals), k = 0..30 (31 intervals; α_31 = 0.1772), then
Σ(α, x*) contains a sector of angle α_k
with the same apex and radius, whose second Dirichlet eigenvalue is j²_{ν,2}/x*² with ν = π/α_k (eigenfunctions J_{mν}(j r/x*)
sin(mν(θ − θ₀)); the second eigenvalue is j²_{ν,2}/x*² because j_{ν,2} < 2ν < j_{2ν,1}). Hence λ₂(Q) ≤ j²_{ν,2}/x*², while
λ₁(Q) ≥ π²/ℓ² ≥ π²/(x* tan α_{k+1})² (chord bound, as in 5.1). So ξ(Q) ≤ (j_{ν,2} tan α_{k+1}/π)². For each k, `lemma_sector.py`
certifies in ball arithmetic (128 bits) an upper bound X_k > j_{ν,2} by the sign pattern J_ν(Y_k) < 0 < J_ν(X_k) with 0 < Y_k < X_k
checked explicitly (J_ν > 0 on (0, j_{ν,1})), checks X_k < 2ν, and checks (X_k tan α_{k+1}/π)² < 5/2; the largest value is
2.489242112… < 2.4893 (k = 30). Log: `lemma_sector.log`. ∎

## 5b. Lemma 5.2 (flat hull triangle) — the corner θ_B → π
Let θ_A = ∠BAC, θ_C = ∠BCA be the angles of the hull T at A and C. If **both are at most arctan(1/10) = 5.7106°** (so both are acute
and τ := max(tan θ_A, tan θ_C) ≤ 1/10), then ξ(Q) ≤ (16/9)(1 + τ²) + 64τ² ≤ 2.4356 < 5/2.
(A bound on the tangents alone would admit an obtuse angle with a negative tangent; the hypothesis is on the angles. The runs use
it through θ_B ≥ 174.29°: θ_B is the hull angle at B, so θ_A + θ_C = π − θ_B ≤ 5.71° < arctan(1/10), and both angles qualify.)

*Proof.* Place C = (0,0), A = (1,0), B above the x-axis. Both base angles are acute, so B projects into (0,1), and D ∈ int T
gives 0 < d_x < 1. With t(x) the upper boundary of T (tent C→B→A) and n(x) the upper boundary of N = CDA (tent C→D→A),
Q = {0 < x < 1, n(x) < y < t(x)}, so every vertical chord of Q is a single interval of length c(x) = t(x) − n(x) ≤ ℓ := max c,
and λ₁(Q) ≥ π²/ℓ² (1.3). Q_L := Q ∩ {x < d_x} = {0 < x < d_x, (d_y/d_x)x < y < t(x)} is convex with apex C (lower boundary
linear, upper concave, both 0 at x = 0); likewise Q_R := Q ∩ {x > d_x} with apex A. Let x* be a maximiser of c. If x* ≤ d_x apply
1.5 to Q_L (whose apex C is at the origin): Q contains a parallelogram with vertical sides (3/4)ℓ, length x*/4 and
shear |s| ≤ max(t(x*), n(x*))/x* ≤ tan θ_C (θ = 1/4); here ℓ ≤ t(x*) ≤ x* tan θ_C. If x* ≥ d_x apply 1.5 to Q_R, reflected by x ↦ 1 − x so that its apex A is at the origin, with x* replaced by
1 − x* and θ_C by θ_A. In both cases 1.4 gives λ₂(Q) ≤ π²(1+τ²)(4/(3ℓ))² + 4π²(4/L)² with ℓ/L ≤ τ, L the distance from x* to the
apex used. Divide. ∎

## 6. Lemma 6.1 (tails of long arms)
Let S be a segment cutting Q into Q₁ and Q₂, where Q₂ is a triangle with S as a side and its other two sides on ∂Q. Then the mixed
problem on Q₂ (Neumann on S) has μ₁ ≥ π²/|S|², hence λ₁(Q) ≥ min(μ₁(Q₁; Neumann on S), π²/|S|²).
*Proof.* Chords of Q₂ parallel to S end on the two other sides (Dirichlet) and are shorter than S (similar triangles); apply 1.3 and
1.2. ∎
Consequences used in the design of the (FEM) far-arm charts, superseded in the final cover by Lemma 11.3 (see §13): (a) for a cut S parallel to BD inside the arm
ABD, |S| ≤ |BD| (so the tail bound is at least the BD-chord bound); (b) the truncated core Q₁ of a dart with a long arm is contained,
with the same cut lines, in the truncated core of the limiting bent waveguide W of Lemma 4.1 (walls = rays from B and from D parallel
to BA, BC), so by 1.2 its mixed eigenvalue is bounded below by that of the waveguide core, a two-parameter compact family.

## 7. The faces of the shape simplex: what is regular for the box certificate
The box certificate (§8) uses the BD chart with each triangle meshed as an affine image of a fixed reference mesh.
* **Shallow face (β → 0, reflex angle → π):** regular — nothing degenerates in the chart, and the certificate is valid for convex
  quadrilaterals as well, so boxes may straddle the face. Prop 5.1 removes β ≤ 0.0339.
* **Faces b_A → 0 / b_C → 0 (one triangle degenerates to a sliver along BD):** NOT regular for the BD-chart certificate: G_X has
  entries ∝ 1/|y|, the reference elements in the sliver degenerate, and box sizes must shrink in proportion to a_y. The part of the
  face with the apex close to D (collapse) is covered by Lemma 3.1. The part with a long thin spike (A far along a ray close to BD;
  one of the critical curves, CBD equilateral, ξ → 7/3) is made regular by Lemmas 9.1, 9.2 and 9.4: the certificates there use the
  charts SPIKE (C B E D with a Neumann cut DE, for tan b_A ≤ t₁) and SPIKE2 (with the stub of Lemma 9.4, for t₁ ≤ tan b_A ≤ t₂),
  whose parameters are (t, C) only — uniform in |A|, including |A| → ∞. Implemented (`charts.spike_problem`,
  `charts.spike2_problem`) and run for C ∈ CBOX (Theorem S) and for C in the KF sector (certs2/kf_*); all other C are covered by
  the families SP_col/SP_spk/SF_col/SF_spk with the BO2 certificate (Lemma 11.3, spike arms as 'far' arms with R = 0) (§13).
* **Faces ∠BAD → 0 / ∠BCD → 0 (infinitely long arm), including the deep notch β → 1:** not compact in the BD chart; Lemma 9.3
  (+ Lemma 6.1) reduces them to compact families in (C, direction of BA, L). Charts implemented (`charts.far_problem`,
  `deep_problem`, `far2_problem`) and test boxes certified, but too weak for steep arms (PILOT §5); not used. The final cover handles
  this face with the BO2 certificate with tails (Lemma 11.3 (ii), families PF/FP/FF, ρ ≥ 20; §13). Lemma 9.5 closes the
  part of this face with b_X ∈ [40°, 140°] analytically (for |BX| ≥ L_k|BD|); outside [40°, 140°] Lemma 9.5 does not apply (its bound cannot reach 5/2 at all when b_X ≤ 39.23° or b_X ≥ 140.77°), and there the face is covered by the BO2 families alone.
* **Corner θ_B → 0:** Lemmas 5.1 and 5.3 (complete for θ_B ≤ 0.1772). **Corner θ_B → π (flat hull):** Lemma 5.2.

## 8. The box certificate (cert.py) — statement and rigour
For a box P = p₀ + r·[−1,1]⁴ in the chart (a_y > 0 > c_y on P), a reference mesh level k (κ = 2 grading toward D) and a number
σ₂, `BoxCert.certify` returns ok = True only if the following chain is verified in ball arithmetic (python-flint arb, 106 bits):
1. Exact reference data: all reference CR and P2 matrices are exact rationals rounded once to doubles (`refexact.py`);
   every float product used later carries a rigorous error bound (γ_n-type bounds + the rounding of the stored entries).
2. Fixed test vectors (any vectors are admissible): CR x(t) = x₀ + Σ t_i x_i + Σ t_it_j x_ij (eigenvector at p₀ and its scaled first
   and second parameter derivatives), P2 y_a(t) = y_a + Σ t_i y_{a,i}, a = 1, 2.
3. The p-dependence enters only through G_A, G_C (3 entries each) and J_A, J_C: G_X = (1/η)[[η²+x², −x],[−x, 1]], J_X = η.
   ρ(p) is computed in the cancelled form λ₀ + Σ_X η_X⁻¹[xᵀP_X x − λ₀η_X² xᵀM̂_X x]/b with exact polynomial brackets; all scalar
   functions are Taylor models in t of order 3 (`tm.ORDER`; order 2 in early probes only) with rigorous remainders (every
   truncated term is folded into the remainder).
4. The residual Kx − ρMx is scaled on each DOF set (A-interior by η_A, C-interior by η_C, BD by η_Aη_C) so that its coefficients are
   exact polynomials; its Taylor coefficient vectors are formed with explicit error vectors; ε² is bounded by
   (Σ_γ‖R_γ‖_w + ‖err‖_w)² per set (weights = inverse diagonal mass, rescaled by sup-bounds of η-powers over P).
5. Temple (1.6) with σ₂ ≤ λ_{2,h} on P from 1.9 at a centre + 1.10 over P; requires σ₂ > ρ on P.
6. Liu (1.7) with h = sup over P of the element diameters (`hmax_box`, monotone bound in x and η).
7. Rayleigh–Ritz (1.8) with the transported P2 vectors: N = K₂ − (5/2)L₁M₂ is shown negative definite on P via N₁₁ < 0 and
   N₂₂ + N₁₂²/|N₁₁| < 0 (Schur complement), with M₂ positive definite.
Then for every p ∈ P: λ₂(Q(p)) ≤ largest pencil eigenvalue < (5/2)L₁(p) ≤ (5/2)λ₁(Q(p)), i.e. ξ(Q(p)) < 5/2.

## 9. Reduction lemmas for the non-compact faces (spikes and long arms)

**Lemma 9.1 (spike head).** In the BD chart let 0 < b_A < π/2 and a_x ≥ 1 (A beyond the abscissa of D). Put t := tan b_A,
E := (1, t) (the point of BA above D), H := triangle BDE, P := int(CBD ∪ H ∪ BD) (a quadrilateral C, B, E, D). Then
  λ₁(Q) ≥ min( μ₁(P; Neumann on DE, Dirichlet elsewhere), π²/t² ),   λ₂(Q) ≤ λ₂(P)  (Dirichlet).
In particular ξ(Q) ≤ λ₂(P)/min(μ₁(P), π²/t²) depends only on (C, t), not on |A|.
*Proof.* DE cuts Q into P and the triangle DEA (if a_x > 1); the chords of DEA parallel to DE end on EA ⊂ BA and AD and have
length ≤ |DE| = t (Lemma 6.1), and P ∪ DEA ∪ DE = Q, so Neumann bracketing (1.2) gives the first claim. H ⊂ ABD because a_x ≥ 1
(H = {0 < x < 1, 0 < y < x t} and BA lies above it for x ≤ a_x), so P ⊂ Q and 1.1 gives the second. ∎
**Lemma 9.2 (monotonicity of the spike head in t).** For 0 < t ≤ t₁: μ₁(P(t)) ≥ μ₁(P(t₁)) and λ₂(P(t)) ≤ λ₂(CBD).
*Proof.* P(t) ⊂ P(t₁), ∂P(t) ∩ P(t₁) = the segment BE(t) ⊂ Dirichlet part, Neumann part DE(t) ⊂ DE(t₁): 1.2 (mixed monotonicity).
P(t) ⊃ CBD: 1.1. ∎
Remark (a_x < 1). If A = a(1, t) with 0 < a < 1 and t ≤ t₁, then A ∈ H̄(t₁) = tri(B, D, (1, t₁)), so ABD ⊂ H̄(t₁), Q ⊂ P(t₁) and
λ₁(Q) ≥ λ₁(P(t₁)) ≥ μ₁(P(t₁); Neumann on DE), while λ₂(Q) ≤ λ₂(CBD) still holds. Hence the collar certificate below covers
EVERY a > 0.
Consequence: the spike collar {a_x ≥ 1, tan b_A ≤ t₁} is covered as soon as λ₂(CBD) < (5/2)·min(μ₁(P(C, t₁)), π²/t₁²) — a
two-parameter (C) certificate; and {a_x ≥ 1, tan b_A ∈ [t₁, t₂]} needs a three-parameter (C, t) certificate. Neither depends on |A|.

**Lemma 9.3 (long arm, cut parallel to BD; design only — the final cover uses the oblique cut of Lemma 9.6).** Let
e_A = (A − B)/|A − B| and S_A := {B + s e_A + τ(D − B): s > 0, 0 < τ < 1} (the half-strip spanned by e_A and BD). Then ABD ⊂ S̄_A,
hence Q ⊂ int(CBD ∪ S̄_A). For L > 0 let Σ_L := {B + L e_A + τ(D − B): 0 ≤ τ ≤ 1} (a cut parallel to BD, |Σ_L| = |BD| = 1) and
suppose (i) 0 < L < |A − B| and (ii) the closed triangle CBD lies in the open half plane on B's side of the line through Σ_L
(i.e. all its points have s < L in the coordinates (s, τ) of the strip). Then the line of Σ_L crosses the segments BA and DA and not
BD, so Σ_L ∩ Q is a single segment cutting off the far part of ABD, and
  λ₁(Q) ≥ min( μ₁(core_L; Neumann on Σ_L), π² ),   core_L := int(CBD ∪ (S̄_A ∩ {s ≤ L})),
which depends on (C, direction of e_A, L) but not on |A|; and for fixed C and direction e_A, λ₂(Q) is non-increasing in |A|
(ABD(R) ⊂ ABD(R') for R < R' because A(R) lies on the segment B A(R')).
*Proof.* ABD = {B + λ(A−B) + μ(D−B): λ, μ ≥ 0, λ+μ ≤ 1} ⊂ S̄_A. The far part ABD ∩ {s > L} is a triangle whose chords parallel to Σ_L
end on BA and DA and are shorter than |BD| (Lemma 6.1); the near part of Q is contained in core_L with Q's piece of Σ_L inside the
Neumann segment Σ_L and the rest of ∂(Q ∩ {s < L}) inside core_L Dirichlet (1.2, mixed monotonicity). ∎
(Both arms long — the deep notch β → 1 — is the same construction applied to both triangles; the core then depends only on the two
directions, i.e. on (b_A, b_C), and L.)

These lemmas make every non-compact face of the shape simplex compact for the box certificate, except that the certificates
themselves (a chart with the pieces CBD, H or core_L, and Neumann DOFs on the cuts) have to be run. (Superseded for long arms by
Lemma 11.3; the spike-head charts are used by the KF runs and as the FEM fallback of the SP families.)

**Lemma 9.4 (long spike with a stub; used by `charts.spike2_problem`).** BD chart, t := tan b_A > 0, A = a_x(1, t) with a_x ≥ 1.
Put w := 3t, E := (1, t), D' := (1 + w, tw), E' := (1 + w, t(1 + w)) and P' := int(CBD ∪ BDE ∪ DD'E'E ∪ common edges), the
stub DD'E'E being the parallelogram between the line through D of slope t and BA. Then
  λ₁(Q) ≥ min( μ₁(P'; Neumann on D'E', Dirichlet elsewhere), π²/t² ),
and for every X₀ ∈ [1, a_x]: λ₂(Q) ≤ λ₂(Q_{X₀}), Q_{X₀} the dart with the same C and A replaced by X₀(1, t).
*Proof.* ABD = BDE ∪ (needle DEA). The line DA has slope t a_x/(a_x − 1) > t, so for 1 ≤ x ≤ min(a_x, 1 + w) the needle lies
above the line DD' and below BA, i.e. inside the stub. If a_x ≤ 1 + w then Q ⊂ P' with ∂Q ∩ P' ⊂ ∂Q (Dirichlet), and 1.2 gives
λ₁(Q) ≥ μ₁(P'). Otherwise cut the needle by σ = {1 + w} × [y_DA(1 + w), t(1 + w)] ⊂ D'E' (y_DA(1 + w) ≥ tw); the far part has vertical
chords ending on DA and BA of length ≤ t (Lemma 6.1), the near part lies in P' with its cut inside the Neumann segment D'E' (1.2).
For λ₂: X₀(1, t) lies on the segment from B to A, so ABD(X₀) ⊂ ABD and Q_{X₀} ⊂ Q (1.1).
*The endpoint a_x = 1.* Then A = E, ABD = BDE and Q = int(CBD ∪ BDE ∪ BD) ⊂ P'; the boundary of Q inside P' is the segment DE (a
common edge of BDE and of the stub), on which the functions of H¹₀(Q) vanish: extending them by zero gives functions in the form
domain of the mixed problem on P' (Dirichlet except on D'E'), so λ₁(Q) ≥ μ₁(P') ≥ min(μ₁(P'), π²/t²) (1.2); and X₀ = 1 is the only
admissible X₀, Q_{X₀} = Q. ∎
Note: in the BD chart every dart with c_x < 1 has a_x > 1 (f_D < 0 gives a_x − 1 > (1 − c_x)a_y/|c_y|).

**Lemma 9.6 (oblique cut of a long arm: the FARL chart; `charts.farl_problem`, validity `fam5.farl_setup`).** BD chart,
B = (0, 0), D = (1, 0). Let the arm apex be A = a (u, 1) with u = cot b_A ≥ 0 and a ≥ R > 0, and let C lie in the open lower half
plane. Fix rationals L > 0 and k₀ > 0 and put
  Q₁ := L (u, 1) (on the ray BA),   Q₂ := D + (L − k₀u)(u, 1) (on the parallel to BA through D),   ℓ := the line Q₁Q₂,
  P := int(CBD ∪ tri(B, D, Q₂) ∪ tri(B, Q₂, Q₁) ∪ common edges)   (the chart domain: pieces C, S1 = BDQ₂, S2 = BQ₂Q₁).
Assume (exactly the conditions verified by `farl_setup`, in ball arithmetic, for every u of the box):
  (H1) R > L > k₀u  (so 0 < L − k₀u, i.e. Q₂ lies beyond D on its ray, and Q₁ lies strictly between B and A on BA);
  (H2) B, D and the closed triangle CBD lie strictly on one side of ℓ (checked at B, D and the corners of the Cartesian box that
       contains C; the triangle is the convex hull of B, D, C);
  (H3) the point R(u, 1) lies strictly on the other side of ℓ.
Then  λ₁(Q) ≥ min( μ₁(P; Neumann on Q₁Q₂, Dirichlet elsewhere), π²/|Q₁Q₂|² ),  |Q₁Q₂|² = (1 − k₀u²)² + (k₀u)²,
and  λ₂(Q) ≤ λ₂(Q_R)  with Q_R the dart with A replaced by R(u, 1).
*Proof.* (a) The ray B + s(u, 1), s ≥ 0, meets ℓ only at Q₁ (s = L); by (H3) the points with s > L, in particular A (s = a ≥ R > L),
lie on the far side of ℓ, while B lies on the near side (H2). Hence ℓ crosses BA at Q₁ and, as D is on the near side and A on the
far side, crosses DA at one point X; it does not meet BD (H2). The closed strip between the line BA and its parallel through D
contains the segment DA, and ℓ meets this strip in the segment [Q₁, Q₂]; so X ∈ [Q₁, Q₂].
(b) The tip T := ABD ∩ {far side of ℓ} is the triangle A Q₁ X; its sides AQ₁ ⊂ BA and AX ⊂ DA lie on ∂Q (the triangle CBD is on the
near side by (H2), so it does not meet T), and its third side is σ := [Q₁, X]. By Lemma 6.1 the mixed problem on T (Neumann on σ)
has μ₁ ≥ π²/|σ|² ≥ π²/|Q₁Q₂|², since σ ⊂ [Q₁, Q₂]; and |Q₁Q₂|² = |(1 + (L − k₀u − L)u, −k₀u)|² = (1 − k₀u²)² + (k₀u)².
(c) The near part Q' := Q ∩ {near side of ℓ} = int(CBD ∪ polygon B D X Q₁ ∪ ...) is contained in P: the polygon BDXQ₁ lies in the
convex quadrilateral B D Q₂ Q₁ (X ∈ [Q₁, Q₂]; BQ₁ and DQ₂ are parallel, of positive lengths L√(1 + u²) and (L − k₀u)√(1 + u²)), which is the union of
the pieces S1 and S2. The boundary of Q' consists of σ ⊂ Q₁Q₂ (the Neumann edge of P) and of parts of ∂Q; the part of ∂Q' inside
P (the segment DX of DA) carries the Dirichlet condition. Extension by zero maps the form domain of the mixed problem on Q'
(functions vanishing on ∂Q' \ σ) into that of the mixed problem on P, so μ₁(Q'; Neumann on σ) ≥ μ₁(P; Neumann on Q₁Q₂) (1.2).
(d) Neumann bracketing along σ (1.2): λ₁(Q) ≥ min(μ₁(Q'; N on σ), μ₁(T; N on σ)), and (b), (c) give the claim.
(e) R(u, 1) lies on the segment BA (a ≥ R), so tri(B, D, R(u, 1)) ⊂ ABD, Q_R ⊂ Q and λ₂(Q) ≤ λ₂(Q_R) (1.1). ∎
Use. The FEM certificate of a FARL box bounds μ₁(P; Neumann on Q₁Q₂) from below over u in the u-interval and C in the Cartesian box
(CR chart: pieces C, S1, S2 with the Neumann cut), uses the tail constant π²/max|Q₁Q₂|² (`charts.cut_tail`: the maximum of the
convex quadratic g(s) = (1 − k₀s)² + k₀²s over the range of s = u², rounded upward, then π²/· rounded downward), and bounds
λ₂(Q_R) from above (P2 chart: the BD chart with A = R(u, 1)). The u-interval is [cot b₂, cot b₁] of the arm's angle interval
[b₁, b₂] ⊂ (0, π/2] (directed: the lower end rounded down, the upper end rounded up), R is a dyadic lower bound of
ρ_min sin b₁ (so a = ρ sin b ≥ R for every ρ ≥ ρ_min), and the chart box (u, C) is an exact-endpoint float box containing these
ranges (`nanguard.enclose`). For PP boxes with one long arm (ρ ≥ 4 on the box), ρ_min is that arm's smallest radius on the box
(the superset family with ρ ≥ ρ_min); for FP/PF, ρ_min = RF = 20; a long C arm is mirrored first (reflection in BD).

**Lemma 9.5 (very long steep arm; closes part of the long-arm face and of the deep notch β → 1).** BD chart, |BD| = 1,
L := |BA|, b := b_A ∈ (0, π). For every ℓ ∈ (0, L − 1):
  ξ(Q) ≤ (L + 1)² / (sin²b · (L − 1 − ℓ)²) + 4/ℓ².
*Proof.* λ₁(Q) ≥ π² (Lemma 4.1). Use coordinates s along e_A = (A − B)/L and τ along the unit normal pointing to D's side, so
B = (0,0), A = (L,0), D = (cos b, sin b). Put s₀ := max(cos b, 0) ≤ 1. For s₀ ≤ s ≤ L the section of the triangle ABD at abscissa s
is {0 ≤ τ ≤ τ_DA(s)} with τ_DA(s) = sin b (L − s)/(L − cos b): the edge BA is τ = 0, the edge BD lies in {s ≤ s₀}, and the third
edge is DA. τ_DA is decreasing, so the rectangle R := (s₀, s₀ + ℓ) × (0, w), w := sin b (L − s₀ − ℓ)/(L − cos b), lies in
ABD ⊂ Q, and w ≥ sin b (L − 1 − ℓ)/(L + 1). By 1.1, λ₂(Q) ≤ λ₂(R) ≤ π²/w² + 4π²/ℓ² (the (1,2) mode
sin(πτ/w) sin(2π(s − s₀)/ℓ) and the (1,1) mode span a 2-dimensional space). Divide. ∎
The right-hand side decreases in L for fixed ℓ and depends on b only through sin b. `lemma_longarm.py` (arb, 128 bits; log
`lemma_longarm.log`) certifies, for each row, a value ℓ_k with RHS(L_k, ℓ_k, b_k) < 5/2, hence:
**ξ(Q) < 5/2 whenever b_A ∈ [b_k, 180° − b_k] and |BA| ≥ L_k·|BD|** (and the same with A ↔ C), for
(b_k, L_k) = (40°, 1318), (41°, 418), (42°, 230), (43°, 155), (45°, 92), (47.5°, 61), (50°, 46), (55°, 31), (60°, 25), (65°, 21),
(70°, 19), (80°, 17), (90°, 16). The bound can never reach 5/2 when sin²b_A ≤ 2/5 (b_A ≤ 39.23° or b_A ≥ 140.77°).
*Deep notch.* In the hull placement D − B = (1 − β)(F − B) (F = foot of the cevian BD on AC), so |BD| = (1 − β)|BF|,
b_A = ∠ABF, b_C = ∠FBC. Hence, for fixed T and F, the limit β → 1 is closed by Lemma 9.5 as soon as one of ∠ABF, ∠FBC lies in
[40°, 140°] — for instance whenever 80° ≤ θ_B ≤ 140° — namely for 1 − β ≤ |BX|/(L_k|BF|) (X the steep arm). It is NOT closed when
both ∠ABF, ∠FBC are < 40° (e.g. T equilateral with F the midpoint: 30°/30°) or one of them exceeds 140°: there the chord bound
λ₁ ≥ π²/|BD|² is too weak (the true λ₁ is near the waveguide threshold π²/(|BD|² max(sin² b_A, sin² b_C))), and a numerical
certificate of the bent-waveguide core (Lemma 9.3; `charts.deep_problem`) would be needed.

## 10. The general certificate (gcert.py) — used for all runs
Charts are unions of triangular pieces q, each the affine image of an exact dyadic reference mesh (`chart.py`); the vertex
positions are polynomials in the box parameters (`charts.py`), so per piece G_q = N_q/d_q with N_q, d_q exact polynomials
(d_q = |det E_q|, sign checked on the box), J_q = d_q/|det R_q|. Steps 1–7 of §8 hold verbatim with: the CR chart may carry Neumann
edges (cuts; Liu's bound holds for mixed problems, 1.7) and its own geometry; the P2 chart (Dirichlet everywhere) may be a different
domain contained in every Q of the box (Lemmas 9.1, 9.4: CBD, or the dart Q_{X₀}); extra Neumann-bracketing constants c (tails)
enter as λ₁ ≥ min(L₁(p), c), so the Schur test is run for each of L₁ and the c's; the residual is scaled on each DOF group by the
product of the d_q of the pieces containing it, which makes all stiffness coefficients exact polynomials; σ₂ comes from ONE verified
LDLᵀ per box of the piecewise Loewner minorant Σ_q μ_q K_q(p_c) − s Σ_q j_q M_q(p_c) (`gsigma.sigma_box`): for every v,
vᵀK_q(p)v ≥ μ_q vᵀK_q(p_c)v with μ_q = min_box (J_q(p)/J_q(p_c)) λ_min(T Tᵀ), T = E_q(p_c)E_q(p)⁻¹ (LEMMAS 1.10), and
M_q(p) ≤ j_q M_q(p_c), j_q = max_box J_q(p)/J_q(p_c); hence K(p) − sM(p) ≥ A_low in the Loewner order and λ₂(K(p) − sM(p)) ≥ λ₂(A_low).

## 11. One-dimensional (Born–Oppenheimer) bounds for families of darts (`bo1d.py`, `bo2.py`, `p2fix.py`, `cover5.py`)

**Lemma 11.1 (1-D reduction with the longest slice component).** Let Ω ⊂ ℝ² be open (bounded or not), u a unit vector,
x = ⟨·, u⟩, and ℓ(x) the length of the LONGEST connected component of the slice Ω ∩ {⟨·,u⟩ = x} (ℓ(x) = 0 if the slice is
empty). Let V be measurable with 0 ≤ V(x) ≤ π²/ℓ(x)² wherever the slice is non-empty. Then for every φ ∈ C_c^∞(Ω)
  ∫_Ω |∇φ|² ≥ ∫ (v′² + V v²) dx,   v(x) := (∫ φ(x, y)² dy)^{1/2}   (y = ⟨·, u^⊥⟩),
and v ∈ H¹ with compact support in the projection of Ω. Consequently λ₁(Ω) ≥ inf σ(−d²/dx² + V) on that projection
(Dirichlet at finite ends).
*Proof.* Fix x. The slice is a union of disjoint open intervals J (its components); φ(x, ·) has compact support in the slice,
so on each J it vanishes near both ends and the 1-D Poincaré inequality gives ∫_J |∂_yφ|² ≥ (π/|J|)² ∫_J φ² ≥ (π/ℓ(x))² ∫_J φ²;
summing over J, ∫ |∂_yφ|² dy ≥ V(x) v(x)². For the x-derivative: v² is smooth, v is smooth where positive, and
d(v²)/dx = 2∫ φ φ_x dy (differentiation under the integral is legitimate because φ has compact support in Ω, so no
boundary terms arise), hence |v′| ≤ (∫ φ_x² dy)^{1/2} by Cauchy–Schwarz where v > 0, and v′ = 0 a.e. on {v = 0}. Integrate
∫|∇φ|² = ∫∫ (φ_x² + φ_y²) ≥ ∫ (v′² + V v²). Density of C_c^∞(Ω) in H¹₀(Ω) gives the eigenvalue statement. ∎

**Lemma 11.2 (geometry of an arm).** Let X be a point in the open upper half plane (lower for C), b = ∠XBD ∈ (0, π),
e(b) = (cos b, ± sin b) and T_X = tri(B, D, X), B = (0, 0), D = (1, 0).
(a) T_X ⊂ S̄(b) := {B + s e(b) + τ(D − B): s ≥ 0, 0 ≤ τ ≤ 1} (the half-strip spanned by e(b) and BD; its width is w(b) = sin b).
(b) A line L = {⟨·, u⟩ = x} (u a unit vector) not meeting the segment BD meets T_X in a segment (or not at all) whose
    endpoints lie on BX and DX; if ⟨e(b), u⟩ ≠ 0, its length is ≤ w(b)/|⟨e(b), u⟩|.
(c) The angle ∠XDB = atan2(ρ sin b, 1 − ρ cos b) (X = ρ e(b)) is strictly increasing in ρ.
(d) If X′ satisfies ∠X′BD ≤ ∠XBD and ∠X′DB ≤ ∠XDB (both on the same side of BD), then T_{X′} ⊂ T_X.
*Proof.* (a) T_X = conv{B, D, X} and X = B + ρ e(b), so every point is B + λρ e(b) + μ(D − B) with λ, μ ≥ 0, λ + μ ≤ 1:
s = λρ ≥ 0, τ = μ ∈ [0, 1]. (b) L ∩ T_X is convex, and its endpoints lie on ∂T_X = BD ∪ BX ∪ DX; L misses BD, so they lie on
BX and DX, hence in S̄(b) and in the full strip between the parallel lines B + ℝe(b) and D + ℝe(b), whose width (the distance
of D from the line B + ℝe(b)) is |e(b) × (D − B)| = sin b. The line L has direction u^⊥; a line with unit direction d crosses a
strip of width w and unit normal n in a segment of length w/|⟨d, n⟩|; with d = u^⊥ and n = e(b)^⊥ (rotation by 90° preserves
inner products) ⟨d, n⟩ = ⟨u, e(b)⟩. (c) d/dρ [ρ sin b/(1 − ρ cos b)] = sin b/(1 − ρ cos b)² > 0 on each branch, and the
angle is continuous through ρ cos b = 1 (where it equals π/2). (d) A point P above BD lies in T_X iff ∠PBD ≤ ∠XBD and
∠PDB ≤ ∠XDB; apply this to P = X′ and use convexity. ∎

**Lemma 11.3 (the BO2 box certificate).** Let 𝓕 be a family of darts given by arm specifications (`bo2.Arm`):
'polar' (X = ρ e(b), b ∈ [b₁, b₂], ρ ∈ [r₁, r₂]) or 'far' (b ∈ [b₁, b₂], ρ ≥ R ≥ 0), one for A and one for C. Fix a unit u.
(i) Compact part: the union over 𝓕 of T_X is contained in conv(P) with P = {B, D, r₁e(b₁), r₁e(b₂), r₂e(b₁), r₂e(b₂),
    (r₂/cos((b₂−b₁)/2)) e((b₁+b₂)/2)} (the annular sector lies in the convex hull of its four corners and of the intersection
    point of the tangents at the ends of its outer arc). (ii) Far arms: if ⟨e(b), u⟩ > 0 on [b₁, b₂], then for any x_cut ≥
    max(⟨B,u⟩, ⟨D,u⟩) the union over 𝓕 of T_X ∩ {⟨·,u⟩ ≤ x_cut} lies in conv{B, D, the four points where the rays
    B + ℝ₊e(b_i), D + ℝ₊e(b_i) (i = 1, 2) meet {⟨·,u⟩ = x_cut}} (by (a): S̄(b) ∩ {x ≤ x_cut} is the quadrilateral with these
    vertices for b = b_i, and the crossing points move monotonically along the cut line as b varies); and for x > x_cut the
    slice of each dart consists only of far-arm pieces, each of length ≤ c := max_{b∈[b₁,b₂]} sin b/⟨e(b), u⟩ by (b).
(iii) For a slab S_k = {x_k ≤ ⟨·,u⟩ ≤ x_{k+1}}, let ext_X(S_k) be the y-extent of conv(P_X) ∩ S_k, computed from the points of
    P_X inside S_k and from the intersections of ALL segments between points of P_X with the two slab lines (a superset of the
    vertices of conv(P_X) ∩ S_k), and, if S_k misses the projection of BD, also ≤ c_X by (b). If S_k meets the projection of BD
    let H_k = ext_A + ext_C, otherwise H_k = max(ext_A, ext_C) (slices not meeting BD split into an A-part and a C-part lying on
    different sides of the line BD, which are different components). Then every slice of every dart of 𝓕 with x ∈ S_k has all
    components of length ≤ H_k.
(iv) With V = π²/H_k² on S_k and V = π²/c² beyond a cut on a side towards which some far arm extends (V_∞ = π²/(max c)²),
    Lemma 11.1 gives λ₁(Q) ≥ μ for every Q ∈ 𝓕, where μ is the bottom of the spectrum of −d²/dx² + V. μ ≥ E is certified by
    Sturm shooting in ball arithmetic (`bo2._no_zero_arb`): the solution with v(x_L) = 0, v′(x_L) = 1 (or, at an unbounded left
    end, v = 1, v′ = √(V_L − E): the decaying exponential) stays positive across all cells (hyperbolic cells: positivity at
    both ends; oscillatory cells: phase condition), and at an unbounded right end v′ + √(V_R − E) v > 0 (then v stays positive
    for all larger x). A positive solution on the whole line/half-line implies E ≤ inf σ (Allegretto–Piepenbrink), since the
    operator is a limit of Dirichlet problems on bounded intervals for which Sturm comparison applies.
(v) Inner dart: by Lemma 11.2(c)–(d), X_in with ∠X_in BD ≤ b₁ and ∠X_in DB ≤ min_{𝓕} ∠XDB (the minimum over b at the smallest
    radius, enclosed by interval arithmetic over subintervals) satisfies T_{X_in} ⊂ T_X for all X of the spec; for spikes
    (A = a(1, t), a ≥ 1, t ≥ t₁) the point (1, t₁) has this property (Lemma 9.4). Q_in = int(T_{A_in} ∪ T_{C_in} ∪ BD) (or T_{C_in}
    alone) is contained in every Q ∈ 𝓕, so λ₂(Q) ≤ λ₂(Q_in) (domain monotonicity), and `p2fix.p2_cert` bounds λ₂(Q_in) by conforming
    P2 Rayleigh–Ritz with rigorous 2×2 forms (§8 step 7 with zero box radius).
**Conclusion:** if the Schur test certifies λ₂(Q_in) < (5/2)·E, then ξ(Q) < 5/2 for every Q ∈ 𝓕. ∎

Remarks. (1) The family may contain non-darts (convex quadrilaterals); the certificate holds for them too and is harmless.
(2) Far arms make the certificate uniform in the arm length: the far family {ρ ≥ R} is one box in ρ. (3) Validation:
`test_bo1d.py` checks the slab bounds of (iii) against brute-force slices of sampled darts (35 370 samples, no violation).
(4) Arm replacement: a polar arm spec (b ∈ [b₁, b₂], ρ ∈ [r₁, r₂]) may be replaced in (i)–(iv) by the far spec (b ∈ [b₁, b₂],
ρ ≥ r₁): the far family contains the polar one, so its lower bound E is a lower bound for the polar family as well (the inner
dart of (v) is still taken from the polar family). BO2 records with key `arms` ∈ {farA, farC, farAC} use this replacement for the
named arm(s) (`bo2.arm_variant`); the cover tries it for arms with r₂ ≥ 4, where the chord and tail bounds of the far arm are
sharper than the slab extents of the hull. (5) Any P2 mesh level k may be used in (v) (the record stores k; the cover tries
k = 4, 5, 6, since the P2 bound of long inner darts converges slowly in k); likewise any slab count K and slicing direction θ
(stored in the record).

**Lemma 11.4 (the INC box certificate: outer and inner darts).** Let 𝓑 be a polar × polar box of darts: b_A ∈ [b₁ᴬ, b₂ᴬ],
ρ_A ∈ [r₁ᴬ, r₂ᴬ], b_C ∈ [b₁ᶜ, b₂ᶜ], ρ_C ∈ [r₁ᶜ, r₂ᶜ], all angles in (0, π). For X = ρ e(b) write d(X) = ∠XDB =
atan2(ρ sin b, 1 − ρ cos b) (for C the mirror image; the formulas are the same).
(a) ∂d/∂ρ = sin b/|X − D|² > 0 and ∂d/∂b = ρ(cos b − ρ)/|X − D|². Hence over a polar arm box, max d = max_{b ∈ [b₁, b₂]} d(b, r₂),
    attained at b₁, at b₂ or at b* = arccos r₂ (only if r₂ < 1 and b* ∈ [b₁, b₂]), and min d = min(d(b₁, r₁), d(b₂, r₁)).
(b) Outer apex: if X⁺ lies in the open upper half plane with ∠X⁺BD ≥ b₂, ∠X⁺DB ≥ max d and ∠X⁺BD + ∠X⁺DB < π, then
    T_X ⊂ T_{X⁺} for every X of the arm box. Inner apex: if ∠X⁻BD ≤ b₁ and ∠X⁻DB ≤ min d, then T_{X⁻} ⊂ T_X for every X.
(c) Consequently Q⁻ := int(T_{A⁻} ∪ T_{C⁻} ∪ BD) ⊂ Q ⊂ Q⁺ := int(T_{A⁺} ∪ T_{C⁺} ∪ BD) for every Q of 𝓑, and
    λ₁(Q) ≥ λ₁(Q⁺), λ₂(Q) ≤ λ₂(Q⁻) (Dirichlet domain monotonicity).
(d) If L₁ ≤ λ₁(Q⁺) is certified by the general certificate (§10: CR + Temple + Liu, σ₂ from the verified LDLᵀ) on a Cartesian
    BD-chart box containing the point Q⁺ (in practice the box Q⁺ ± 2⁻³⁰), E is a float ≤ the exact lower end of the L₁ enclosure,
    and `p2fix.p2_cert` certifies λ₂(Q⁻) < (5/2)E (§8 step 7 at zero box radius), then ξ(Q) < 5/2 for every Q of 𝓑.
*Proof.* (a) With x = 1 − ρ cos b, y = ρ sin b, d(atan2(y, x)) = (x dy − y dx)/(x² + y²). For ρ: dy = sin b, dx = −cos b, so
x dy − y dx = (1 − ρ cos b) sin b + ρ sin b cos b = sin b > 0. For b: dy = ρ cos b, dx = ρ sin b, so x dy − y dx =
ρ cos b − ρ² = ρ(cos b − ρ). For fixed ρ the sign of ∂d/∂b is that of cos b − ρ: if ρ ≥ 1, d is decreasing on (0, π); if ρ < 1
it increases up to b* = arccos ρ and decreases after. So on [b₁, b₂] the maximum of d(·, ρ) is at b₁, b₂ or b*, and the minimum
(of a function increasing then decreasing) at b₁ or b₂. Since d increases in ρ, the box maximum is attained at ρ = r₂ and the box
minimum at ρ = r₁. (b) The closed triangle T_{X⁺} on BD is the set of points P of the closed upper half plane with
∠PBD ≤ ∠X⁺BD and ∠PDB ≤ ∠X⁺DB (intersection of its angular sectors at B and at D). Every X of the arm box has ∠XBD = b ≤ b₂ and
∠XDB ≤ max d, so X ∈ T_{X⁺}, and T_X = conv{B, D, X} ⊂ T_{X⁺} by convexity. The inner statement is Lemma 11.2(d).
(c) Unions and interiors are monotone; Dirichlet eigenvalues decrease under domain inclusion. (d) λ₂(Q) ≤ λ₂(Q⁻) < (5/2)E ≤
(5/2)L₁ ≤ (5/2)λ₁(Q⁺) ≤ (5/2)λ₁(Q). ∎
Implementation (`inc.py`): `d_range` evaluates (a) in ball arithmetic (if r₂'s enclosure straddles 1, the maximum is enclosed by
evaluating d over the whole interval [b₁, b₂]); the apexes are floats constructed from the target angles with a small outward
(outer) or inward (inner) pad and re-checked in ball arithmetic (`check_outer`, `check_inner`); apex radius ≤ 8. The replay
(`replay3.py`) re-checks the inclusion for the recorded apexes, recomputes L₁ at the recorded Q⁺ (must be ≥ the recorded E) and
re-runs the P2 test of Q⁻ with the recorded E. Use: family PP, where the 1-D bound of Lemma 11.3 is about 3× off for 'S-shaped'
short darts (A short and near D, C beyond D) and the box FEM certificate needs very small boxes (its strain factors μ_q penalise
the whole stiffness of a thin piece); the INC loss depends only on the spread of the four apex angles over the box.

## 12. A lemma for darts that are thin along BD

**Lemma 12.1 (thin stem).** BD chart (B = (0,0), D = (1,0)). For 0 ≤ x ≤ 1 let h_A(x), h_C(x) be the lengths of the vertical
slices of the triangles ABD, CBD at abscissa x, and c* := max_{0≤x≤1} (h_A(x) + h_C(x)). Then for every θ ∈ (0, 1)
  ξ(Q) ≤ (1 + 4c*²)/(1 − θ)² + 16 c*²/θ²;   in particular ξ(Q) ≤ 2.4356 < 5/2 if c* ≤ 0.05 (θ = 1/4).
*Proof.* (i) Chords. For 0 < x < 1 the vertical slice of Q is the single interval (−h_C(x), h_A(x)) (it crosses BD). For x > 1
the slice consists of the part in ABD (between DA and BA, non-empty only if a_x > 1, of length decreasing linearly from h_A(1)
to 0) and the part in CBD (likewise ≤ h_C(1)), which are different components (the ray beyond D is not in Q); so every component
has length ≤ max(h_A(1), h_C(1)) ≤ c*. For x < 0 likewise ≤ max(h_A(0), h_C(0)) ≤ c*. Every vertical chord ends on ∂Q, so
λ₁(Q) ≥ π²/c*² by 1.3.
(ii) Parallelogram. h_A is concave on [0, 1] (the upper boundary of ABD over [0, 1] is the lower envelope of the lines BA and DA,
or one of them), likewise h_C; let x* maximise c = h_A + h_C on [0, 1] and let e ∈ {0, 1} be the endpoint farther from x*
(|e − x*| ≥ 1/2). For x = λx* + (1 − λ)e, λ ∈ [1 − θ, 1], concavity gives h_A(x) ≥ λh_A(x*) + (1−λ)h_A(e) and the same for h_C,
so the slice of Q at x contains the interval centred at m(x) = [λ(h_A − h_C)(x*) + (1 − λ)(h_A − h_C)(e)]/2 (affine in x) of length
λc* + (1 − λ)c(e) ≥ (1 − θ)c*. Hence Q contains the parallelogram {x between x* and x* + θ(e − x*), |y − m(x)| < (1 − θ)c*/2} of
horizontal length L = θ|e − x*| ≥ θ/2, vertical side (1 − θ)c* and shear |m′| ≤ (|h_A − h_C|(x*) + |h_A − h_C|(e))/(2|e − x*|) ≤ 2c*.
By 1.4 and 1.1, λ₂(Q) ≤ π²(1 + 4c*²)/((1 − θ)c*)² + 16π²/θ². Divide by (i). For c* ≤ 0.05, θ = 1/4:
(1.01)/(0.5625) + 16·0.0025/0.0625 = 1.7956 + 0.64 = 2.4356. ∎
*Use.* Any box whose darts all have c* ≤ 0.05 is discarded ('Lemma 12.1'); c* is bounded on a box by the vertical-slicing slab
extents of Lemma 11.3(iii) over [0, 1] (rigorous upper bound of h_A + h_C). In particular every dart with tan b_A + ρ_C ≤ 0.05
(h_A ≤ x tan b_A, h_C ≤ |c_y| ≤ ρ_C), e.g. the corner {tan b_A ≤ 0.02, |C| ≤ e⁻⁴ = 0.0183} left by the spike/collar grids.

## 13. The global decomposition (every dart is covered)

Notation: BD chart, b_A = ∠ABD, b_C = ∠CBD ∈ (0, π), θ_B = b_A + b_C, ρ_A = |BA|, ρ_C = |BC| (|BD| = 1), t = tan b_A,
A = ρ_A e(b_A), C = ρ_C ē(b_C). β is the notch depth.

**Lemma 13.1 (dart condition in polar form).** The line AC meets the x-axis at abscissa f = sin θ_B/(sin b_A/ρ_C + sin b_C/ρ_A),
and Q is a dart with reflex vertex D iff θ_B < π and f > 1; then 1 − β = 1/f, i.e.
  1 − β = (sin b_C/ρ_A + sin b_A/ρ_C)/sin θ_B.
*Proof.* f = (a_x|c_y| + c_x a_y)/(a_y + |c_y|) (the zero of the affine interpolation of y along AC); insert the polar
coordinates: numerator ρ_Aρ_C sin(b_A + b_C), denominator ρ_A sin b_A + ρ_C sin b_C. D = (1, 0) is interior to the hull iff it
lies strictly between B and the point F = (f, 0) of AC on the ray BD, i.e. f > 1; and D = (1 − β)F by similar triangles. ∎
Consequences used as discards: (a) if sin b_C/ρ_A + sin b_A/ρ_C − 0.9661 sin θ_B ≥ 0 on a box, every dart of the box has
β ≤ 0.0339 (Prop 5.1) or is not a dart (for b_A, b_C ∈ [0, π] the left side is non-increasing in ρ_A and in ρ_C, so on a polar box
it suffices to check it at the largest radii, over a subdivision of the two angle intervals: `fam5.PP` predicate; for the
spike / collar families, where ρ_A is unbounded, the term sin b_C/ρ_A ≥ 0 is dropped and sin b_A/ρ_C − 0.9661 sin θ_B ≥ 0 is
checked at the largest ρ_C: `fam5.SP` predicate); (b) a dart has ρ_A > sin b_C and ρ_C > sin b_A (both terms are < sin θ_B ≤ 1).

**Theorem 13.2 (cover).** Every dart Q satisfies ξ(Q) < 5/2 provided the runs listed in README §2 are complete with FAIL 0
(each verified by its fail-closed replay). *Proof.* By the reflection in BD we may assume b_A ≤ b_C (the reflection swaps the
roles of A and C and preserves ξ); then b_A < π/2 (as b_A + b_C < π). Every such dart is in one of the following sets, each
covered by a lemma or by a run (runs cover their parameter grids exactly; their discards are re-checked predicates naming one of
these lemmas or one of the other sets):
 (1) θ_B ≤ 0.1772: Lemma 5.3.   (2) θ_B ≥ 174.29°: the hull angles at A and C sum to ≤ 5.71°, Lemma 5.2.
 (3) β ≤ 0.0339: Prop 5.1.      (4) max vertical stem chord ≤ 0.05: Lemma 12.1.
 (5) t ≤ t₁ := float(0.02) = 0.0200000000000000004163… ('collar', every ρ_A; covers t < 0.02): C polar with
     ρ_C ∈ [e⁻⁴, e³]: run SP_col (BO2 and FEM collar chart, Lemmas 9.1–9.2, 11.3)
     and, for C in the sector [38°, 100°] × [0.5, 1.8], the FEM run certs2/kf_collar; ρ_C ≥ 20: run SF_col; ρ_C < e⁻⁴:
     c* ≤ t + ρ_C < 0.02 + 0.0184 < 0.05, Lemma 12.1. (The polar and far grids overlap: e³ > 20.)
 (6) t₁ ≤ t ≤ 3/10 and a_x ≥ 1 ('spike', every ρ_A): ρ_C ∈ [e⁻⁴, e³]: run SP_spk (BO2 and FEM spike chart, Lemma 9.4) and
     certs2/kf_spike for C in the KF sector; ρ_C ≥ 20: run SF_spk; ρ_C < e⁻⁴: no dart (Lemma 13.1(b): ρ_C > sin b_A ≥ 0.0199).
 (7) all other darts with b_A ≤ b_C (t ≥ 0.3, or 0.02 ≤ t ≤ 0.3 with a_x < 1): ρ_A < e^{−2.5}: then sin b_C < ρ_A < 0.083 by
     13.1(b), so b_C < 4.8° (and θ_B ≤ 2b_C < 0.1772: set (1)) or b_C > 175.2° (set (2)); ρ_C < e⁻⁴: no dart (ρ_C > sin b_A ≥
     0.0199); ρ_A ∈ [e^{−2.5}, e³] and ρ_C ∈ [e⁻⁴, e³]: run PP (BO2, INC (Lemma 11.4) and FEM: BD chart, or for a box with exactly one long arm (ρ ≥ 4 on the box) the FARL chart (Lemma 9.6) of
     the far families applied to the superset family in which that arm is far with ρ ≥ its smallest radius on the box; the boxes inside the shallow box of
     Theorem H or in the region of Theorem M are discarded and covered by certs/shallow_crit, certs/main_cbox); ρ_A ≥ 20 and/or
     ρ_C ≥ 20 (the other polar): runs FP, PF, FF (BO2 with tails, Lemma 11.3 (ii)).
In every run the remaining discards are: θ_B ≥ π (not a dart), b_C < b_A on the box (the mirror image is in the domain), and the
sets (1)–(7) themselves. ∎
*Represented boxes (exactness of the grid; adversarial review 2026-10-02).* A box is stored as floats (p₀, r) and denotes the closed
box [p₀ − r, p₀ + r] in exact reals. Every cell (root or skipped) of the cover5 family manifests has exact float ends (the float
expressions p₀ ∓ r are exact; spike t-cells, whose nominal ends 0.02·1.5ⁱ are not dyadic, are encoded outward by
`run_fam.exact_cell`: the cell [0.02, 0.03] becomes [t₁, 0.030000000000000002], the last cell [0.2278125…, 0.3] becomes
[0.22781249999999997, 0.30000000000000004] ⊇ [·, 3/10], etc.). The processed roots of the legacy and KF runs (certs/, certs2/) are
stored with nearest-rounded encodings that are in general NOT exact-endpoint boxes (4 of 4 crit_collar, 52 of 52 crit_spike_X1,
63 of 75 main_cbox, 256 of 256 shallow_crit and 378 of 882 kf_spike roots are not). No predicate or certificate ever uses a
processed root's own encoding: every processed root, legacy or family, is represented by its outward dyadic snapping
(`grid.dyadic_root`, exact rational arithmetic), which its leaves tile exactly and which `accounting.py` checks to contain the
root's nominal cell. Every leaf box, every skipped cell and every chart box has exact float ends, checked by `nanguard.exact_box`
before a predicate or certificate forms the float expressions p₀ ∓ r, which are then the true ends. `accounting.py` checks, in rational arithmetic, that the represented boxes of each
family cover its parameter ranges (contiguous cells per axis; the col and spk t-cells together cover [0, 3/10]), that the KF runs
cover the t-ranges [0, t₁] (collar) and [1310/65536, 19661/65536] ⊇ [t₁, 0.30000000000000004] (spike, snapped slabs) used by the
SP predicate, and that the legacy manifests equal their generators. Decimal thresholds (0.0339, 0.05, 0.1772, 174.29°, 0.3, the
KF sector 38°, 100°, 0.5, 1.8, …) are exact decimals in every predicate (`nanguard.dec`, compared with directed ball endpoints).
