# Addendum: affine-transport certificate for the remaining region (2026-10-01)

## Why
Fixed physical-coordinate test polynomials lose accuracy across a box when a near-flat vertex moves (your pilot: κ 3.4–3.9 at
radius 0.005). Transporting the box-centre test functions WITH the vertices removes that loss. Float prototype:
../explore/transport_proto.py — near-flat boxes (η = c_i/c_{i+2} ≈ 0.05–0.07) give bounds 1.63–1.91 at radius 0.02 (true ≤ 1.78);
generic boxes 1.88 at radius 0.05 (true 1.84). Predicted ~2×10⁴ boxes in total.

## The certificate on a box B = p0 + Π[−r_i, r_i]
Reference domain Q0 = Q_{p0}, split along v1v3 into K1 = (v1,v2,v3), K2 = (v1,v3,v4). For p ∈ B let A_K(p) = M_K(p) M_K(p0)⁻¹ with
M_K(p) = [v_b(p) − v_a(p), v_c(p) − v_a(p)] (2×2, affine in p, A_K(p0) = I), and Ψ_p affine on each K mapping K(p0) → K(p) (continuous).
- Test functions u_k(p) = û_k ∘ Ψ_p⁻¹, û_k = b_{Q0} P_k (first three RR eigenvectors at p0; optional first-order corrections in δ
  computed on Q0). u_k(p) ∈ H¹₀(Q_p) (continuous, piecewise polynomial).
- K_ij(p) = Σ_K (1/det A_K) Σ_{a,b} (adj A_K adj A_Kᵀ)_{ab} I^K_{ab,ij},  I^K_{ab,ij} = ∫_{K(p0)} ∂_a û_i ∂_b û_j.
  M_ij(p) = Σ_K det A_K · ∫_{K(p0)} û_i û_j.
- Flux (Lehmann–Goerisch): reference ŵ = det A_K · τ̂ + curl ψ̂ on K(p0), τ̂ = R(x) ê with R the antiderivative of û_1 along ê ∥ v3 − v1
  from the line v2v4 (div τ̂ = −û_1; τ̂·n = 0 on the diagonal v1v3), ψ̂ a fixed polynomial (least squares of ‖τ̂ + curl ψ̂‖² on Q0).
  Physical flux = contravariant Piola transform w = A_K ŵ / det A_K. Then div w = −u_1(p) EXACTLY for every p and w ∈ H(div)(Q_p)
  (normal flux continuous across the diagonal because τ̂·n = 0 there and curl ψ̂ is continuous).
  A2(p) = ‖w‖² = Σ_K [ det A_K ⟨G_K, TT^K⟩ + 2⟨G_K, TC^K⟩ + ⟨G_K, CC^K⟩ / det A_K ],  G_K = A_Kᵀ A_K,
  TT^K_{ab} = ∫τ̂_aτ̂_b, TC^K_{ab} = ∫τ̂_a(curl ψ̂)_b, CC^K_{ab} = ∫(curl ψ̂)_a(curl ψ̂)_b over K(p0).
- All reference integrals: Arb, exact Gauss on the two triangles of Q0 (integrand degrees: û ~ n+4, τ̂ ~ n+5, curl ψ̂ ~ kψ−1).
- L = (β A1 − A0)/(β A2 − A1), A0 = K11, A1 = M11, β a rigorous constant lower bound for λ₂ over B (your β machinery), β A1 > A0.
- Criterion exactly as before (cleared denominators: multiply by det A_1(p) det A_2(p) > 0 to make every entry polynomial in δ, or
  use Taylor models with 1/det series); bound over B with order-2/3 Taylor models or second-order centred forms. All p-dependence is
  through the two 2×2 matrices A_K(p), so the box evaluation is cheap; the cost is the per-centre setup.
- Required side conditions: det A_K > 0 on B (no fold), M ≻ 0, A11 < 0, β A1 > A0.

## Region
Cover F ∩ {convex} minus E1 (local cube, 0.049 is fine) — you may re-use your certified shell 0.049 < |p|∞ ≤ 0.225 and start the new
cover at |p|∞ ≥ 0.2 (overlap) — minus E2 (c_i ≤ 0.035 c_{i+2}) and E3 (width/diam ≤ 0.1). Use a COMPOUND discard test (a box is
discarded if every point satisfies at least one exclusion; implement by refining a box only where no single condition covers it).
Exact rational box endpoints throughout. Pilot first (box count, wall-clock prediction with units), then run if ≤ 12 h on 8 threads.
