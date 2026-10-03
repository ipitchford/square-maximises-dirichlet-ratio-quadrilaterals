
## 10. The affine-transport cover (TRANSPORT_SPEC.md) — the rest of the fundamental domain

### 10.1 Certificate
Reference quadrilateral Q0 = Q_{p0} (box centre), split along a diagonal into K1, K2. For p in the box,
Ψ_p is affine on each K with Ψ_p(K(p0)) = K(p): x = v_{a}(p) + A_K(p)(x̂ − v_{a}(p0)),
A_K(p) = M_K(p) M_K(p0)⁻¹ (affine in p, A_K(p0) = I); the two maps agree on the diagonal, so Ψ_p is a
continuous piecewise-affine homeomorphism Q0 → Q_p mapping ∂Q0 onto ∂Q_p whenever d_K = det A_K > 0 (certified
on the whole box: "no fold").
* **Test functions** u_k(p) = û_k(·;δ) ∘ Ψ_p⁻¹ with û_k(x̂;δ) = b_{Q0}(x̂)(P_k⁰(x̂) + Σ_i δ_i P_kⁱ(x̂)), δ = p − p0:
  continuous, piecewise polynomial, zero on ∂Q_p ⇒ u_k(p) ∈ H¹₀(Q_p). P_k⁰ = Rayleigh–Ritz eigenvectors at p0;
  P_kⁱ = central differences (step 1e-4) of the optimal eigenvectors of the transported N×N pencil at p0 ± h e_i
  (aligned: sign of u₁, Procrustes for (u₂,u₃)). The transported pencil at any p is assembled from full-basis
  reference matrices at p0 without new quadrature. Thin shapes make these derivatives noisy (|dC| ~ 10⁷), so the
  screen also tries P_kⁱ = 0 ("plain transport") and keeps the better; the choice is logged (`fo`).
* **Flux** (Lehmann–Goerisch): reference ŵ = d_K τ̂(δ) + curl ψ̂(δ) on K(p0), physical w = A_K ŵ/d_K (contravariant
  Piola). τ̂(δ) = φ̂(δ)·e, e = diagonal direction of Q0, e·∇φ̂(δ) = −û₁(·;δ), φ̂ = 0 on the other diagonal of Q0
  (a global polynomial; computed at the Duffy nodes by exact Lagrange interpolation along the Duffy lines, which are
  parallel to e); ψ̂(δ) = ψ̂⁰ + Σ δ_i ψ̂ⁱ (least-squares optimum of ‖w‖² at p0 and its p-derivatives). Then on each
  K: div_x w = (1/d_K) div_x̂ ŵ = div τ̂ = −û₁ = −u₁(p); the normal flux across the diagonal is (curl ψ̂·n̂)dŝ from
  both sides (Piola preserves normal flux, τ̂·n̂ = 0 on the diagonal), so w ∈ H(div, Q_p) and div w + u₁(p) = 0
  exactly for every p in the box.
* **Base quantities** (all p-dependence through A_K(p) and δ):
  K_ij = Σ_K (1/d_K) Σ_ab (adj A_K adj A_Kᵀ)_ab I^K_ab,ij(δ), M_ij = Σ_K d_K M^K_ij(δ),
  A2 = Σ_K [d_K⟨A_KᵀA_K, TT^K(δ)⟩ + 2⟨A_KᵀA_K, TC^K(δ)⟩ + ⟨A_KᵀA_K, CC^K(δ)⟩/d_K], where the reference integrals
  I^K_ab,ij(δ) = ∫_{K(p0)} ∂_a û_i(δ) ∂_b û_j(δ) etc. are quadratic polynomials in δ whose 15 + 15 + … coefficient
  integrals are computed **exactly** (Duffy–Gauss, m = n+6, integrand degree ≤ 2n+10 with deg ψ̂ = n+6) in Arb
  as ball-valued 15×15, 5×5, 5×10, 10×10 matrices (weighted `arb_mat_mul` products).
  Denominators are cleared with D = d_1 d_2 > 0 (homogeneity), giving polynomials in δ, represented as **order-4**
  Taylor models (exact full products; every coefficient of degree ≥ 5 folded into the remainder with |t^γ| ≤ 1).
* **Pointwise lower bounds** (bounds_tm.jl): instead of box constants (which must take the worst corner — for thin
  shapes λ₁ ≈ π²/w² varies by tens of per cent over a box), the criterion may use Taylor models β(p) ≤ λ₂(Q_p),
  L0(p) ≤ λ₁(Q_p) valid at every p: strip bounds π²|e_i|²/c_i(p)²; rectangles along fixed directions with affine
  side bounds (π²(1/Â²+1/B̂²), π²min(4/Â²+1/B̂², 1/Â²+4/B̂²) on a certified branch); Krahn–Szegő 2πj₀₁²/(1−a²−d²);
  Andrews–Clutterbuck L0(p) + 3π²/D̂(p)². Maxima over vertices/pairs are bounded by the dominant term plus a
  certified constant slack; reciprocals by the geometric series with remainder U^{K+1}/(1−U). The criterion is
  unchanged (β(p), L0(p) enter Ã as Taylor-model factors; homogeneity in (K, M, A2) is not affected).
* **Diagonal choice**: a vertex with angle close to π must be an end of the split diagonal (as the apex of a split
  triangle it makes A_K ill-conditioned over the box); the screen tries both diagonals and keeps the better.
* **Compound discard**: if no single exclusion holds on a box whose centre is excluded, the box is bisected
  internally (depth ≤ 6, widest side) and discarded if every leaf satisfies one single condition; the leaf tree is
  logged (e.g. `b(E2,a(E2,E3))`) and re-checked by the replay.
