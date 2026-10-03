# Darts: certified proof that λ₂/λ₁ < 5/2 for every non-convex (simple) quadrilateral — spec, 2026-10-01

## Setting
A dart Q has a reflex vertex D inside the triangle T = ABC of the other three vertices; boundary A→B→C→D→A (ccw),
notch N = conv(C,D,A), Q = T \ N = ABD ∪ BCD (two triangles sharing the segment BD).
Parameters (up to similarity): the angles of T (θ_A, θ_C; θ_B = π − θ_A − θ_C) and the barycentric coordinates
D = αA + βB + γC (α,β,γ > 0, α+β+γ = 1). β = notch depth δ/h (δ = dist(D, CA), h = dist(B, CA)).
Symmetry: the reflection exchanging A and C maps darts to darts (θ_A ↔ θ_C, α ↔ γ).

## Already known (do not redo)
- Prop 5.1 (paper/PAPER.md §5): β ≤ 0.0339 ⇒ ξ ≤ (7/3)/(1−β)² < 5/2 (homothety at B inside Q, Q ⊂ T, triangle theorem).
- FEM adversarial search (explore/dart_max.py): sup ξ over darts ≈ 2.3334 ≈ 7/3, approached only as β → 0 or as
  α → 0 / γ → 0 (thin spike) with a near-equilateral remaining triangle. Margin to 5/2 ≥ ~1/6 everywhere seen.
- Deep darts (β ≥ 0.5) have ξ ≤ ~2.0 (explore/dart_deep_angles.py; FEM unconverged for β ≥ 0.93 — use better meshes).
- Dart-bubble test functions u = ℓ_AB ℓ_BC max(ℓ_CD, ℓ_DA) P(x) (continuous, piecewise polynomial; kink on ℓ_CD = ℓ_DA)
  with Lehmann–Goerisch flux w = (−∫_{−∞}^x ũ(t,y)dt, 0) + curl ψ (ũ = u extended by 0; div w = −u exactly):
  explore/dart_lg_proto.py. Float results: essentially exact for β ≤ 0.1 and for the spike limit (α → 0); ratio bound
  < 2.45 for β ≤ ~0.3; fails for β ≥ 0.6–0.7 (LG loses 25–90% of λ₁ at the strong reflex singularity).
  At α = 0 or γ = 0 or β = 0 the dart bubble reduces exactly to a triangle bubble, so those faces are not singular.
- Liu's Crouzeix–Raviart bound λ_k ≥ λ_{k,h}/(1 + (0.1893 h)² λ_{k,h}) (X. Liu, Appl. Math. Comput. 267 (2015)), float
  prototype explore/cr_liu.py: robust at the reflex corner; ratio bound 1.79 at β = 0.9 with 8192 elements (true 1.46).
  The CR mass matrix is diagonal on triangles. A rigorous lower bound for the discrete λ_{1,h} needs a verified positive-
  definiteness test of K − σM (e.g. RCM-ordered banded Cholesky in floats + backward-error bound |E| ≤ γ_{b+1}|L||Lᵀ|,
  with the matrix entries themselves enclosed by interval assembly).
- Upper bounds for λ₂ are easy (Rayleigh–Ritz with any fixed H¹₀ functions evaluated with interval arithmetic).

## Tasks
1. Decide the region decomposition and the certificate per region, e.g.:
   (a) β ≤ 0.0339: Prop 5.1.  (b) moderate notches: Taylor models in the 4 parameters with the dart bubble + LG (the convex
   global step in ../global2/ is building the same machinery for convex quadrilaterals — reuse its code when available;
   coordinate through files only, do NOT modify ../global2/).  (c) deep notches: verified CR-FEM λ₁ lower bounds with
   inclusion boxes (Q_in ⊂ every dart of the box ⊂ Q_out; λ₂(Q) ≤ RR on Q_in, λ₁(Q) ≥ CR bound on Q_out), or a
   Taylor-model version of the discrete Temple bound.
2. Make the uncovered set compact by proving lemmas (with complete proofs, in darts/LEMMAS.md) for:
   (i) very deep notches β → 1 (the dart rescales to a V-shaped waveguide; consider Neumann cuts across the thin arms far
   from the junction — the arm tips then have μ ≥ π²/(cut length)² by the transverse Poincaré inequality — reducing to a
   compact family of "core" problems); (ii) thin hull triangles T (some angle of T → 0 or → π) — e.g. a 1-D
   Born–Oppenheimer lower bound λ₁ ≥ μ₁(−f'' + π²/h_max(x)²) along the long axis plus explicit upper bounds; (iii) check
   that the spike faces α → 0, γ → 0 and the shallow face β → 0 are covered by the box certificate (they are regular).
   If a limit cannot be closed, say so precisely and quantify the remaining set.
3. Pilot first (float): box counts and times per regime with units, BEFORE any long run.
4. CPU courtesy: while ../global2/ runs (check `pgrep -f julia`), use at most 2 threads at `nice -n 10`; run heavy jobs only
   after it finishes. One heavy job at a time (16 GB machine).
5. Deliverables in darts/: code, PILOT.md, LEMMAS.md, certificate logs, replay script, README.md stating exactly what is
   proved (rigorous), what is only numerically supported, and what remains open.
