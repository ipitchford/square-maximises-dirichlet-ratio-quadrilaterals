# global2 — certified box cover for ξ = λ₂/λ₁ < 5/2 on convex quadrilaterals (global step)

## 0. Status (read this first)

Two certified pieces and one piece handed to a separate lemma. What has been replayed, exactly:
* **S (shell):** exhaustive replay of all 8 216 certificates, the E1 leaf and the coverage on 2026-10-01 05:59 with the
  replay.jl of that time (REPLAY OK, `runs/shell0225_replay.log`; that replay.jl did not yet reject unknown statuses —
  an independent scan in the external review found only recognised statuses and no duplicate paths). With the
  fail-closed replay.jl (17:47): strict loading, coverage and the E1 leaf of the whole log, plus 300 sampled
  certificates — PARTIAL REPLAY PASSED, 0 failures (`runs/shell0225_replay_failclosed_sample.log`). Log SHA-256
  `6e38aa632d23e4c6…`.
* **T (transport):** one clean, hash-bound, exhaustive replay of all 50 182 non-SHELL leaves of the frozen log
  (15:37–17:42, REPLAY OK, 0 failures; receipt `runs/RECEIPT_transport.md` and `runs/transport_replay_clean_RECEIPT.json`).

1. **Shell 0.049 < |p|∞ ≤ 0.225** (fixed-polynomial certificate, `cover.jl`): 8 216 certified boxes, 0 failures,
   smallest margin 2.36×10⁻⁵ (§7). With the local theorem and D4: λ₂/λ₁ < 5/2 for all convex Q_p with
   0 < |p|∞ ≤ 0.225.
2. **Rest of the fundamental domain** (affine-transport certificate, TRANSPORT_SPEC.md, `cover_t.jl`, §10): the
   a-priori box of F is tiled by 50 183 leaves: **31 426 certified**, 13 933 E2, 138 E4, 3 117 COMPOUND (every
   sub-leaf one of E2/E4/E3T), 1 568 E3T, and 1 SHELL leaf (the inner root [0,1/5]³×[−1/5,1/5], inside piece 1);
   0 FAIL/ERROR. Log structure verified (`check_log_t.py`); **clean, hash-bound replay of all 50 182 non-SHELL leaves: REPLAY OK** (receipt `runs/RECEIPT_transport.md`, `runs/transport_replay_clean_RECEIPT.json`; log SHA-256 `cb187ebb05553ae4…`). Smallest certified margin 5/2 − κ* = 5.05×10⁻⁶ (every certificate holds at κ = 5/2 with strict inequalities).
3. **The band 1/10 < width/diam ≤ 2/5 ("E3T") is handed to Lemma B′**, which is being proved separately in
   `../bandB/` (BAND_SPEC.md). The transport cover logs E3T leaves (and E3T sub-leaves of COMPOUND trees) without a
   certificate; the replay re-checks that width/diam ≤ 2/5 holds at every point of each of them. **Until Lemma B′ is
   proved the global step is incomplete** (§11).

Everything in a certified box's chain is ball arithmetic (Arb, 128 bit) except *choices* (test-function
coefficients and their δ-derivatives, stream function, split directions, degree n, diagonal, bound candidate,
rectangle directions), which are Float64 and are admissible whatever their values. Non-rigorous items: §8.

### Earlier finding (fixed-polynomial certificate, kept for the record)
The first design (fixed P_k in physical coordinates) was infeasible for the whole domain: pilot ≥ 1.2×10⁶ boxes,
realistic 10⁷–10⁸ (obtuse/near-flat vertices; §5–6). The transport certificate removes that obstacle.

## 1. Files

| file | content |
|---|---|
| `geom.jl` | Endo–Osting parametrisation (vertices, c_i, area, D4 action), Float64 width/diam/angles, point classification |
| `rules.jl` | rigorous Gauss–Legendre rules on [0,1] (Arb balls) and the line-integration matrix S |
| `fsolve.jl` | Float64 solver (NOT rigorous; only *chooses* test functions): Rayleigh–Ritz, LG flux with least-squares ψ, first-order p-derivatives (`fsolve_d`) |
| `tm.jl` | Taylor models of order 3 in 4 variables: Arb balls + Mag remainders (rigorous) or Float64 (pilot); norm-bounded products and exact "full" products |
| `assemble.jl` | Taylor models of the 13 base integrals over a box (exact quadrature, Taylor-centred evaluation of the fixed polynomials) |
| `criterion.jl` | LG / domain-monotonicity criterion, Haynsworth test, rigorous β and λ₁ scalar bounds, exclusion tests E1–E4 |
| `cover.jl` | adaptive bisection driver (threads, JSONL log, resumable) |
| `pilot.jl`, `pilot_summary.jl` | Float64 pilot of the whole domain (stratified sample of coarse cells) and its statistics |
| `calib.jl` | Arb timing and Float64→Arb agreement on pilot-passing boxes |
| `replay.jl` | independent replay of a log: coverage of the region, exclusions (separate implementation), certificates (recomputed) |
| `validate_points.jl`, `fem_compare.py`, `fem_refine.py` | Float64 validation against P2 FEM; D4 invariance |
| `PILOT_REPORT.md` | the pilot reports: calibration, box-count and time predictions (with units), §D for the transport cover |
| `transport.jl`, `transport1.jl` | affine-transport certificate: Float64 choice (incl. first-order δ-corrections), Arb reference integrals, box Taylor models |
| `bounds_tm.jl` | pointwise Taylor-model bounds β(p) ≤ λ₂, L0(p) ≤ λ₁ (strip, rectangles, Krahn–Szegő, Andrews–Clutterbuck, LG bootstrap) |
| `bound1d.jl` | rigorous 1-D reduction lower bound for λ₁ (box constant; used only in §11's analysis, not in the certificates) |
| `cover_t.jl`, `pilot_t.jl`, `pilot_t_summary.py`, `calib_t.jl`, `replay_t.jl` | transport cover driver, pilot, calibration, replay (fail-closed loading, results bound to record/log/verifier SHA-256, resumable, receipt; `rebuild_cert_t` rebuilds any logged certificate) |
| `margin_refine_t.jl`, `log_stats_t.py` | refined certified κ* of the tightest transport certificates (exact-rational bisection); statistics of a transport log |
| `strictjson.jl`, `logcheck.jl` | fail-closed log loading for both replays (strict JSON subset parser; status whitelist and per-status schemas; duplicates; SHA-256 of records, logs and verifier sources) |
| `MANIFEST.sha256` | SHA-256 of every source, specification, report and key run file in this directory (`shasum -a 256 -c MANIFEST.sha256`) |
| `runs/RECEIPT_transport.md`, `runs/transport_replay_clean_RECEIPT.json` | frozen-log record and the receipt of the clean, hash-bound transport replay |
| `check_log_t.py`, `xcheck_cf.jl` | structural check of a transport log (Python, independent of the Julia code: JSON, duplicates, tree, bitwise boxes, volumes, crash prefix); centred-form exclusion checks used by the replay |
| `tests_t.jl` | transport self-tests (point identity with the RR solver, Taylor-model enclosure vs direct assembly, D4, lower bounds below RR values) and semantic checks of a log (random points located in the tree) |
| `tests.jl` | consolidated self-tests (square, point box vs solver, TM vs exact, Arb ⊇ Float64, criterion, exclusions, D4, exact box containment) |
| `log_stats.py`, `margin_refine.jl` | statistics of a cover log; high-resolution certified κ* of the tightest boxes |
| `scratch/` | the exploratory experiments quoted in this README (flux variants, TM order, true-vs-certified loss, near-flat boxes); historical — some predate later API changes (e.g. the order-2 index tables) and no longer run unmodified |
| `runs/` | all logs (see §7) |

## 2. Mathematics

**Parametrisation** exactly as in GLOBAL_SPEC.md (vertices affine in p, c_i, |Q_p| = 1 − a² − d², D4 action,
fundamental domain F = {a, b, c ≥ 0}).

**Test functions.** u_k(x;p) = b(x;p)·P_k(x;p), k = 1,2,3, with bubble b = Π ℓ_k, ℓ_k(x) = (v_{k+1}−v_k)×(x−v_k), and

  P_k(x;p) = P_k⁰(x) + Σ_i (p_i − p0_i) P_kⁱ(x)                                    (**deviation 1, see below**)

with P_k⁰, P_kⁱ fixed polynomials (Legendre products in coordinates scaled by the bounding box of Q_{p0}, total
degree ≤ n). P_k⁰ = first three Rayleigh–Ritz eigenvectors at p0 (degree-n basis); P_kⁱ = central finite
differences (step 1e-4, same scaling, eigenvectors aligned to those at p0: sign of u₁, Procrustes rotation of
the (u₂, u₃) pair). For every p, u_k(·;p) is a polynomial vanishing on ∂Q_p, so u_k ∈ H²∩H¹₀(Q_p). Any choice
of these polynomials is admissible; the choice only affects tightness.
*Deviation 1 (allowed: the certificate holds for any p-dependent polynomial choice):* the spec fixes P_k; with
fixed P_k the true loss over a box is O(|δ|²) and the pilot box count was 5.6×10⁷ on a Monte Carlo sample of R
even with perfect bounding (`scratch/t_truebox.jl`). The first-order terms make the defect O(|δ|²) and the
quadratic-in-defect losses O(|δ|⁴).

**Upper side (Rayleigh–Ritz).** K = (∇u_i,∇u_j), M = (u_i,u_j), 3×3.

**Lower side (Lehmann–Goerisch, n = 1).** Flux w = ∇u₁/ρ0 + τ + curl ψ with ρ0 = θ₁(p0) (a fixed number),
τ = φ·e, e = v3 − v1, and

  φ(x) = −(1/ρ0) ∫_{σ(x)}^{0} (Δu₁ + ρ0 u₁)(x + s e) ds,  x + σ(x)e ∈ line(v2, v4)       (**deviation 2**)

i.e. the antiderivative of the residual along the diagonal direction, starting on the *other* diagonal. φ is a
global polynomial (σ is affine in x), e·∇φ = −r/ρ0, so div w = Δu₁/ρ0 − r/ρ0 = −u₁ exactly for every p.
ψ(x;p) = ψ⁰ + Σ_i (p_i − p0_i)ψⁱ is a polynomial of degree ≤ n+6 (Tikhonov-regularised least squares
minimiser of ‖w‖² at p0, and its finite-difference p-derivatives). *Deviation 2:* the spec's x-direction τ
needs an inner Gauss rule per quadrature node (≈ 8× the node count); the diagonal τ is computed from values
already available at the nodes (below). With deg ψ = n+6 every global polynomial particular solution of
div τ = −r/ρ0 of degree ≤ n+5 spans the *same* affine flux space, so the LG value equals the spec's one; the
Float64 test (`scratch/t_tau.jl`) reproduces the prototype and improves it (157° test angle, n=14: LG 0.056 %
below FEM λ₁ vs 0.093 % for lg_proto.py). A piecewise version (start lines on the edges) was rejected: kinked,
1–10 % worse.

**Exact identity for A2.** Since u₁ ∈ H¹₀, div τ = −r/ρ0 and div curl ψ = 0:
  A2 = ‖w‖² = (2ρ0 A1 − A0)/ρ0² + ‖τ + curl ψ‖²   (A0 = K11, A1 = M11).
The code assembles T = ‖τ + curl ψ‖² and adds (2ρ0 M11 − K11)/ρ0² (checked to 1e-11 against the direct ‖w‖²).
This keeps the large part ∇u₁/ρ0 out of the remainder chain.

**LG theorem used** (proof in a few lines): for u ∈ H¹₀, w ∈ H(div) with div w + u = 0, and 0 < β ≤ λ₂ with
βA1 > A0:  λ₁ ≥ L := (βA1 − A0)/(βA2 − A1). (Expand u = Σ a_k φ_k; A2 ≥ Σ a_k²/λ_k by the minimal-norm property of
∇(−Δ)⁻¹u; Σ a_k²(λ_k−β)(λ_k−λ₁)/λ_k ≥ 0 gives λ₁(βA2* − A1) ≥ βA1 − A0, hence the claim, and βA2 − A1 > 0.)

**β (a constant ≤ λ₂ on the whole box), all rigorous:**
max of (i) Krahn–Szegő 2πj₀₁²/max_B|Q_p|; (ii) enclosing rectangles π²·min(4/A²+1/B², 1/A²+4/B²) for FIXED
directions u (axis, the four edges and the two diagonals of Q_{p0}; the extents along u, u⊥ are maxima of affine
functions of p, so their sup over the box is exact); (iii) the Andrews–Clutterbuck gap λ₂ ≥ λ₁ + 3π²/D² (convex
domains) with λ₁ ≥ L0 := max(π²(1/A²+1/B²) [rectangles], π²/(4r²) [Hersch, convex, r = inradius upper bound from
strips and certified-bounded drop-triangles]); (iv) the same gap with λ₁ ≥ L_LG(β₀) (bootstrap, recomputed in Arb).
(iii)/(iv) and the generalised rectangle directions extend the spec's β (they were needed: β ≈ λ₁ for elongated
shapes, where LG is impossible).

**Two modes of the criterion** (κ = 5/2 for the theorem; smaller κ gives the logged certified bound ξ < κ):
* LG mode:  Ã = K·(βA2 − A1) − κ(βA1 − A0)·M  (= (βA2 − A1)(K − κ L M), βA2 − A1 > 0);
* MONO mode (thin/elongated shapes, where β ≤ λ₁): Ã = K − κ L0·M with the scalar L0 ≤ λ₁ above.
If Ã₁₁ < 0 and S′ = −Ã₁₁Ã_cc + Ã_c1Ã_1c (c = {2,3}) has a negative eigenvalue (any of tr S′ < 0, det S′ < 0,
S′₁₁ < 0, S′₂₂ < 0), Ã has ≥ 2 negative eigenvalues (Haynsworth inertia additivity), so K − κLM is negative on a
2-dimensional coefficient subspace; there the u_z are non-zero, hence λ₂ ≤ max R(u_z) < κL ≤ κλ₁. Side conditions
certified: βA1 − A0 > 0 and βA2 − A1 > 0 (LG mode), M ≻ 0 (M11 > 0, leading minor > 0, det M > 0).

**Homogeneous normalisation.** All criterion polynomials are homogeneous in (K, M, A2) (Ã of degree 2, S′ of
degree 4, det S′ of degree 8). Before evaluating them, all 13 base Taylor models are multiplied by one polynomial
q(t) ≈ 1/M11(t) with exactly representable coefficients, and q > 0 is certified on the box. This changes no sign
(Ã → q²Ã, S′ → q⁴S′) but removes the common-mode dilation variation of the integrals (b ∝ L⁸ under dilation,
M11 varied ±21 % over a 0.02-box), which otherwise destroys the bounds of the degree-8 criterion.

## 3. Rigour

* **Quadrature exactness.** Q split along v1v3 into (v1,v2,v3) [Jacobian c₁] and (v1,v4,v3) [−Jacobian = c₃];
  m×m Duffy–Gauss–Legendre product rule, m = n + 6. A polynomial of total degree d on a triangle becomes, after
  the Duffy map, degree ≤ d+1 in X and ≤ d in Y, so the rule is exact for d ≤ 2m − 2 = 2n + 10. Degrees: M: 2n+8,
  K: 2n+6, T = |τ + curl ψ|²: 2(n+5) (φ has degree n+5, curl ψ degree n+5). Nodes/weights are Arb balls from
  `arb_hypgeom_legendre_p_ui_root`. At every p with c₁, c₃ > 0 (in particular every convex p) the quadrature
  sum *equals* the integral; the Taylor models enclose the sum for every p in the box.
* **Line integration (τ).** On a Duffy line (X = X_i fixed) the residual r is a polynomial of degree ≤ n+4 in Y,
  so ∫_0^{Y_j} r = Σ_k S_jk r(Y_k) and ∫_0^{Y*} r = Σ_k I_k(Y*) r(Y_k) exactly (Lagrange basis on the m = n+6
  nodes; S_jk, and the Taylor coefficients of I_k at Y*(p0), are Arb balls computed at 256 bits with the same Gauss
  rule, exact for the degree-(m−1) integrands). The line v2v4 is the line Y = Y* in both triangles (Y* = the
  fraction of v1v3 at which the diagonals cross; Y* = cross(v4−v2, v2−v1)/cross(v4−v2, v3−v1) is a rational
  function of p, enclosed by a Taylor-model reciprocal with remainder U⁴/(1−U)).
* **Taylor models** (order 3 in t = δ/r ∈ [−1,1]⁴, Arb coefficients, Mag remainders): vertices are exact linear
  models; node positions are exact linear models (fixed barycentric weights); the fixed polynomials are evaluated
  at moving nodes by their finite Taylor expansion at the thin node centre with the derivative table computed in
  ball arithmetic (tail Σ_{a+b>3} |D_ab| H_u^a H_v^b/(a!b!), H = Σ|h_i|), so no dependency blow-up; every product
  bounds the dropped terms of degree > 3 by Σ_{d1+d2>3} N_{d1}N_{d2} (N_d = sum of |coefficients| of degree d,
  using |t^γ| ≤ 1); per-box products (normalisation and criterion) are computed in full (degree ≤ 6) and the
  degree-4..6 coefficients are folded into the remainder with their exact values. No truncation without folding;
  no order reduction of an existing remainder (both pitfalls flagged by the local-step review are avoided).
* **Range bounds** of a Taylor model over the unit box: exact 1-D ranges of c_i t_i + c_ii t_i², |c_ij| for cross
  terms, |c_γ| for degree 3, plus the remainder (ball arithmetic throughout).
* **Box geometry is exact.** Boxes are [lo, hi] with Float64 endpoints (exact binary rationals); children share
  the Float64 midpoint; the Taylor-model box p0 ± r is verified in exact rational arithmetic to contain [lo, hi].
  Thresholds are exact rationals: cube = 49/1000 (read as the decimal written in config.toml), η = 35/1000,
  s₀ = 1/10; the cube split point is Float64(49/1000, RoundDown).
* **Exclusions** (whole box, one fixed condition, Arb Taylor models): E1 exact rational comparison; E4 c_i ≤ 0;
  E2 c_i − η c_{i+2} ≤ 0; E3 cross(e_i, v_j − v_i)² ≤ s₀²|e_i|²|v_k − v_l|² for both j ∉ {i, i+1} and a fixed pair
  (k,l) (then width ≤ max_j dist(v_j, line_i) ≤ s₀|v_k − v_l| ≤ s₀·diam).
* **Scalar bounds** (β, L0, D, r) are Arb computations rounded outward; sums of Float64 bounds are re-done in
  Arb with downward rounding.

## 4. Validation (Float64 pipeline, before any rigorous run)

* **Square (p = 0):** θ₁ = 2π² to 8e-13 and θ₂ = θ₃ = 5π² to 1.4e-10 (n = 10); LG(β = 5π²) − 2π² = −4e-14 (n = 8), 4e-15 (n = 10).
* **LG flux vs the prototype** (`scratch/t_tau*.jl`; triangle + bump test shapes of `../explore/lg_proto.py`, FEM λ₁
  in brackets): δ = 0.1 (157° angle) n = 10/12/14: L = 44.749 / 44.777 / 44.791 (44.8164), rel. gap 0.15 / 0.089 /
  0.056 %; prototype 44.679 / 44.753 / 44.774. δ = 0.03 (173°): 49.607 / 49.749 / 49.831 (49.9826).
* **FEM comparison** (`validate_points.jl` → `validation_points.csv`, `fem_compare.py` → `validation_fem.log`):
  25 points of R (7 near the square with |p|∞ ∈ [0.04, 0.1], 6 near-triangles with min c_i/c_{i+2} < 0.06,
  6 with width/diam ∈ [0.1, 0.2], 6 random), n = 8 and 12, P2 FEM k = 6: **all 50 rows satisfy
  L ≤ FEM λ₁ ≤ θ₁**. θ₂ ≥ FEM λ₂ failed in 10 rows by ≤ 8.4e-7 relative — all smooth shapes where the
  degree-12 Rayleigh–Ritz value is *more* accurate than P2 FEM at k = 6; with k = 8 (`fem_refine.py` →
  `validation_fem_refine.log`) FEM λ₂ drops below θ₂ at all 13 re-run points. For the 12 rows where β < θ₁
  (thin shapes; LG not applicable) the domain-monotonicity bound gives θ₂/L0 ∈ [1.32, 2.31].
* **D4 invariance** (`validate_points.jl`, printed to `runs/validation_d4.log`): θ₁..θ₃ and A2 agree over the 8
  images of 5 points to 1e-15 … 1e-7 (relative; the larger values for thin/near-triangle shapes come from the
  conditioning of the basis).
* **Taylor models vs exact values:** at 10–30 random points of boxes of radius 0.01–0.02, |TM(t) − exact| is
  0.2–5 % of the remainder bound (`scratch/t_asm.jl`, `scratch/t_order.jl`); point boxes reproduce the solver
  (M = I, K = diag θ, A2 to 1e-15).
* **Arb vs Float64:** coefficients agree to within ball radius + 1e-12 relative, remainders to 1e-9 relative
  (`scratch/t_time.jl`); 20/20 pilot-passing boxes are certified in Arb (`runs/calib.txt`).
* **Certified κ vs the true maximum over a box** (fixed test functions, 150 sampled points incl. corners): equal to
  4 digits where the shape is moderate (1.8408 vs 1.8405); the gap opens only near obtuse corners (§6).

## 5. Pilot of the whole domain (Float64, mimics the rigorous test exactly up to rounding)

Method (`pilot.jl`, `pilot_summary.jl`): the a-priori box (a ≤ 0.905, b, c ≤ 4.36, |d| ≤ 0.905, split at the
cube faces) is refined with the rigorous exclusion tests alone down to the attempt level (all widths ≤ 0.3):
**867 coarse cells** survive (total volume 1.66; 937 E4, 137 E2, 125 E3, 1 E1 cells discarded). 25 cells chosen
uniformly at random were run through the full adaptive Float64 pilot with a cap of 2500 processed boxes per cell.
The Float64 test is the same Taylor-model computation as the rigorous one (only rounding differs), so no extra
remainder model is needed; the calibration (20 boxes, `runs/calib.txt`) confirms 20/20 Float64 passes certify in
Arb, at 0.52 s Arb + 0.09 s Float64 per certified box (n = 8, one thread).

| sampled cells | finished | capped (≥ 2500 boxes) | of which: no certificate at all |
|---|---|---|---|
| 25 | 13 (mean 337 processed, 143 certified-candidate leaves) | 12 (48 %) | 4 |

**Prediction.** Counting capped cells *at the cap* gives a lower bound of 1.19×10⁶ processed boxes and 1.03×10⁵
certificates; at the measured throughput (60 007 processed boxes in 977 s on 8 threads, `runs/shell025_pilot.log`)
this is ≥ 5.4 h of Float64 screening + ≥ 1.9 h of Arb, *before* the capped cells are finished. The capped cells
do not finish: (a) in the near-flat-vertex layer the true loss with polynomial×bubble test functions exceeds the
margin already at box radius 0.005 (κ_true = 3.4–3.9 at r = 0.005 for three shapes with 174–176° angles whose point
κ is 1.85–2.13; `scratch/t_e2box.jl`); this layer (min c_i/c_{i+2} < 0.08) is 16 % of R (volume 0.049), i.e.
≈ 8×10⁷ boxes of width 0.005; (b) boxes straddling the junction of two different exclusion conditions (e.g. E2 for
two indices, or E2/E4) cannot be discarded by one fixed condition and are refined to the minimum width
(4 of 25 sampled cells: hundreds of exclusion leaves, zero certificates). **Predicted total: ≥ 10⁷ boxes, ≫ 12 h on
8 threads → the full run was not started** (as the spec requires).

## 6. Why the full cover is infeasible with this certificate, and what would make it feasible

1. **Obtuse and near-flat vertices (the dominant cost).** Any smooth function vanishing on two edge lines through a
   vertex has a critical point there, so b·P behaves like r² at every vertex, while the eigenfunctions behave like
   r^{π/α}. For α > 90° the polynomial P must approximate φ/b ~ r^{π/α−2}, which is singular (r^{−0.78} at 148°,
   ~r^{−1} at 176°). At a fixed p this only costs accuracy (point κ stays good: 1.85–2.13 at 174–176°, LG within
   0.1–0.8 % of FEM at n = 10–14). Over a box the vertex moves, the sharp feature of the fixed polynomial does not,
   and the defect grows like |δ|·r^{π/α−3} near the corner. Measured: κ_true (fixed test functions + first-order
   corrections, sampled) exceeds 5/2 at box radius 0.005 for 174–176° shapes, and the ‖τ + curl ψ‖² remainder
   (cancellation between τ and curl ψ) blows up in the corner directions already at 148°. First-order
   p-dependence of P and ψ (implemented) helps for moderate shapes (κ_true 1.96 → 1.88 at r = 0.02) but cannot
   follow a moving singularity. What would fix it:
   * corner-adapted test functions (add r^{π/α} sin(πθ/α) at the obtuse vertex, with α = α(p) handled by Taylor
     models of the power function), or transported test functions u(x;p) = b·P(G_p x) with G_p the
     piecewise-affine map fixing the vertices (then the singular feature moves with the corner; LG needs a flux
     absorbing the normal-derivative jump of u₁ across the diagonal);
   * for the near-flat layer (min c_i/c_{i+2} ∈ [0.035, ~0.1]), a sharper near-triangle comparison: for convex Q,
     the homothety of T = conv(v_i, v_{i+2}, v_{i+3}) from v_{i+3} with ratio 1 + c_i/c_{i+2} contains Q (the
     barycentric coordinates of v_{i+1} with respect to (v_i, v_{i+2}, v_{i+3}) are (c_{i+1}, c_{i+3}, −c_i)/c_{i+2}, and
     c_{i+1} + c_{i+3} = c_i + c_{i+2} = 2|Q|), hence
     ξ(Q) ≤ ξ(T)·(1 + c_i/c_{i+2})². With the *certified ξ(T) of the actual triangle* instead of 7/3, this covers
     every near-flat shape whose triangle is not near-equilateral (polynomial bounds on triangles have no
     near-flat vertex). The residual set (near-equilateral T + small bump) is small and has margin ≈ 6 %.
2. **Exclusion junctions.** A box whose points satisfy different exclusion conditions (E2 for two indices, E2 and
   E4, …) cannot be discarded by one fixed polynomial inequality; the cover refines such boxes down to the minimum
   width (5 of 25 sampled cells). Fix: a compound exclusion test "every point satisfies at least one of the
   conditions", certified by an internal subdivision of the box with a fixed index per sub-box (cheap: no
   eigenvalue work), logged with its sub-box list so that the replay can re-check it.
3. **Constant factors.** Arb (128-bit balls) costs 0.52 s per certified box at n = 8 with order-3 Taylor models
   (Float64: 0.02 s). All balls stay tiny (radii ~1e-35), so midpoint–radius Float64 arithmetic with rigorous
   rounding-error sweeping (Rump-style) would be ~10× faster. The Float64 screening (choice of n and of the
   first-order test functions: ~12 small eigen-solves per attempt) dominates on boxes that are later split.

## 7. Rigorous result: the shell 0.049 < |p|∞ ≤ 0.225

Command: `julia -t 8 cover.jl runs/shell0225_cert.jsonl --region 0,0.225,0,0.225,0,0.225,-0.225,0.225`
(roots = the region box split at the cube faces ±Float64(49/1000, RoundDown); the inner root is discarded by E1).
Log: `runs/shell0225_cert.jsonl` (JSONL, one line per certified/discarded box and per split), progress
`runs/shell0225_cert.log`, statistics `python log_stats.py runs/shell0225_cert.jsonl`.

| quantity | value |
|---|---|
| processed boxes | 16 410 (8 193 splits) |
| certified | **8 216** (LG mode 8 209, domain-monotonicity mode 7; degree n = 8 for all) |
| discarded | 1 (E1, the local cube); no E2/E3/E4 boxes in this region |
| FAIL / ERROR | **0 / 0** |
| negativity test that certified | zᵀS′z 4 723, det S′ 3 415, tr S′ 78 |
| wall time | 1 171 s on 8 threads (Float64 screening + Arb); 0.99 s per certified box on one thread incl. its κ* search |
| certified volume | 5.1143×10⁻³ = volume of the region minus the cube (coverage cross-check) |
| certified bound κ* on ξ per box | median 2.342, min 1.886; 459 boxes with κ* ≥ 2.49 |
| **smallest certified margin** | **5/2 − κ* = 2.36×10⁻⁵** (box `R10.d1.d1.b1.c0.a0.d0.d1.a0`, refined by 40-step bisection in Arb: `runs/margin_refine.txt`) |

Every certified box satisfies the criterion at κ = 5/2 with strict inequalities; κ* is the smallest κ for which the
same box certificate still holds (logged with a 14-step bisection; the 15 tightest boxes refined to 40 steps).
Consequence (with D4 invariance and the local theorem for 0 < |p|∞ ≤ 1/20): **λ₂/λ₁ < 5/2 for every convex
quadrilateral Q_p with 0 < |p|∞ ≤ 0.225.** (Monte Carlo: every p with |p|∞ ≤ 0.25 is convex and in R.)
The Float64 pilot of the 0.25-cube (`runs/shell025_pilot.jsonl`, stopped at 60 007 processed boxes) already
covered 100 % of |p|∞ ≤ 0.225 and 99.88 % of |p|∞ ≤ 0.25; the uncovered remainder sits at the outer corners
(a ≈ |d| ≈ 0.25, angles up to 159°), i.e. the obtuse-corner obstacle of §6.

**Replay** (`julia -t 8 replay.jl runs/shell0225_cert.jsonl --region 0,0.225,0,0.225,0,0.225,-0.225,0.225`,
output `runs/shell0225_replay.log`, 1 076 s on 8 threads): **REPLAY OK** — the tree rebuilt from the logged splits has
8 217 leaves, 0 uncovered nodes, 0 bound mismatches, 0 unused records; the E1 leaf re-checked by the independent
implementation; all 8 216 certificates recomputed from the log alone (Float64 choice at the logged degree, rigorous
β/L0, Arb Taylor models, Arb criterion; logged β accepted only if ≤ the rigorously recomputed admissible shift):
0 failures. Self-tests: `julia tests.jl` → 480/480 pass (`runs/tests_output.log`).
That replay used the replay.jl of 05:59, which did not yet reject unknown statuses. The fail-closed replay.jl
(strict JSON, status whitelist and schemas, duplicates, unreached records; §8) was run on the same log at 17:47 with
`--certs-sample 300`: log accepted, coverage complete, E1 re-checked, 300/300 sampled certificates re-verified —
PARTIAL REPLAY PASSED (`runs/shell0225_replay_failclosed_sample.log`). An exhaustive fail-closed re-run is
`nice -n 10 julia -t 3 replay.jl runs/shell0225_cert.jsonl --region 0,0.225,0,0.225,0,0.225,-0.225,0.225` (≈ 1 h).

## 8. What is not rigorous / caveats (loud)

* **The thin band 1/10 < width/diam ≤ 2/5 is NOT covered here** (E3T leaves of the transport cover and E3T sub-leaves
  of its COMPOUND trees; §11). Every such box is re-checked by the replay to satisfy width/diam ≤ 2/5 at every point,
  so the gap is contained in that band. It is handed to Lemma B′ (`../bandB/`, BAND_SPEC.md), which is not yet
  proved; until it is, the global step is incomplete.
* Choices made in Float64 are not rigorous and need not be — the certificate is valid for any values: test-function
  coefficients and their first-order δ-derivatives, ψ̂ and its derivatives, the degree n, the split diagonal, the
  lower-bound candidate, rectangle directions, the normalising polynomial q, split directions, compound split points.
* Classical results used without re-proof (stated from the literature; no web search in this session):
  Krahn–Szegő (λ₂ ≥ 2πj₀₁²/|Ω|), domain monotonicity (strips, rectangles), rectangle eigenvalues, Hersch's
  inequality λ₁ ≥ π²/(4ρ²) (convex), the Andrews–Clutterbuck fundamental-gap theorem λ₂ − λ₁ ≥ 3π²/D² (convex),
  Haynsworth inertia additivity, the Lehmann–Goerisch bound for n = 1 (proof sketch §2; its β-monotonicity
  follows from A₁² ≤ A₀A₂), and for the transport certificate the elementary facts that a continuous piecewise-
  affine homeomorphism preserves H¹₀ and that the contravariant Piola map preserves H(div) and normal fluxes (§10).
* j₀₁ is enclosed by the decimal bounds given in the spec (not recomputed).
* Independence of the replays: certificates are recomputed from the logs with the *same* kernel code (assemble /
  transport / criterion / bounds_tm); coverage (tree rebuild, exact child boxes, root boxes) and the exclusion tests
  are separate implementations (in the transport replay: centred forms of the exact polynomials with adaptive
  bisection, `xcheck_cf.jl`, then a uniform ball subdivision; the cover's own Taylor-model test only as a counted
  "TM-only" fallback, §10.2). The replays are deterministic recomputations, not independent implementations of the
  certificate; independent exact-rational re-certification of sampled leaves is done separately (`../verify/`).
  Both replays fail closed (`logcheck.jl`, `strictjson.jl`): every line is parsed by a strict JSON parser, statuses are
  whitelisted with per-status schemas (exact key sets, types, ranges, compound-tree grammar), duplicate paths,
  unknown/FAIL/ERROR/pending records, uncovered nodes and unreached records give a nonzero exit, and the unqualified
  "REPLAY OK" is printed only by exhaustive runs. Transport replay results are bound to the SHA-256 of the exact log
  line of their leaf, of the log and of the verifier (all sources, config, Julia/Arblib versions, precision); cached
  results are reused only when all hashes and the kind match. Software correctness of Julia/Arblib/FLINT is assumed.
* The logged `kappa_cert` of a transport certificate is a ladder value (first accepted rung above the Float64
  estimate, 2.5 if none) and a Float64 rendering of the rational actually tested (κ = 5/2 itself is always tested
  first, exactly); the quoted margins come from the replay's exact-rational bisection (§10.2).
* Deviations from GLOBAL_SPEC.md and TRANSPORT_SPEC.md (all admissible, documented in §2 and §10): first-order
  p-dependent test polynomials and stream function; diagonal-direction τ; extra lower-bound candidates (rectangles,
  Taylor-model bounds, Andrews–Clutterbuck bootstraps); a domain-monotonicity (MONO) mode; Taylor models of order 3
  (shell) / 4 (transport) with exact full products; homogeneous normalisation by q ≈ 1/M11; diagonal choice; and the
  deferred band E3T (a deviation in scope, not in rigour).

## 9. Reproduction

```
cd global2
julia tests.jl                                                    # self-tests of the shell code (~2 min)
julia validate_points.jl > runs/validation_d4.log                 # Float64 validation points + D4 check
~/.venvs/analysis-py313/bin/python fem_compare.py                 # FEM comparison (validation_fem.log)
~/.venvs/analysis-py313/bin/python fem_refine.py                  # FEM k=6,7,8 near the square
julia -t 8 pilot.jl 32 2500 runs/pilot1                           # whole-domain pilot, fixed polynomials (stopped after 25 cells)
julia pilot_summary.jl runs/pilot1 867 2500 0.6 0.95 1.4          # its statistics
julia -t 8 cover.jl runs/shell025_pilot.jsonl --pilot --maxboxes 60000 --region 0,0.25,0,0.25,0,0.25,-0.25,0.25
mkdir runs/calib_src && cp runs/shell025_pilot.jsonl runs/pilot1/cell_{5,17,25}.jsonl runs/calib_src/
julia -t 1 calib.jl runs/calib_src 60 runs/calib.txt              # Arb calibration (20 boxes)
julia -t 8 cover.jl runs/shell0225_cert.jsonl --region 0,0.225,0,0.225,0,0.225,-0.225,0.225   # rigorous shell run (resumable)
julia -t 8 replay.jl runs/shell0225_cert.jsonl --region 0,0.225,0,0.225,0,0.225,-0.225,0.225  # independent replay
julia -t 1 margin_refine.jl runs/shell0225_cert.jsonl 15          # smallest certified margins (shell)
python log_stats.py runs/shell0225_cert.jsonl
# transport cover (§10–§12)
julia -t 8 pilot_t.jl 48 3000 runs/pilot_t6                       # Float64 pilot on 48 random coarse cells
python pilot_t_summary.py runs/pilot_t6 1737 3000                 # extrapolation (PILOT_REPORT §D)
julia -t 1 calib_t.jl runs/pilot_t6 20 runs/calib_t.txt           # Arb calibration of transport certificates
julia -t 8 cover_t.jl runs/transport_cert.jsonl > runs/transport_cert.log 2>&1   # rigorous run (resumable)
python check_log_t.py runs/transport_cert.jsonl runs/transport_cert.jsonl.bak_prereboot   # log structure (+ crash prefix)
nice -n 10 julia --heap-size-hint=1500M -t 3 replay_t.jl runs/transport_cert.jsonl --results runs/transport_replay_clean.jsonl >> runs/transport_replay_clean.log
                                                                  # clean hash-bound replay (resumable: re-run the same command); receipt runs/transport_replay_clean_RECEIPT.json
julia -t 8 margin_refine_t.jl runs/transport_cert.jsonl 200       # refined κ* of the tightest certificates
julia -t 1 tests_t.jl runs/transport_cert.jsonl > runs/tests_t_output.log       # transport self-tests + log semantics
python log_stats_t.py runs/transport_cert.jsonl
python scratch/band_volume.py runs/transport_cert.jsonl           # composition of the deferred band
```
Julia 1.12.7, Arblib 1.8.1 (default environment); BLAS is set to one thread (deterministic Float64 choices).
The machine was shared with another agent (≤ 2 threads at nice 10) during the transport run and its replay.

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
* **Lower-bound candidates actually used**: the box constants (β₀ = max(KS, rectangles, L0+AC), L0 = max(rectangle,
  Hersch)), the constant AC bootstrap from the Arb LG bound, the Taylor-model candidates above, and the pointwise
  bootstrap "bootTM": with a valid β₁(p), the LG bound L₁(p) = g₂/g₁ (Taylor-model quotient, g₁, g₂ > 0 certified)
  satisfies L₁ ≤ λ₁ pointwise, so β₂(p) = L₁(p) + 3π²/D̂(p)² ≤ λ₂(Q_p) (Andrews–Clutterbuck; iterated once more as
  "bootTM2"). The screen tries candidates in order of observed usage and keeps the first that certifies; the Arb
  certificate rebuilds exactly the logged candidate in ball arithmetic (the replay does the same).
* **Deviations from TRANSPORT_SPEC.md** (all admissible; the spec's certificate is the special case "plain
  transport, box-constant β"): first-order δ-corrections; pointwise bounds and bootstraps; MONO mode; order-4
  Taylor models with exact full products; diagonal choice; the deferred thin band E3T.

### 10.2 The run of 2026-10-01: statistics, crash/resume, replay

Command: `nice -n 5 julia -t 8 cover_t.jl runs/transport_cert.jsonl` (started 07:21:22; progress `runs/transport_cert.log`).
**Crash/resume.** The machine crashed and rebooted at ≈ 08:45; the process received SIGTERM after 96 929 processed
boxes (elapsed 5 004 s at the last progress line). It was resumed at 08:52:41 with the identical command;
`run_cover` skips logged leaves, re-expands logged splits and processes every node without a record, and finished in
163 s (2 879 boxes). The pre-reboot log is kept as `runs/transport_cert.jsonl.bak_prereboot`. `python check_log_t.py
runs/transport_cert.jsonl runs/transport_cert.jsonl.bak_prereboot` (independent of the Julia code) verifies: the backup
is a line prefix of the final log; no path is recorded twice (in particular none after the resume); every line is
valid JSON; the tree rebuilt from the 24 roots (children recomputed in IEEE double exactly as `children`) reaches every
record exactly once, every split has both children, every leaf's lo/hi equals the recomputed box bitwise, and the
leaf volumes sum to the volume of the a-priori box (31.13868, |Δ| = 3×10⁻¹³). Output: LOG STRUCTURE OK.

| quantity | value |
|---|---|
| records | 100 342 = 50 159 splits + 50 183 leaves (max depth 34) |
| certified | **31 426** — mode LG 31 421, MONO 5; lower bound: box constants 28 424, Krahn–Szegő TM 3 002; first-order transport 31 383, plain 43; n = 8 for all; diagonal 0/1: 16 491/14 935; negativity test zᵀS′z 22 264, tr 9 056, det 106 |
| E2 / E4 | 13 933 / 138 (no E1 or E3 leaves occurred) |
| COMPOUND | 3 117; their sub-leaves: E2 11 592, E3T 2 770, E4 441 (depth ≤ 6, mean 3.4) |
| SHELL | 1 (inner root [0,1/5]³×[−1/5,1/5], certified by §7 + local theorem) |
| E3T (→ Lemma B′) | 1 568 leaves + the 2 770 E3T sub-leaves of COMPOUND trees |
| FAIL / ERROR | **0 / 0** |
| leaf volume | certified 0.1416, E3T 0.1135, E2 0.7798, E4 1.2584, COMPOUND 28.842 (mostly non-convex), SHELL 0.0032 |
| wall time | ≈ 5 030 s + 163 s ≈ 1.44 h on 8 threads (machine shared; the job had ≈ 1.5–4 of the 8 cores) — prediction 1.4 h (0.7–2.5 h); predicted 1.5×10⁵ boxes / 4.6×10⁴ certificates, actual 1.0×10⁵ / 3.1×10⁴ |
| time per certified box | mean 0.67 s on one thread incl. the Float64 screen (0.25 s), median 0.47 s; 21 304 thread-s over all leaves |
| logged κ (ladder) | median 2.314, 5 % / 95 % quantiles 1.986 / 2.500; 3 262 certificates hit the ladder fallback 2.5 (only κ < 5/2 known from the cover) |

| **smallest certified margin** | 5.051e-06 (exact-rational bisection in the clean replay; largest κ* = 2.499994948580528 at `T8.b0.b0.b1.b0.d1.b0.a1.c0.d0.b0.a1.c1.d0.b0.d1.a1.b0.c0`); every certificate holds at κ = 5/2 with strict inequalities; the ladder-logged certificates have margin ≥ 4.487e-02 |

**About the logged `kappa_cert`.** Every certificate was first tested at exactly κ = 5/2 (an exact rational in Arb);
the ladder rungs were then tested at r = rationalize(kk; tol = 1e-12), but the run log records the Float64 rung kk,
a Float64 rendering of the tested rational within 10⁻¹² of it (e.g. logged 2.267625 for the tested 18141/8000, the
Float64 lying 2.2×10⁻¹⁶ below). Since the run, `cover_t.jl` logs Float64(r, RoundUp) instead. The margins quoted in
this README come from the replay's exact-rational bisection of κ, not from these renderings.

**Replay** (`nice -n 10 julia --heap-size-hint=1500M -t 3 replay_t.jl runs/transport_cert.jsonl --results runs/transport_replay_clean.jsonl`): the **clean** replay (no cached results; started 2026-10-01 15:37 from an empty results file, 3 threads at nice 10) re-verified **every one of the 50 182 non-SHELL leaves** of the frozen log (SHA-256 `cb187ebb05553ae4a45ad273e4e1304400fdc8246999e2234c1ea835eddd7bd7`): 31 426 certificates recomputed from the log alone (Float64 choice at the logged degree/diagonal/first-order flag, Arb reference integrals, Arb Taylor models, the logged lower-bound candidate rebuilt in Arb, Arb criterion at exactly κ = 5/2; for the 3 262 ladder-fallback certificates additionally a certified κ* by exact-rational bisection), and 18 756 exclusion/compound leaves re-checked by the centred-form checker (TM-only fallbacks: 0); the SHELL leaf by exact containment. Failures: **0**. Coverage rebuilt from the 24 roots: 0 uncovered, 0 mismatches, 0 unreached records. Every result line is bound to the SHA-256 of its log line, of the log and of the verifier (`4322dcc65eff1ccd…`: 16 source/config files + Julia/Arblib versions + precision). Wall time 7481 s over 1 invocation(s) (2026-10-01 15:37:42 → 2026-10-01 17:42:23 (7481 s, 50182 leaves, 3 threads)). Results file `runs/transport_replay_clean.jsonl`, SHA-256 `c9fa7f8b6dca757d1fa0c2c0d99c77f3e51928156677d53e1ee346c621422c69`; output `runs/transport_replay_clean.log`; receipt `runs/transport_replay_clean_RECEIPT.json` (summary in `runs/RECEIPT_transport.md`).

**Self-tests** (`julia -t 1 tests_t.jl runs/transport_cert.jsonl`, `runs/tests_t_output.log`): 10 886 passed, 0 failed —
incl. the exact-rational two-triangle Piola test (unequal determinants d₁ ≠ d₂, both affine maps agree on the
diagonal, equal physical normal traces while tangential traces differ, div_x w = −u), Taylor-model enclosures
against ball-arithmetic point evaluations of the transported pencil at 180 random points, the point identity with
the Rayleigh–Ritz solver at the 50 FEM validation rows, D4 invariance, lower bounds below Rayleigh–Ritz values,
and random points of the log located in the tree with their leaf condition checked.

*Superseded:* the earlier replay invocations of 09:00–10:35 (`runs/transport_replay.log`, results `runs/superseded_unbound_replay_results.jsonl`) also ended with REPLAY OK and 0 failures, but they resumed from results keyed only by path while the verifier was being edited; they are not used as evidence.

## 11. The thin band 1/10 < width/diam ≤ 2/5 (E3T) — handed to Lemma B′ (not covered here)

**Status (2026-10-01).** The E3T leaves of the transport cover (1 568 leaves, plus the 2 770 E3T sub-leaves of
COMPOUND trees) are now the subject of **Lemma B′** — every convex quadrilateral Q with 1/10 < width(Q)/diam(Q) ≤ 2/5
satisfies λ₂(Q) < (5/2)λ₁(Q) — which is being proved separately in `../bandB/` (specification `../bandB/BAND_SPEC.md`:
a rigorous cover of two 4-parameter charts normalised by the diameter, with one-dimensional Born–Oppenheimer bounds
along the diameter). Nothing in `../bandB/` was modified from here. What this directory contributes to that
hand-over: the replay re-checks, for every E3T leaf and every E3T sub-leaf, that width/diam ≤ 2/5 holds at every
point of the box (centred-form check of the polynomial condition, §10.2), so each of them lies in the domain of
Lemma B′ ∪ Lemma B (width/diam ≤ 1/10). Until Lemma B′ is proved, the global step remains incomplete.

The analysis below records why this band was not covered by the transport certificate itself.


**Definition.** E3T = {p : width(Q_p)/diam(Q_p) ≤ 2/5}, tested on a box exactly like E3 (one edge i and one vertex
pair (k,l) with cross(e_i, v_j − v_i)² ≤ (2/5)²|e_i|²|v_k − v_l|² for both j ∉ {i, i+1}, Arb Taylor models).
Boxes on which E3T holds at every point are logged with status `E3T` and are **not certified**.

**Why.** With the transport certificate, every remaining cost of the full cover sits at small width/diam: in the
pilot with first-order corrections (`runs/pilot_t4`), 98 % of the 4 686 certificates have w/D < 0.35, with median
box volumes 1.6×10⁻⁸ … 8.6×10⁻⁸ (versus 4.7×10⁻⁶ for 0.35 ≤ w/D < 0.40), and κ ≈ 2.1–2.5 even on those tiny
boxes; extrapolated ≥ 5.1×10⁵ certificates (≥ 20 h). The test functions are not the problem there; the rigorous
eigenvalue bounds are:
* no available λ₂ lower bound (Krahn–Szegő, enclosing rectangles, Andrews–Clutterbuck) is clearly above λ₁ for
  elongated shapes, so Lehmann–Goerisch is unavailable or extremely sensitive (β/λ₁ ≈ 1.01 at
  p = (0.597, 0.15, 0.15, −0.15)); the pointwise bootstrap is exact at a point (κ 2.45 → 1.47 there) but its
  Taylor-model quotient g₂/g₁ blows up already at box radius 0.003;
* the λ₁ lower bounds usable in the domain-monotonicity mode (strip, rectangle, Hersch) are 35–80 % of λ₁ for
  tapered elongated shapes, so θ₂/L0 ≈ 2.1–2.5 although the true ξ is 1.2–1.5;
* for thin shapes λ₁ ≈ π²/w² genuinely varies by tens of per cent across a box, so box constants lose twice.

**The 1-D reduction** (bound1d.jl, rigorous, box constant): λ₁(Q) ≥ μ₁(−d²/dx² + π²/h(x)²) along a direction,
with a slab-wise upper envelope of the slice length h and a Sturm shooting certificate. Measured (`scratch/t_1d.jl`):

| p0 | w/D | λ₁ (RR) | 1-D bound at p0 | r = 0.005 | r = 0.01 | strip/Hersch at p0 |
|---|---|---|---|---|---|---|
| (0.597, 0.15, 0.15, −0.15) | 0.296 | 57.73 | 56.32 (98 %) | 50.65 (88 %) | 45.86 (79 %) | 40.00 (69 %) |
| (0.321, 0.606, 1.013, 0.668) | 0.121 | 126.24 | 120.47 (95 %) | 102.71 (81 %) | 89.38 (71 %) | 97.13 (77 %) |
| (0.3, 0.5, 0.8, 0.5) | 0.226 | 55.47 | 51.91 (94 %) | 46.32 (84 %) | 41.82 (75 %) | 38.46 (69 %) |
| (0.74, 0.025, 0.15, 0.222) | 0.116 | 189.92 | 189.27 (100 %) | 160.06 (84 %) | 140.47 (74 %) | 182.09 (96 %) |
| (0.5, 0.3, 0.9, −0.2) | 0.220 | 81.67 | 34.35 (42 %) | 31.05 | 28.29 | 54.84 (67 %) |

So at a point the 1-D reduction is nearly sharp for elongated shapes (it is exact for rectangles), but as a box
constant it degrades quickly; a pointwise version (e.g. slices transported with the vertices) is what a dedicated
thin-band cover needs.

**Two cheap variants were measured and do not close it** (`scratch/t_1dbox.jl`, `scratch/t_t1d_r.jl`, 22 E3T leaf
centres of the run, w/D 0.27–0.39, ξ_RR = θ2/θ1 = 1.45–1.86):
* the 1-D reduction as a *box constant* (bound1d.jl over the box): L/θ1 = 0.82–0.92, 0.78–0.88, 0.71–0.81 at box
  radius r = 0.0025, 0.005, 0.01, i.e. a MONO-mode estimate θ2(p0)/L = 1.60–2.22, 1.68–2.31, 1.84–2.51 before the
  variation of θ2 over the box; certification therefore needs r ≈ 0.005 boxes, ≈ 6×10⁶ of them for the band
  (`scratch/band_volume.py` on the final log: the 1 568 E3T leaves have volume 0.1135 in p-space, of which 52.6 % is
  convex, outside E2 and has w/D > 1/10 — volume ≈ 0.060 to certify, not counting E3T sub-leaves of COMPOUND trees);
* a *transported* 1-D bound λ1(Q_p) ≥ c·μ1D(Q_{p0}) with c = min_K inf d_K/σ_max(A_K)² / max_K sup d_K (the
  Rayleigh-quotient comparison through the piecewise-affine map of §10): c = 0.44–0.93 at r = 0.005 and 0.17–0.85 at
  r = 0.01 (the distortion ‖A_K − I‖ is ≈ 6–20·r for these shapes), no better than the box constant.
The obstruction is structural: λ1 ≈ π²/w² varies by ±(8r/w) over a box, so any *box-constant* lower bound loses that
factor twice (once in L0, once in θ2); only a *pointwise* λ1 bound that follows the 1/w² scaling (e.g. a Taylor-model
version of the 1-D reduction) allows boxes larger than r ≈ 0.005 in this band.

**What closes the gap** (choose one):
1. a **strengthened Lemma B′**: ξ(Ω) < 5/2 for convex Ω with width/diam ≤ 2/5. The earlier thin study
   (`../explore/thin_REPORT.md`) reports that the 1-D reduction along the diameter gives this up to s ≈ 0.47
   (with an interval Sturm certificate); with B′ the present two covers + Lemma A′ + the local theorem would complete
   the convex case;
2. a dedicated cover of the band with a pointwise 1-D λ₁ bound (MONO mode suffices there, since ξ ≲ 1.5).
