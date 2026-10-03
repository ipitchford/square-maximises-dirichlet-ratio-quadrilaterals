# Lemma B′: the thin band 1/10 < width/diam ≤ 2/5 (convex quadrilaterals)

Author: Anonymous · 2026-10-01 · directory `quad-ppw-2026-09-30/bandB/` · spec: `BAND_SPEC.md`

## 0. Status (read this first)

**STATUS (2026-10-01 16:53): PROVED, computer-assisted. Re-certified with the REPAIRED code and replayed
fail-closed.**

> **Lemma B′.** Every convex quadrilateral Q with 1/10 < s(Q) := width(Q)/diam(Q) ≤ 2/5 satisfies
> λ₂(Q) < (5/2)·λ₁(Q) (Dirichlet Laplacian).

**History (LOUD).** Sol's adversarial review (`~/sol_reviews/quad_ppw_bandB_2026-10-01_review.md`) REFUTED the first
certification (cover 09:24–09:40) as a rigorous proof. The FATAL finding: in `gram_arb` the nodal products were
formed in Float64 before conversion to Arb, so the Gram balls did not enclose N and D of the chosen test functions.
The defect was verified here, repaired, and tested by exact-dyadic enclosure tests (§0a). All 44 268 CERT leaves were
then re-certified with the repaired code: none failed, and no leaf needed re-subdivision. A full fresh fail-closed
replay of the new logs passed (§12). The status read CANDIDATE between the review and that replay. The original logs
`runs/cover_chart{1,2}.jsonl` were produced by the defective code and are kept only as history; the evidence is
`runs/recert_chart{1,2}.jsonl` plus the replay results.

| | chart I (diameter = diagonal) | chart II (diameter = side) |
|---|---|---|
| leaves (exact tiling of the chart root) | 45 169 | 8 048 |
| CERT (re-certified with the repaired code) | 36 696 | 7 572 |
| exact exclusions (LEMMA_A / LEMMA_B / OUT_OF_BAND / OUTSIDE_CHART) | 110 / 141 / 8 094 / 128 | 161 / 0 / 0 / 315 |
| FAIL_FLOOR (refined below width 2⁻²⁴) / re-subdivided in re-certification | **0** / **0** | **0** / **0** |
| worst certified margin max R₂/((5/2)E) (rigorous upper bound, fresh replay) | 0.9999985269 | 0.9999679033 |

* **Replay (fresh, repaired code; §12): PASSED, exit 0** (16:43–16:53, 3 threads, 583 s). Both charts tile exactly;
  8 949 exclusions were re-derived and 44 268 certificates re-run with 0 errors. Results are bound to the record and
  source hashes. The fail-closed test passes 24/24.
* **Inputs taken as proved, not re-proved here:** paper Lemma A (near-triangles; it rests on the triangle theorem of
  Siudeja 2010 and Arbon–Mannan–Psenka–Ragavan 2020) and paper Lemma B (s ≤ 1/10, analytic); both in
  `../paper/PAPER.md` §2.2–2.3.
* **Trusted computing base:** Julia 1.12.7, FLINT/Arb ball arithmetic via Arblib.jl (128-bit), GMP rationals
  (`Rational{BigInt}`), and the code in this directory. The code is supported by a symbolic check of every formula,
  exact-rational and exact-dyadic enclosure tests (`tests.jl`, T1–T7), a 2-D FEM comparison and the replay. The
  replay re-executes the same certificate code (`cert.jl`); it is not a second implementation of the certificate (the
  exclusions *are* re-derived by separately written code). The coordinator's checker (`../verify/`) is the independent
  check. It is not part of this directory or its MANIFEST, and the Python-checker bug in the same review
  (`verify/xband.py`) is the coordinator's and is being fixed there.
* **Not rigorous, and not used by the proof:** the floating-point choices of E, the meshes and the Ritz test
  functions (every choice is admissible; only the verified inequalities matter), and the FEM validation.
* **Consequence:** together with Lemma B (s ≤ 1/10), every convex quadrilateral with s ≤ 2/5 has ξ < 5/2. This
  discharges the boxes that `../global2` logs as `E3T` (README_global2 §11).
* **Deviations from BAND_SPEC (both proved in §3):** chart I `OUT_OF_BAND` uses the exact width formula (the spec's
  bound s ≥ y3 cannot exclude anything, because y3 ≤ 2/5 on the whole chart-I root, while chart I contains points
  with s up to ≈ 0.62 where the 1-D bound fails); chart I `LEMMA_B` also uses an edge-normal width upper bound.

### 0a. Response to the 2026-10-01 adversarial review (every claim re-checked here before acting)

| finding (Sol) | verified here? | repair | evidence |
|---|---|---|---|
| **FATAL**: Float64 nodal products in `gram_arb` (cert.jl:191–195) | **yes**: Sol's shifts reproduce exactly on leaf `AbcdAbcDABDabDAbd`. N₁₁ +4.3762·10⁻¹⁶, N₁₂ −9.4730·10⁻¹⁵, N₂₂ −5.5239·10⁻¹⁵, D₁₁ −2.2671·10⁻²⁰ against radii 10⁻³⁴…10⁻³⁸; the old and corrected balls are disjoint. Corrected det((5E/2)D − N) = 5.9652755024·10⁻⁵ (chart I) and 1.5591293491·10⁻³ (chart II, `ABCdABcDADcCbaCDAc`), as Sol found | nodal values are promoted to Arb before any product; element integrals are factored into `element_forms`, which converts every Float64 input before any arithmetic. Audit of the whole acceptance chain: no other Float64 arithmetic enters a ball (H_k and Ritz node values are directed-rounded bounds; E, meshes and f are choices) | `tests.jl` T7: 2 896 STRICT containment checks against exact `Rational{BigInt}` references (polynomial parts, constant-h and end elements) and 2048-bit references (log and π parts). Cases: random elements of every type, `Iints_arb` on both branches with ~200 cases straddling ρ = 1/4, random meshes, and 12 real certificate meshes. The pre-review assembly (kept as `gram_arb_v1`) fails containment in 72/72 meshes |
| **SERIOUS**: Float64 remainder radius in `Iints_arb` (cert.jl:120–143) | **yes**: ρ̄, K and the tail were Float64 with empirical cushions | branch test ρ̄ ≤ 1/4 by exact Arf comparison; K fixed per bracket (16/24/32/48 for ρ̄ ≤ 1/64, 1/16, 1/8, 1/4); tail ρ̄^K/(1 − ρ̄)·(1/h_m)·∫w computed in Arb and added as an Arf upper bound | T7(b) |
| **SERIOUS**: `run_cover.sh` masks failures | **yes** (prints exit codes, exits 0) | `set -euo pipefail` plus an explicit failure exit; the same for `run_recert.sh` and `run_replay.sh` | — |
| **MINOR**: summary regex ends in `.*` (replay.jl:38) | **yes** (extra fields and duplicate count keys passed) | every summary field enumerated in order; count keys whitelisted and strictly sorted; subtree count and max_depth cross-checked | fail-closed test 24/24, incl. extra field, duplicate key, wrong subtree count, wrong max_depth |
| SERIOUS/MINOR: Python checker `h_out` bug; unbound checker evidence; checker ignores metadata | not this directory | the coordinator's code (`../verify/`), being fixed there; not edited here | — |

**Re-certification** (`recert.jl`, `run_recert.sh`, 16:29–16:36, 3 threads). Every leaf path of the original logs
was re-processed by the cover's own `process_subtree` with the repaired `cert.jl`: exclusions were re-derived, CERT
leaves re-certified, and anything that failed would have been bisected further. Result: 45 169 → 45 169 and
8 048 → 8 048 leaves, **0 re-subdivided, 0 status changes**; worst margins 0.9999985269 and 0.9999679033 (unchanged
to ~10⁻¹⁵, the size of the defect). The defective `cert.jl` is archived as `runs/cert_v1_defective.jl.txt`.

**Field order of `lo`/`hi` (for independent checkers): chart I = (y3, t, x1, x3); chart II = (y2, t, x2, x3).**
The full log format is given in §9.

## 1. Notation

Q is a convex quadrilateral (non-degenerate: four vertices, all interior angles < π). D = diam Q, w = minimal width,
s = w/D, ξ = λ₂/λ₁. For vertices labelled counter-clockwise, c_i = (P_{i+1} − P_i) × (P_{i+2} − P_i) (indices mod 4)
is twice the area of the triangle P_iP_{i+1}P_{i+2}, so c_i > 0 for all i. The quantities s, ξ and the sign
pattern of the c_i are invariant under orientation-preserving similarities.

## 2. Chart reduction (BAND_SPEC §1)

**Charts.** Chart I: p = (y3, t, x1, x3), P0 = (0,0), P1 = (x1, −t·y3), P2 = (1,0), P3 = (x3, y3), root
R_I = [1/20, 2/5] × [0,1] × [0,1] × [0,1/2]. Chart II: p = (y2, t, x2, x3), P0 = (0,0), P1 = (1,0), P2 = (x2, y2),
P3 = (x3, t·y2), root R_II = [1/10, 2/5] × [0,1]³. Q(p) = conv(P0, …, P3). A point p is *valid* if Q(p) is convex
with c_i > 0 in the order P0P1P2P3, all six vertex distances are ≤ 1, and (chart II) x3 < x2.
Closed forms (`sympy_check.py` 4a–4g): chart I c = (t·y3, y3(1 + t − x1 − t·x3), y3, y3(x1 + t·x3));
chart II c = (y2, y2(1 − x3) − t·y2(1 − x2), t·y2·x2 − x3·y2, t·y2).

**Lemma 1.** Every convex quadrilateral with 1/10 < s ≤ 2/5 is similar to Q(p) for a valid p ∈ R_I (chart I) or
a valid p ∈ R_II (chart II).

*Proof.* The diameter of a polygon is attained at a pair of vertices {P, P′}. Label the vertices counter-clockwise.

*Chart I (P, P′ opposite).* Name them P0 = P, P2 = P′ and the others P1, P3 in counter-clockwise order. The
orientation-preserving similarity z ↦ (z − P0)/(P2 − P0) (complex notation) sends P0 ↦ (0,0), P2 ↦ (1,0), keeps the
order counter-clockwise and the c_i positive, and makes D = 1, so all six distances are ≤ 1. Now c0 = −y(P1) > 0
and c2 = y(P3) > 0: write P1 = (x1, −y1), P3 = (x3, y3) with y1, y3 > 0. From |P_k P0| ≤ 1 and |P_k P2| ≤ 1:
1 − √(1 − y²) ≤ x ≤ √(1 − y²) for both other vertices, so x1, x3 ∈ (0, 1).
WLOG y1 ≤ y3: the reflection y ↦ −y fixes P0, P2 and, after relabelling counter-clockwise from P0, exchanges
(x1, y1) ↔ (x3, y3). WLOG x3 ≤ 1/2: the reflection x ↦ 1 − x exchanges P0 ↔ P2; relabelling counter-clockwise from
the new origin gives (x1, x3) ↦ (1 − x1, 1 − x3) with y1, y3 unchanged (so the first WLOG survives).
*Width bounds.* Q ⊇ T = P0P2P3. The minimal width of a triangle is its smallest altitude (Lemma 3 below: the
minimum is attained normal to an edge, where the width is the altitude), i.e. 2·area/(longest side) = y3 because
P0P2 (length 1) is the longest side. Width is monotone under inclusion (width in each direction is), so s ≥ y3.
Q lies in the strip −y1 ≤ y ≤ y3, so s ≤ y1 + y3 ≤ 2y3. Hence s ≤ 2/5 gives y3 ≤ 2/5 and s > 1/10 gives
y3 > 1/20. With t = y1/y3 ∈ (0, 1], p = (y3, t, x1, x3) ∈ R_I and p is valid.

*Chart II (P, P′ adjacent).* Name them so that P0 = P, P1 = P′ are consecutive counter-clockwise, then P2, P3. The
similarity z ↦ (z − P0)/(P1 − P0) gives P0 = (0,0), P1 = (1,0), D = 1; c0 = y(P2) > 0 and c3 = y(P3) > 0.
As above, x2, x3 ∈ (0,1). WLOG y3 ≤ y2: the reflection x ↦ 1 − x with the counter-clockwise relabelling
(P0, P1, P2, P3) ← (P1, P0, P3, P2) gives P2′ = (1 − x3, y3), P3′ = (1 − x2, y2). Then x3 < x2, because
c2 = x2·y3 − x3·y2 > 0 gives x3·y2 < x2·y3 ≤ x2·y2. *s = y2 exactly:* Q ⊇ P0P1P2, whose longest side is P0P1 = 1,
so w ≥ 2·area/1 = y2; and Q lies in the strip 0 ≤ y ≤ max(y2, y3) = y2. So y2 = s ∈ (1/10, 2/5], t = y3/y2 ∈ (0,1],
and p ∈ R_II is valid. □

## 3. Exclusions (BAND_SPEC §2), decided exactly in rationals on the whole closed box

**Lemma 2 (multilinear ranges).** If F is affine in each coordinate separately, then its minimum and maximum over a
box are attained at corners. (Fix all coordinates but one: F is affine in it, so the extremes lie at endpoints;
induct.) All vertex-coordinate differences, all c_i and all c_i − (7/200)c_{i+2} are multilinear in both charts
(`sympy_check.py` 4d, 4h). For |P_j − P_k|² = g₁² + g₂² with multilinear g's, the exact ranges [gᵢ⁻, gᵢ⁺] give the
bounds Σ min_{[gᵢ⁻,gᵢ⁺]} g² ≤ |P_j − P_k|² ≤ Σ max_{[gᵢ⁻,gᵢ⁺]} g².

**Lemma 3 (width of a convex polygon).** w(K) = min over the edges e of the largest vertex distance from the line
through e. For a quadrilateral this is w = min_i max(c_i, c_{i+3})/|e_i| with e_i = P_{i+1} − P_i, because
e_i × (P_{i+2} − P_i) = c_i and e_i × (P_{i+3} − P_i) = c_{i+3} (`sympy_check.py` 4i).
*Proof.* Let W(θ) be the width in direction u(θ). The edge normals and their negatives cut the circle of directions
into closed arcs. On each arc the maximising vertex v and the minimising vertex v′ of ⟨·, u⟩ are fixed, so
W(θ) = ⟨v − v′, u(θ)⟩ = |v − v′| cos(θ − θ₀) ≥ 0. A cosine is concave where it is non-negative, so W attains its
minimum over the arc at an endpoint, i.e. in a direction normal to an edge. There one supporting line contains the
edge and W is the largest vertex distance from it. □

**Rules** (`geom.jl: exclusion`; re-derived in `replay.jl: ex_verify` from the closed forms). Each rule is a single
inequality that holds on the whole closed box (Lemma 2):

| status | test (chart) | why no valid in-band point of the box needs a certificate |
|---|---|---|
| OUTSIDE_CHART | min_B \|P_jP_k\|² > 1 for some j < k; (II) max_B c_i ≤ 0 for some i, or x3⁻ ≥ x2⁺ | no point of B is valid |
| LEMMA_A | max_B (c_i − (7/200)c_{i+2}) ≤ 0 for some i | ξ ≤ (7/3)(1 + 7/200)² = 2.4995… < 5/2 (paper Lemma A) |
| LEMMA_B | (I) y3⁺(1 + t⁺) ≤ 1/10, or for some i: M := max(max_B c_i, max_B c_{i+3}) ≥ 0 and M² ≤ (1/100)·min_B \|e_i\|²; (II) y2⁺ ≤ 1/10 | s ≤ y1 + y3, resp. s ≤ max(c_i, c_{i+3})/\|e_i\| (Lemma 3), resp. s = y2: so s ≤ 1/10 (outside the band; ξ < 5/2 by Lemma B anyway) |
| OUT_OF_BAND | (I) y3⁻ > 2/5, or for every i: m_i := max(min_B c_i, min_B c_{i+3}) > 0 and m_i² > (4/25)·max_B \|e_i\|²; (II) y2⁻ > 2/5 | s > 2/5 by Lemma 3 / s ≥ y3 / s = y2 |

## 4. Profiles and containment (BAND_SPEC §3)

Write a convex Q(p) as {0 < x < 1, α_Q(x) < y < β_Q(x)} (its interior) and h_Q = β_Q − α_Q (projections lie in
[0,1] because all vertex abscissae lie in [0,1] on both roots). On a box B, a coordinate's range is [·⁻, ·⁺];
y1 ∈ [t⁻y3⁻, t⁺y3⁺] (chart I) and y3 ∈ [t⁻y2⁻, t⁺y2⁺] (chart II), since t, y ≥ 0.

Tents. Let τ(x; a) = x/a on [0,a] and (1 − x)/(1 − a) on [a,1]. For fixed x, a ↦ τ(x; a) increases for a ≤ x and
decreases for a ≥ x. Hence
τ⁺(x; l, u) := max_{a∈[l,u]} τ(x; a) = min(1, x/l, (1 − x)/(1 − u)) and
τ⁻(x; l, u) := min_{a∈[l,u]} τ(x; a) = min(x/u, (1 − x)/(1 − l)), the latter with its kink at u/(1 + u − l)
(`sympy_check.py` 5a; case analysis x < l, l ≤ x ≤ u, x > u). Both are minima of affine functions, hence concave
and piecewise linear.

**Lemma 4 (chart I).** For p ∈ B: α_Q = −y1·τ(·; x1) and β_Q = y3·τ(·; x3). The edges P0P1, P1P2 lie below the axis
and P2P3, P3P0 above it, and with x1, x3 ∈ [0,1] every vertical line meets one lower and one upper edge. Hence
h_Q ≤ h_out := y1⁺τ⁺(·; x1⁻, x1⁺) + y3⁺τ⁺(·; x3⁻, x3⁺), and α_Q ≤ a_in := −y1⁻τ⁻(·; x1⁻, x1⁺) <
b_in := y3⁻τ⁻(·; x3⁻, x3⁺) ≤ β_Q on (0,1), i.e. Q_in := {a_in < y < b_in} ⊂ Q. Here h_out is concave, a_in convex,
b_in concave, and all are piecewise linear with exact rational breakpoints {x1^±, x3^±} and the two τ⁻ kinks.

**Lemma 5 (chart II).** For convex Q(p), p ∈ B: Q ⊂ {y ≥ 0}, its lower boundary is P0P1, and β_Q is the upper
concave envelope of the four vertices. *Outer:* P_k lies in the convex hull of the 4 corners of its rectangle
[x_k⁻, x_k⁺] × [y_k⁻, y_k⁺], so Q ⊂ conv{(0,0), (1,0), 8 corners} and h_Q ≤ h_out := the upper concave envelope of
these 10 points (its lower boundary is y = 0). *Inner:* the triangle P0P1P_k equals {0 ≤ y ≤ y_k·τ(x; x_k)}, so
T_k := ∩_{p∈B} P0P1P_k = {0 ≤ y ≤ y_k⁻τ⁻(x; x_k⁻, x_k⁺)}, which is the triangle with apex
(x_k⁺, y_k⁻)/(1 + x_k⁺ − x_k⁻) (`sympy_check.py` 5b). Q ⊇ T_2 ∪ T_3 and Q is convex, so
Q ⊇ conv(T_2 ∪ T_3) = {0 ≤ y ≤ b_in} with b_in the upper concave envelope of (0,0), (1,0), apex₂, apex₃; a_in = 0.
(The envelopes are computed exactly by a monotone-chain hull in rationals; `geom.jl: upper_envelope`.)

**Lemma 6 (the data the code actually uses).** (i) *Sturm cells.* For float nodes 0 = ξ₀ < … < ξ_K = 1, H_k is the
upward-rounded ball upper bound of max(h_out(ξ_k), h_out(ξ_{k+1}), h_out(X) for breakpoints X ∈ (ξ_k, ξ_{k+1})). A
piecewise linear function attains its maximum over an interval at an endpoint or at a breakpoint, so
H_k ≥ max_{[ξ_k, ξ_{k+1}]} h_out ≥ h_Q there. Breakpoints are located by exact Float64/Rational comparisons
(Julia `Base` compares these exactly; `tests.jl` T6). (ii) *Ritz profile.* For float nodes ζ_j, take
a_j = (ball upper bound of a_in(ζ_j)) rounded up and b_j = (ball lower bound of b_in(ζ_j)) rounded down, with
a = b = 0 exactly at x = 0, 1. Since a_in is convex, its chords lie above it, and raising the chord ends keeps them
above; so the P1 interpolant a_int ≥ a_in, and symmetrically b_int ≤ b_in. The code checks b_j − a_j > 0 at interior
nodes (balls), so Q_int := {0 < x < 1, a_int < y < b_int} ⊂ Q_in ⊂ Q.

## 5. Lower bound: 1-D reduction and Sturm shooting

**Lemma 7 (1-D reduction; cf. `../global2/bound1d.jl`).** Let Q be convex with projection ⊂ [0,1] and let V be
bounded with V ≤ π²/h_Q² wherever h_Q > 0. Then λ₁(Q) ≥ μ₁(V) := inf_{v∈H¹₀(0,1)} ∫(v′² + Vv²)/∫v².
*Proof.* Take φ ∈ C_c^∞(Q) and v(x) = (∫φ(x,y)²dy)^{1/2}. On each slice, φ(x,·) is compactly supported in an
interval of length h_Q(x), so ∫|∂_yφ|²dy ≥ (π/h_Q)²∫φ²dy. The functions v_ε = (v² + ε²)^{1/2} − ε are smooth with
|v_ε′| = |∫φ∂_xφ|/(v² + ε²)^{1/2} ≤ (∫|∂_xφ|²dy)^{1/2}, so v is Lipschitz, compactly supported in (0,1), and
|v′| ≤ (∫|∂_xφ|²dy)^{1/2} a.e. Hence ∫_Q|∇φ|² ≥ ∫₀¹(v′² + Vv²)dx ≥ μ₁(V)∫₀¹v² = μ₁(V)∫_Qφ². Conclude by density. □

The code uses V = π²/H_k² on [ξ_k, ξ_{k+1}], which satisfies V ≤ π²/h_Q² by Lemma 6(i).

**Lemma 8 (Sturm criterion).** Let V be bounded and piecewise constant on [0,1], E ∈ ℝ, and v the solution of
−v″ + (V − E)v = 0, v(0) = 0, v′(0) = 1. If v > 0 on (0,1], then μ₁(V) > E.
*Proof.* Suppose μ₁ ≤ E and let ψ be the first eigenfunction: ψ > 0 on (0,1), ψ(0) = ψ(1) = 0, and ψ′(1) < 0 (if
ψ′(1) = 0, ODE uniqueness would give ψ ≡ 0). The Wronskian W = vψ′ − v′ψ is continuous and piecewise C¹ with
W′ = vψ″ − v″ψ = (E − μ₁)vψ ≥ 0 and W(0) = 0, so W(1) ≥ 0. But W(1) = v(1)ψ′(1) < 0. Contradiction. □

**Lemma 9 (cell rules, `cert.jl: sturm_ok`).** Let q_k = V_k − E on J_k = [ξ_k, ξ_{k+1}]. Suppose that for every k,
v(ξ_{k+1}) > 0; v(ξ_k) > 0 for k ≥ 1 (for k = 0, v′(0) = 1 > 0); and either q_k > 0, or q_k < 0 with
√(−q_k)·|J_k| < π. Then v > 0 on (0,1].
*Proof.* On J_k, v″ = q_k v. If q_k > 0, a non-trivial solution has at most one zero and changes sign there, so
positivity at both ends excludes a zero; on J₀, v = sinh(κs)/κ > 0. If q_k < 0 with ω = √(−q_k), then
v = R cos(ωs − φ) has simple zeros spaced π/ω > |J_k| apart, so it has at most one zero, with a sign change, and
none if it is positive at both ends; on J₀, v = sin(ωs)/ω > 0 for 0 < s ≤ |J₀|. □

The transfer over J_k uses (cosh, sinh/κ; κ sinh, cosh) or (cos, sin/ω; −ω sin, cos), evaluated in 128-bit Arb
balls. The cell length is exact (difference of doubles), H_k is an exact double, and π is an Arb enclosure. A cell
whose q_k-ball contains 0 is treated as failure; E is then lowered and the check repeated.
E = (1 − 10⁻⁷)·μ₁^float(step potential); on failure E is multiplied by (1 − 10⁻⁴), (1 − 10⁻³), (1 − 10⁻²) in up
to three retries.

## 6. Upper bound: weighted Ritz on Q_int

**Lemma 10.** For a P1 function f on the nodes ζ with f(0) = f(1) = 0, set u_f = f(x)·sin(π(y − a_int)/h_int),
h_int = b_int − a_int. Then u_f ∈ H¹₀(Q_int), with
∫|∇u_f|² = N(f) := ∫₀¹ [f′²h/2 + f f′h′/2 + f²(π²/h)(1/2 + a′b′/2 + h′²(1/6 + 1/(4π²)))] dx and
∫u_f² = D(f) := ∫₀¹ f²h/2 dx (a = a_int, b = b_int, h = h_int; `sympy_check.py` 1a–1f re-derives this for general
smooth a, b). For f₁, f₂ with D_ij ≻ 0: λ₂(Q) ≤ λ₂(Q_int) ≤ R₂ := the larger root of det(N_ij − R·D_ij) = 0.
*Proof.* On interior elements h_int ≥ min(h_j, h_{j+1}) > 0, so u_f is smooth with bounded derivatives there. On the
end elements f, a_int and h_int are linear and vanish at the end, so f/h_int is constant there, and with it ∂_y u_f
and the term f·π cos(·)(a′ + ηh′)/h of ∂_x u_f. So u_f is Lipschitz on the closure of Q_int and vanishes on
∂Q_int, hence u_f ∈ H¹₀(Q_int). The identities follow by integrating in y on each slice (Fubini). Domain monotonicity
and min–max on span(u_{f₁}, u_{f₂}) finish the proof. □

The coefficient 1/2 + a′b′/2 + h′²C₆ (C₆ = 1/6 + 1/(4π²)) is in fact ≥ 1/2 (`sympy_check.py` 5c–5d). On each
element (`cert.jl: gram_arb`; formulas `sympy_check.py` 2a–2c, 3a–3f):
∫f′g′h/2 = Δf·Δg(h_a + h_b)/(4L); ∫(f′g + fg′)h′/4 = (h′/4)(f_b g_b − f_a g_a);
∫fgh/2 = (L/24)[f_a g_a(3h_a + h_b) + (f_a g_b + f_b g_a)(h_a + h_b) + f_b g_b(h_a + 3h_b)];
∫fg/h = L[f_a g_a I₀₀ + (f_a g_b + f_b g_a)I₀₁ + f_b g_b I₁₁], with I₀₀, I₀₁, I₁₁ = ∫₀¹(1−s)², s(1−s), s² over
h_a(1−s) + h_b s. Let ρ̄ be the exact Arf upper bound of |h_b − h_a|/(h_a + h_b) = |ε|/2 (ε = (h_b − h_a)/h_m). If ρ̄ ≤ 1/4 (exact
Arf comparison), the integrals use the series of 1/(1 + εσ) (σ = s − 1/2) with exact rational coefficients, a fixed
order K = 16/24/32/48 for ρ̄ ≤ 1/64, 1/16, 1/8, 1/4, and the tail ρ̄^K/(1 − ρ̄)·(1/h_m)·∫w computed in Arb and added as
an Arf upper bound. Otherwise they use the closed forms with log(h_b/h_a), in Arb. On the two end elements h and f vanish together, so
∫fg/h = L f g/(2h) exactly (no log). f₁, f₂ are float Ritz vectors (lumped-mass eigenvectors plus one step of
block inverse iteration); only N_ij, D_ij are enclosed rigorously. **All nodal values, profile values and node
abscissae are converted to Arb before any arithmetic** (`element_forms`, `gram_arb`; the pre-review code multiplied
nodal values in Float64, §0a).

## 7. The certificate lemma

**Lemma 11.** For a box B, suppose that (i) `sturm_ok` holds at E for V_k = π²/H_k²; (ii) D ≻ 0 (D₁₁ > 0,
det D > 0 as balls); and (iii) (5/2)E·D − N ≻ 0 (leading minors positive as balls). Then
λ₂(Q) ≤ R₂ < (5/2)E < (5/2)λ₁(Q) for every convex Q(p), p ∈ B.
*Proof.* λ₁(Q) ≥ μ₁(V) > E by Lemmas 6(i), 7, 8 and 9. For a generalized eigenpair Nw = RDw,
wᵀ((5/2)E·D − N)w = ((5/2)E − R)·wᵀDw > 0, so every Ritz value, in particular R₂, is < (5/2)E. Lemmas 6(ii) and 10
give λ₂(Q) ≤ R₂. □

(The code also encloses R₂ = (b + √(b² − 4 det D det N))/(2 det D), b = N₁₁D₂₂ + N₂₂D₁₁ − 2N₁₂D₁₂, for the log.
The certificate decision itself is (iii).)

## 8. The cover and the proof of Lemma B′

Each chart root is pre-split into 256 subtrees (8 round-robin bisections in coordinate order 1,2,3,4,1,2,3,4) and
refined depth-first by exact midpoint bisection. A box is excluded if one rule of §3 holds on it, certified if
Lemma 11 holds, and otherwise bisected in the coordinate whose collapse to its midpoint lowers the float ratio
R₂/E most (the widest relative coordinate if the box centre itself has ratio ≥ 2.45). A box that would need
bisection below width 2⁻²⁴ is logged `FAIL_FLOOR`; none occurred (maximum depth 21).

*Proof of Lemma B′.* Let Q be convex with 1/10 < s ≤ 2/5. By Lemma 1, Q is similar to Q(p) with p valid in R_I or
R_II. The leaves of that chart tile its root exactly (§12), so p lies in a leaf L, and:
CERT ⇒ ξ < 5/2 (Lemma 11); LEMMA_A ⇒ ξ ≤ 2.4995… < 5/2 (paper Lemma A); LEMMA_B ⇒ s ≤ 1/10, impossible;
OUT_OF_BAND ⇒ s > 2/5, impossible; OUTSIDE_CHART ⇒ p not valid, impossible. There are no FAIL_FLOOR leaves. ξ is
similarity-invariant. □

## 9. Implementation, rigour map and log format

| step | file / function | arithmetic |
|---|---|---|
| exclusions | `geom.jl: exclusion` (replay: `replay.jl: ex_verify`, separate code) | exact `Rational{BigInt}` |
| profiles h_out, a_in, b_in, envelopes | `geom.jl: profiles, upper_envelope` | exact rationals |
| meshes, E, f₁, f₂, bisection choice | `cert_float.jl`, `cert.jl: float_eval` | Float64 (choices only) |
| H_k, Ritz node values a_j, b_j | `cert.jl: sturm_cells, inner_nodes` | Arb balls from exact rationals, rounded outward |
| Sturm check | `cert.jl: sturm_ok` | Arb, 128 bits |
| Gram matrices, D ≻ 0, (5/2)ED − N ≻ 0 | `cert.jl: element_forms, gram_arb, Iints_arb, certify` | Arb, 128 bits; every Float64 input converted before any arithmetic; series tail by Arb/Arf upper bounds |
| tiling | `tiling.jl: check_tiling` | exact (paths) + Σ exact volumes |

Mesh parameters (`cert.jl: Params`): Sturm cells equidistribute |d log h|·(h/h_max)³ with weight ≤ 0.002 and length
≤ 0.02 (50–372 cells per box); Ritz elements weight 0.02, p = 2, length ≤ 0.008 (125–144 elements). Breakpoints are
nodes (rounded to the nearest double; Lemma 6 makes this irrelevant to rigour).

**Log format** (authoritative logs `runs/recert_chart{1,2}.jsonl`; the historical `runs/cover_chart{1,2}.jsonl` have
the same format; one JSON object per line, fields in exactly this order):
* exclusion leaf: `{"chart":c,"path":P,"status":S,"rule":R,"lo":[4 × "p/q"],"hi":[4 × "p/q"]}` with
  S ∈ {OUTSIDE_CHART, LEMMA_A, LEMMA_B, OUT_OF_BAND} and R from a fixed whitelist (`replay.jl: RULES`);
* certified leaf: `{"chart":c,"path":P,"status":"CERT","lo":[…],"hi":[…],"E":…,"R2":…,"margin":…,"nS":…,"nR":…,"sec":…}`.
  E is the next double below the E used, so the decimal is a valid lower bound for λ₁. R2 and margin = R₂/((5/2)E) are
  rounded up. All are decimal strings, and all are informational: the replay ignores them;
* subtree marker `{"subtree":P8,"done":true,"leaves":n,"sec":…}` after the leaves of each completed depth-8 subtree;
  one final summary line with exactly these fields in this order: `chart, leaves, counts{status:n, keys sorted},
  worst_margin_R2_over_2.5E, max_depth, max_nS, max_nR, cert_seconds, subtrees, subtree_thread_seconds,
  wall_seconds_this_session, threads, finished`.
* rule whitelist (exact strings; indices are the paper's c_i, P_i, i = 0..3):
  OUTSIDE_CHART: `|PjPk| > 1` (0 ≤ j < k ≤ 3, both charts), and in chart II also `ci <= 0` (i = 0..3) and
  `x3 >= x2`; LEMMA_A: `ci <= 7/200 c(i+2 mod 4)`; LEMMA_B: chart I `y1+y3 <= 1/10` or
  `width normal to edge i <= 1/10`, chart II `y2 <= 1/10`; OUT_OF_BAND: chart I `y3 > 2/5` or
  `width lower bound > 2/5`, chart II `y2 > 2/5`. Each string names the inequality of §3 that holds on the whole box.
* lo/hi order: **chart I (y3, t, x1, x3); chart II (y2, t, x2, x3)**. Rationals are canonical `p/q` (q > 0,
  gcd = 1; zero is `0/1`).
* path: one character per bisection; `'a'+(d−1)` is the lower half and `'A'+(d−1)` the upper half of coordinate d
  (1-based, order as above), split at the exact midpoint. Root boxes R_I, R_II as in §2. Every leaf path starts with
  one of the 256 depth-8 prefixes (bisection coordinates 1,2,3,4,1,2,3,4).

## 10. Run statistics

The leaf tree of the re-certification is identical to the original cover's (0 re-subdivisions). The cover statistics
below are for that tree; the margins are the re-certified ones.

| | chart I | chart II |
|---|---|---|
| leaves / CERT | 45 169 / 36 696 | 8 048 / 7 572 |
| LEMMA_A / LEMMA_B / OUT_OF_BAND / OUTSIDE_CHART / FAIL_FLOOR | 110 / 141 / 8 094 / 128 / 0 | 161 / 0 / 0 / 315 / 0 |
| volume fractions CERT / OUT_OF_BAND / LEMMA_B / OUTSIDE_CHART / LEMMA_A | 0.820 / 0.116 / 0.048 / 0.013 / 0.003 | 0.398 / 0 / 0 / 0.572 / 0.030 |
| margin R₂/((5/2)E): max / 99 % / 90 % / median / min | 0.9999985 / 0.9995 / 0.9953 / 0.9798 / 0.8533 | 0.9999679 / 0.9993 / 0.9937 / 0.9744 / 0.8437 |
| depth min / median / max | 8 / 17 / 21 | 8 / 16 / 21 |
| Sturm cells / Ritz elements per box | 140–368 / 125–144 | 50–372 / 125–144 |
| original cover: wall time, threads (nice 10) | 771 s, 6 | 164 s, 2 |
| re-certification (repaired code): wall time, threads | 327 s, 3 | 90 s, 3 |
| fresh replay phase 2 (3 threads) | 455 s | 122 s |

Chart I's original cover ran with 6 threads after `cover_t.jl` had finished, as BAND_SPEC §6 allowed at the time;
everything after the review ran with 3 threads (coordinator limit). Margins close to 1 are expected: a box is
certified as soon as Lemma 11 holds, so many certificates have little slack. Every certified inequality is still
strict and decided in ball arithmetic. Statistics: `log_stats.py`.

## 11. Validation (done before the rigorous run; logs kept)

* **(a) 1-D bounds vs 2-D FEM** (`validate_fem.py` → `validation_fem.csv`; `validate_bounds.jl` →
  `validation_bounds.log`). 130 shapes: 60 random valid in-band points per chart, plus the thin study's adversarial
  maximisers, three shrinks of the cap-0.40 maximiser into the band, and three near-Lemma-A triangles (62 chart I,
  68 chart II). Rigorous point bounds (E, R₂) were compared with P2-FEM λ₁, λ₂ (k = 6; k = 7 on every 10th shape,
  |k6 − k7| ≤ 1.1·10⁻⁴ relative). **Bracket violations: 0.** The minimum of λ₁^FEM/E − 1 is 1.0·10⁻³ and the minimum
  of R₂/λ₂^FEM − 1 is 3.0·10⁻⁴, both at a very elongated shape (s = 0.1002, ξ = 1.046, where the 1-D reduction is
  nearly exact); typical margins are several per cent. FEM λ₂ is itself an upper approximation, so R₂ ≥ λ₂^FEM is
  on the safe side. Max R₂/E = 2.2305; all 130 point certificates pass.
* **(b) Symbolic check** (`sympy_check.py` → `sympy_check.log`): ALL PASS. It covers the weighted form for general
  a, b, the four elementary integrals, all element formulas, the log closed forms, the series coefficients and their
  convergence, the end elements, every c_i formula and the multilinearity of all exclusion quantities, the width
  identities, the τ⁻ kink and apex, and c_e ≥ 1/2.
* **(c) Unit tests** (`tests.jl 400` → `tests.log`, run with the repaired code): ALL PASS. T1 makes 87 254 exact-rational slice comparisons at
  random points of random cover-shaped boxes (both charts), checking h_Q ≤ h_out, a_Q ≤ a_in < b_in ≤ b_Q, that the
  inward-rounded Ritz profile lies inside Q_in, and that h_Q ≤ H_k on Sturm cells. T2 makes 20 210 exact point checks
  of every exclusion status. T3 checks 621 enclosures of I₀₀, I₀₁, I₁₁ against BigFloat references, including
  ρ ≈ 1/4, Δ = 0 and tiny Δ. T4 runs Sturm on a constant potential (exact μ₁) and compares the Gram matrices with
  float assembly. T5 runs the tiling checker on valid and corrupted trees. T6 checks exact Float64/Rational
  comparisons. **T7** (added after the review) is the exact-dyadic enclosure test of §0a: 2 896 strict containment
  checks of `element_forms`, `Iints_arb` and the full `gram_arb` against exact `Rational{BigInt}` and 2048-bit
  references, with no tolerance. The pre-review assembly fails it on 72/72 meshes. T4's float-versus-Arb comparison
  (tolerance 10⁻⁹) is only a sanity check; it could not detect the defect, which is why T7 exists.
* (a) was re-run with the repaired code (`validation_bounds.log`): unchanged, 0 bracket violations, 130/130 point
  certificates. The pre-fix log is kept as `runs/validation_bounds_v1_defective_code.log`.

## 12. Replay (fail-closed)

`replay.jl` (run: `./run_replay.sh` → `replay.log`; bound results `runs/replay_chart{1,2}.jsonl`). It trusts only the
leaf paths and statuses, and re-verifies those:
* **strict schemas:** every line must match exactly one full-line schema: a leaf of a whitelisted status with its
  field list and order, a whitelisted rule string, canonical rationals, the depth-8 subtree marker, or the single
  final summary line, whose every field is enumerated (§9). Unknown statuses, unknown rules, malformed lines, extra
  fields, duplicate or unsorted summary keys and trailing lines all FAIL; a `FAIL_FLOOR` record is itself a failure;
* the box is **recomputed from the path** and must equal the logged lo/hi exactly;
* **duplicate lines and duplicate paths** FAIL; the paths must **tile the root exactly** (complete prefix-free trie
  with paired children, plus Σ exact volumes = root volume);
* **no pending work:** the markers must be exactly the 256 depth-8 prefixes, each once; every leaf must lie in a
  marked subtree; marker and summary counts (including subtrees and max_depth) must match the records;
* every exclusion is re-derived by `ex_verify` (closed-form c_i, separate code), and every CERT is re-certified from
  scratch with the current `cert.jl` (logged numbers ignored). Both charts are required;
* **binding:** `runs/replay_chart<c>.jsonl` holds a header with the SHA-256 of every verifier source and of the input
  log, one result line per leaf with the SHA-256 of that exact log line, and a footer with the SHA-256 of the header
  plus all result lines, so every result is bound to its log record, to the input log and to the verifier sources;
* exit 0 only if everything passes. Any failure gives exit 1; test modes (`--partial`, `--cert-sample`) give exit 3,
  which is never a proof.

**Fail-closed test** (`replay_failclosed_test.sh` → `replay_failclosed_test.log`; corruptions by `corrupt_log.py`),
run on `runs/recert_chart2.jsonl` with the final `replay.jl`. Result: **PASS, 24/24 cases**. The intact log in test
mode exits 3. Exit 1, each with a specific FAIL line, for: unknown status; malformed line; extra field; duplicate
line; duplicate path; missing leaf; missing subtree marker (pending); a FAIL_FLOOR record; a wrong but canonical box
value; a non-canonical rational; a non-whitelisted rule; a false exclusion (a CERT leaf relabelled with a whitelisted
but false LEMMA_A rule, counts kept consistent); summary count mismatch; summary extra field; summary duplicate count
key; wrong summary subtree count; wrong summary max_depth; missing summary; mixed chart; trailing garbage; empty log;
a **forged parent certificate** (two certified sibling leaves replaced by one CERT record for their parent, with
tiling, marker and summary all kept consistent, so only re-running the certificate catches it); and a single-chart
run without `--partial`.

**Result of the fresh replay of the re-certified logs (the evidence):**
`./run_replay.sh` → `replay.log`, `runs/run_replay.out` (exit 0), 2026-10-01 16:43:37–16:53:30, 3 threads, nice 10.

| | chart I | chart II |
|---|---|---|
| log / lines / leaf records | `runs/recert_chart1.jsonl` / 45 426 / 45 169 | `runs/recert_chart2.jsonl` / 8 305 / 8 048 |
| phase 1 (schemas, box = path box, duplicates, 256 markers, strict summary, exact tiling) | ok | ok |
| exclusions re-derived (`ex_verify`) | 8 473 | 476 |
| certificates re-run (`certify`, repaired code) | 36 696 | 7 572 |
| re-verification failures / errors | 0 / 0 | 0 / 0 |
| worst margin R₂/((5/2)E), recomputed (upper bound) | 0.9999985269147025 | 0.9999679032558332 |
| phase 2 time | 455 s | 122 s |
| input log SHA-256 | `4c0a099b931c2978bbe3745577c59af9994aff93f24a404d658148a647767d04` | `c000783b12be553c84ad47f0c243b7066b5be7a88352d41b632c6b3d465670f1` |
| results digest (`results_sha256` in the footer) | `afc14a38c44ab008d7d36737499ac5817e07ec6a1b5390be93e874d72a86b65f` | `42d6fb96df52072c2a69a93e5bed4936fb6c02703b935574e80ec9d246507637` |

Verifier sources (SHA-256, recorded in both result headers and equal to the files in `MANIFEST.sha256`):
`replay.jl` 712ea04940d8faca64cb07808001c772f382b3c88fec7e07771b742c0d16277b,
`tiling.jl` b3300997d500b659240d97c525a80fd71bd3de4a38039dee7b54782869a24ad2,
`cert.jl` 6ba48ef36a0e4c11a227663137efc8b6a61c09e79e57b32e6dc105aecbb021b3,
`cert_float.jl` 2245450d6d217a7e026aa19eae733c0042e7eb3e3d118ca007ed9cac34924bfe,
`geom.jl` 551b37f94511f2c5264660841aa3f36db06476a1aea1fa851905f05e41f241af.
Each result file has one header, one `"verdict":"PASS"` line per leaf (45 169, resp. 8 048; no FAIL or NOT_REPLAYED
line), and a `"verdict":"PASSED"` footer.

**Superseded replays (history, NOT evidence).** (1) 10:01–10:55: a fail-open version of `replay.jl` on the original
logs with the defective `cert.jl`; its log is kept as `runs/replay_v1_failopen_superseded.log`. (2) 15:38–15:52: the
fail-closed `replay.jl` (before the strict summary schema), also with the defective `cert.jl`; it passed with the
same counts, and its console log and result files were overwritten by the fresh replay above. Both reproduced the
defective Gram enclosures and therefore do not count as evidence.

## 13. Caveats (LOUD)

1. **Review history:** the first certification was refuted as rigorous (Float64 nodal products, §0a). Lemma B′ is
   claimed only on the strength of the re-certified logs and the fresh replay with the repaired code.
2. **Inputs:** paper Lemma A (and through it the triangle theorem) and paper Lemma B are used as proved
   (`../paper/PAPER.md` §2.2–2.3); this directory does not re-prove them.
3. **Trusted computing base:** Arb/FLINT, GMP, Julia, and this code. The replay re-executes the same certificate
   implementation (`cert.jl`, `cert_float.jl`, `geom.jl: profiles`); an error in that code that makes a false
   certificate pass would be reproduced, not caught. The Float64-product defect was exactly such an error, so this
   caveat is not hypothetical. Mitigations: T7's exact-dyadic enclosure tests of the Gram assembly, the symbolic
   check of every formula, exact containment tests, the FEM bracket on 130 shapes, and the separately written
   exclusion re-derivation. The coordinator's separate checker (`../verify/`) is the independent check; it is not
   covered here, and its reported profile bug is in the coordinator's code (`verify/xband.py`).
4. **Multilinearity** (exact corner extrema) is verified by sympy for the closed forms. `geom.jl` computes the c_i
   from vertex coordinates by the cross-product definition, while the replay uses the closed forms; they agree on
   every one of the 8 949 re-derived exclusion leaves.
5. **Small slack:** the worst certified margin is 1 − 1.5·10⁻⁶ (chart I). This is rigorous (128-bit balls with
   relative radii ~10⁻³⁰ after the fix), but it leaves little tolerance for an error in the trusted components
   (item 3).
6. **FEM validation is not a proof**, only a sanity check of the bound values (§11 a).
7. **Chart roots contain non-valid points** (non-convex or with a distance > 1). Certificates are claimed only for
   convex Q(p), and Lemma 1 only uses valid points, so this is harmless. Degenerate faces (t = 0; x = 0 or 1) lie in
   closed boxes that are certified or excluded like any other; no limiting argument is used.
8. **Spec deviations:** chart I OUT_OF_BAND and LEMMA_B use the edge-width formula of Lemma 3 (§3). The Ritz domain
   is the inward-rounded Q_int ⊂ Q_in (Lemma 6(ii)). The 5/2 test is (5/2)E·D − N ≻ 0, which is equivalent to
   R₂ < (5/2)E.

## 14. Files and reproduction

Code: `geom.jl` (charts, exclusions, profiles), `cert_float.jl` (float layer), `cert.jl` (rigorous certificate,
repaired), `cover.jl` (cover driver), `recert.jl` (re-certification of a cover log), `tiling.jl`, `replay.jl`
(fail-closed replay), `run_cover.sh`, `run_recert.sh`, `run_replay.sh`, `replay_failclosed_test.sh` +
`corrupt_log.py`, `tests.jl`, `sympy_check.py`, `validate_fem.py`, `validate_bounds.jl`, `log_stats.py`. Float
prototypes from the spec: `proto.py`, `estimate.py`, `estimate2.py`.
Evidence: `runs/recert_chart{1,2}.jsonl` (+ `.out`, `runs/run_recert.out`), `runs/replay_chart{1,2}.jsonl`,
`replay.log`, `runs/run_replay.out`, `replay_failclosed_test.log`, `tests.log`, `sympy_check.log`,
`validation_fem.csv`, `validate_fem.log`, `validation_bounds.log`.
History (not evidence): `runs/cover_chart{1,2}.jsonl` (+ `.out`, `runs/run_cover.out`; produced by the defective
code), `runs/cert_v1_defective.jl.txt`, `runs/replay_v1_failopen_superseded.log`,
`runs/validation_bounds_v1_defective_code.log`.
Checksums: `MANIFEST.sha256` (`shasum -a 256 -c MANIFEST.sha256`).

```
~/.venvs/analysis-py313/bin/python sympy_check.py          # (b)
nice -n 10 julia -t 1 tests.jl 400                          # (c), incl. T7 exact-dyadic enclosure
~/.venvs/analysis-py313/bin/python validate_fem.py 60 && julia validate_bounds.jl   # (a)
./run_cover.sh                                              # original cover (both charts)
./run_recert.sh                                             # re-certify every leaf with the current cert.jl
./run_replay.sh                                             # fail-closed replay of the re-certified logs (exit 0 = PASSED)
./replay_failclosed_test.sh                                 # corrupted logs must all FAIL
```
