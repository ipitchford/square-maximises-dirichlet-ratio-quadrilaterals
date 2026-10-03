# BAND_SPEC — Lemma B′: the thin band 1/10 < width/diam ≤ 2/5 (convex quadrilaterals)

Author line in every document: Anonymous.  Working directory: `quad-ppw-2026-09-30/bandB/`.

## 0. Why
The transport cover (`../global2/`) discards boxes on which s := width(Q)/diam(Q) ≤ 2/5 holds everywhere (status `E3T`,
also inside `COMPOUND` trees) and relies on a lemma that is NOT yet proved.  Lemma B (paper) proves ξ < 5/2 for s ≤ 1/10.
This task proves

> **Lemma B′.** Every convex quadrilateral Q with 1/10 < s(Q) ≤ 2/5 satisfies λ₂(Q) < (5/2) λ₁(Q)  (Dirichlet).

by a rigorous, replayable computer cover of two 4-parameter charts, using two one-dimensional (Born–Oppenheimer)
bounds along the diameter.  Float feasibility study: `../explore/thin_REPORT.md` (worst 1-D ratio ≈ 2.23 for s ≤ 0.4).
Float prototype of exactly this design: `proto.py`, `estimate.py`, `estimate2.py` (Monte-Carlo estimate ≈ 2·10⁴
optimally-shaped boxes in total; expect 10⁵ with bisection overhead).  The prototype is NOT rigorous (dense float
eigensolves, Gauss quadrature, no minorant) — reuse ideas, not code.

## 1. Charts (prove each statement in BAND_README.md)
Let P, P′ be a diameter pair of vertices; normalise by a similarity so P = (0,0), P′ = (1,0).  All six vertex
distances are ≤ 1.

**Chart I (P, P′ opposite).**  P0 = (0,0), P2 = (1,0), P1 = (x1, −y1), P3 = (x3, y3), y1, y3 > 0 (convexity puts the
other two vertices on opposite sides of the diagonal; conversely every such configuration with x1, x3 ∈ (0,1) is
convex because the diagonals cross).  WLOG y1 ≤ y3 (reflect y ↦ −y) and x3 ≤ 1/2 (reflect x ↦ 1−x, swapping P0, P2).
Facts: 1−√(1−y²) ≤ x ≤ √(1−y²) for both other vertices; s ≥ y3 (Q ⊇ triangle P0P2P3, whose longest side is P0P2 = 1,
so its minimal altitude is y3; width is monotone under inclusion); s ≤ y1 + y3 (vertical extent).  Hence in-band ⇒
y3 ∈ (1/20, 2/5].  Coordinates: (y3, t, x1, x3) with t = y1/y3 ∈ (0,1].
Cross products (paper notation c_i = (v_{i+1}−v_i)×(v_{i+2}−v_i), order P0,P1,P2,P3 is counter-clockwise):
c0 = y1, c2 = y3, c1 = (P2−P1)×(P3−P1), c3 = (P0−P3)×(P1−P3).

**Chart II (P, P′ adjacent).**  P0 = (0,0), P1 = (1,0), P2 = (x2, y2), P3 = (x3, y3), y2, y3 > 0, counter-clockwise,
convex (all c_i > 0).  Then 0 < x3 < x2 < 1 (upper chain of a convex polygon between its leftmost and rightmost points
is x-monotone; the diameter constraint gives x3 ≥ 1−√(1−y3²) > 0, x2 ≤ √(1−y2²) < 1).  WLOG y3 ≤ y2 (reflect
x ↦ 1−x, relabel P0↔P1, P2↔P3).  Then **s = y2 exactly** (width ≥ minimal altitude y2 of triangle P0P1P2, whose
longest side is P0P1 = 1; width ≤ vertical extent y2).  In-band ⇒ y2 ∈ (1/10, 2/5].  Coordinates (y2, t, x2, x3),
t = y3/y2 ∈ (0,1].  c0 = y2, c3 = y3, c1 = y2(1−x3) − y3(1−x2), c2 = x2·y3 − x3·y2.

## 2. Exclusions (decided exactly in rationals, or by rigorous intervals, on the WHOLE box)
* `LEMMA_A`: c_i ≤ (7/200)·c_{i+2} for some i at every point of the box (paper Lemma A ⇒ ξ ≤ 2.4995).
* `LEMMA_B`: an upper bound of s on the box is ≤ 1/10 (chart I: s ≤ y1+y3; chart II: s = y2) ⇒ paper Lemma B.
* `OUT_OF_BAND`: a lower bound of s on the box is > 2/5 (chart I: y3; chart II: y2).
* `OUTSIDE_CHART`: on the whole box some vertex distance is > 1, or (chart II) some c_i ≤ 0 / x3 ≥ x2, or (chart I) a
  WLOG inequality fails everywhere.  Points outside the chart are represented elsewhere by §1.
Boxes that straddle an exclusion boundary are certified (the certificate holds for every Q in the box) or bisected.

## 3. The box certificate
For a box B, construct (polylines on [0,1], kinks known exactly):
* an OUTER slice length h_out(x) ≥ h_Q(x) for every Q ∈ B (h_Q = length of Q ∩ {x = const}):
  chart I: h_out = y1⁺·T⁺(x; x1⁻,x1⁺) + y3⁺·T⁺(x; x3⁻,x3⁺), T⁺(x;l,u) = min(1, x/l, (1−x)/(1−u));
  chart II: h_out = upper concave envelope of (0,0), (1,0) and the 8 corners of the P2- and P3-boxes;
* an INNER domain Q_in = {0<x<1, a_in(x) < y < b_in(x)} ⊂ Q for every Q ∈ B:
  chart I: a_in = −y1⁻·T⁻(x; x1⁻,x1⁺), b_in = y3⁻·T⁻(x; x3⁻,x3⁺), T⁻(x;l,u) = min(x/u, (1−x)/(1−l));
  chart II: a_in = 0, b_in = upper concave envelope of (0,0), (1,0), apex₂, apex₃, apex_k = (u_k, y_k⁻)/(1+u_k−l_k)
  with [l_k,u_k] the x-range of P_k  (∩_B triangle(P0,P1,P_k) is that triangle; ∩_B Q is convex and contains both).
  (y1⁻ etc. are rigorous lower/upper bounds of the coordinates over the box.)

**Lower bound (prove; cf. header of `../global2/bound1d.jl`).**  For convex Q and v(x) = (∫φ(x,y)²dy)^{1/2}:
λ₁(Q) ≥ μ₁(−d²/dx² + V, Dirichlet on [0,1]) for every 0 ≤ V ≤ π²/h_Q².  Take V = π²/H_k² on cells J_k with
H_k ≥ max_{J_k} h_out (h_out is linear between kinks, so the max is at a cell end if kinks are nodes).  Certify
μ₁ > E by Sturm oscillation: the solution of −v″ + (V − E)v = 0, v(0)=0, v′(0)=1 has no zero in (0,1]
(exact 2×2 transfer matrices per cell in ball arithmetic; hyperbolic cells: positivity at both cell ends suffices;
oscillatory cells: require k·|J_k| < π and positivity at both ends).  Use a mesh graded geometrically towards x = 0
and x = 1 (h ~ x there, so the minorant needs cells of relative size ~1 %), kinks as nodes; choose E by float
bisection then verify, e.g. E = (1 − 10⁻⁶)·(float μ₁ of the minorant).  `no_zero` in `../global2/bound1d.jl` is a
working ball-arithmetic implementation of this check (Arblib.jl); reuse it if you work in Julia.

**Upper bound.**  For u = f(x)·sin(π(y − a)/h), h = b − a (a = a_in, b = b_in), f(0) = f(1) = 0:
  ∫|∇u|² = N(f) := ∫₀¹ [ f′²h/2 + f f′ h′/2 + f² (π²/h)(1/2 + a′b′/2 + h′²(1/6 + 1/(4π²))) ] dx,
  ∫u²   = D(f) := ∫₀¹ f² h/2 dx.
(Derivation: ∫₀¹sin²=1/2, ∫₀¹η sin(πη)·π cos(πη) dη = −1/4, ∫₀¹cos²(πη)=1/2, ∫₀¹η cos² = 1/4, ∫₀¹η² cos² = 1/6+1/(4π²);
note a′²/2 + a′h′/2 = a′b′/2.  Re-derive it symbolically yourself and keep the script.)  By min–max on the 2-D
space {f(x) sin(π(y−a)/h): f ∈ span(f1,f2)} ⊂ H¹₀(Q_in): λ₂(Q) ≤ λ₂(Q_in) ≤ R₂ := larger root of
det(N_ij − R·D_ij) = 0.  f1, f2 may be ANY functions (use float P1 eigenvectors of the discretised weighted problem on a
kink-aligned graded mesh); only N_ij, D_ij must be enclosed rigorously.  Per element f, a, b are linear, so every
integrand is a polynomial except f²/h, which has the closed form ∫(α+βt)²/t dt with t = h (log term) — use it in ball
arithmetic; on the two end elements h vanishes at 0 or 1 together with f, so f²/h is a polynomial there (no log).
Check D ≻ 0 rigorously.  The coefficient (1/2 + a′b′/2 + …) can be negative on steep elements: exact integration
handles the sign; interval shortcuts must respect it.

**Certificate:** R₂ < (5/2)·E  ⇒  ξ(Q) < 5/2 for every Q ∈ B.

## 4. Cover
Adaptive bisection of each chart domain (dyadic rational box coordinates; bisect the coordinate with the largest
estimated loss, or the largest width relative to the prototype's sensitivity).  Roots: chart I
y3 ∈ [1/20, 2/5], t ∈ [0, 1], x1 ∈ [0, 1], x3 ∈ [0, 1/2]; chart II y2 ∈ [1/10, 2/5], t ∈ [0,1], x2, x3 ∈ [0, 1]
(x3 < x2 enforced by OUTSIDE_CHART).  Degenerate faces (t → 0, x → 0 or 1, x3 → x2) must end in exclusions
(LEMMA_A / OUTSIDE_CHART) or certificates — report any box that refines below a floor (e.g. width < 2⁻²⁴) LOUDLY as
a failure; do not drop it.

## 5. Logging, replay, validation
* JSONL, one record per leaf: chart, box (exact rationals as strings), status, and for CERT: E, R₂ (decimal strings,
  rounded outward), mesh sizes, seconds.  A summary line with counts and wall time.
* `replay.*`: re-derives every exclusion exactly and re-runs every certificate from the logged box (it must not trust
  logged numbers), and checks that the leaves tile both chart domains exactly (sum of exact volumes + no overlaps by
  tree structure).  Write a MANIFEST.sha256.
* Validation before the rigorous run (keep logs): (a) 1-D bounds vs 2-D FEM (`../explore/quadfem.py`, `xi(V,k)`; check
  μ₁ ≤ λ₁ and R₂ ≥ λ₂) at ≥ 50 random band shapes of both charts; (b) symbolic check of N, D; (c) unit tests of the
  inner/outer containments at random points of random boxes (exact rationals).
* Anything not rigorous: say so LOUDLY in BAND_README.md.  Do not claim more than was replayed.

## 6. Resources and conduct
* 16 GB Mac shared with other heavy work.  While `pgrep -f cover_t.jl` is non-empty use ≤ 2 threads at `nice -n 10`;
  afterwards ≤ 6 threads.  One heavy job at a time from you.
* Do not modify `../global2/`, `../local/`, `../darts/`, `../paper/`.  Python: `~/.venvs/analysis-py313/bin/python`
  (python-flint 0.9, numpy, scipy, sympy).  Julia with Arblib is available (`/opt/homebrew/bin/julia`, see ../global2).
* Deliverables in `bandB/`: code, `BAND_README.md` (status first; theorem; chart proofs; containment proofs; certificate
  lemma with proof; numbers; caveats), JSONL logs, replay script + its log, MANIFEST.sha256.
