# Global step: certified box cover for xi = lambda_2/lambda_1 < 5/2 on convex quadrilaterals (spec, 2026-10-01)

## Parametrisation (Endo–Osting; use exactly this)
p = (a,b,c,d). Vertices (counter-clockwise; p = 0 is the unit square [0,1]^2):
 v1 = (0,0) + ½(a+d−b, −a+d+c), v2 = (1,0) + ½(−a+d+b, −a−d−c), v3 = (1,1) + ½(−a−d−b, a−d+c), v4 = (0,1) + ½(a−d+b, a+d−c).
c_i(p) = (v_{i+1}−v_i)×(v_{i+2}−v_i) (indices mod 4; = twice the area of triangle v_i v_{i+1} v_{i+2}); convex non-degenerate iff all c_i > 0.
|Q_p| = 1 − a² − d². D4 acts by σh:(a,b,c,d)↦(a,b,−c,−d), σv:(a,b,c,d)↦(a,−b,c,−d), ρ:(a,b,c,d)↦(−a,−c,−b,d); xi is invariant.
Fundamental domain: F = {a ≥ 0, b ≥ 0, c ≥ 0} (d free).

## What must be covered
Every p in F with all c_i(p) > 0 that is NOT in one of the regions already handled by separate proofs:
 E1 local cube: |p|_∞ ≤ CUBE (CUBE = 0.039; read it from config.toml — it may become 0.059);
 E2 near-triangle (Lemma A'): c_i(p) ≤ 0.035·c_{i+2}(p) for some i;
 E3 thin (Lemma B): width(Q_p)/diam(Q_p) ≤ 0.1 (width = min over edges of the max vertex distance from that edge line);
 E4 non-convex: c_i(p) ≤ 0 for some i.
A box may be DISCARDED only if one of E1–E4 holds at every point of the box (certify this with interval arithmetic;
"one fixed condition holds on the whole box" is the intended test). Every other box must be CERTIFIED or bisected.
A certified box may contain points of E1–E4; that is fine (the certificate is valid wherever Q_p is convex).
Rigorous a-priori box (derived from |Q| ≥ width·diam/2, i.e. diam² ≤ 20|Q| when width/diam ≥ 0.1):
 a ∈ [0, 0.905], b ∈ [0, 4.36], c ∈ [0, 4.36], d ∈ [−0.905, 0.905]  (outside this box every convex p is in E3).
Monte Carlo says the uncovered region R has volume ≈ 0.29 and extent a ≤ 0.83, b ≤ 0.85, c ≤ 1.27, |d| ≤ 0.8; the
large a-priori box is cheap to prune with the E3/E4 tests.

## Certificate for a box B = p0 + Π[−r_i, r_i] (all four parts rigorous)
Test functions (fixed for the box, from a Float64 solve at p0): u_k = b(x;p)·P_k(x), k = 1,2,3, with bubble
b = Π_k ℓ_k, ℓ_k(x) = (v_{k+1}−v_k)×(x−v_k), and P_k fixed polynomials (Legendre products in coordinates scaled by the
bounding box of Q_{p0}, total degree ≤ n): coefficient vectors = the first three Rayleigh–Ritz eigenvectors at p0.
 (1) Upper side: the 3×3 pencil K = (∇u_i,∇u_j), M = (u_i,u_j). Rayleigh–Ritz: λ_2 ≤ θ_2 (pencil), M ≻ 0 required.
 (2) Lower side, Lehmann–Goerisch with n = 1 (NOT Kato–Temple — KT fails for angles ≳ 140°; see
     ../explore/lg_proto.py, which gives λ_1 within 0.1% at a 157° angle where KT is 40% low):
     flux w = ∇u_1/ρ0 + τ(p) + curl ψ, with ρ0 = θ_1(p0) (a fixed number), τ(p) = −(1/ρ0)·(∫_{x_ref}^{x}(Δu_1+ρ0u_1)(t,y)dt, 0)
     (x_ref a fixed number), ψ a fixed polynomial (least-squares minimiser of ‖w‖² at p0; any fixed ψ is admissible).
     Then div w + u_1 = 0 exactly for every p. With A0 = K11, A1 = M11, A2 = ‖w‖² and a shift β ≤ λ_2(Q_p) with
     β·A1 > A0:   λ_1 ≥ L := (β A1 − A0)/(β A2 − A1)      [Lehmann–Goerisch; β A2 − A1 > 0 follows from Cauchy–Schwarz].
     β = constant lower bound for λ_2 over the box: max of Krahn–Szegő 2π j01²/max_B|Q_p| (j01 ∈ [2.4048255576957727,
     2.4048255576957729]) and the bounding-rectangle bound π²·min(4/A²+1/B², 1/A²+4/B²) for rectangles containing all
     Q_p, p ∈ B (axis-aligned: A = 1 − a + max(|b|,|d|), B = 1 + a + max(|c|,|d|); also the four edge-aligned ones).
 (3) Criterion without divisions: Ã := K·(βA2 − A1) − (5/2)(βA1 − A0)·M = (βA2 − A1)(K − (5/2)L M).
     If Ã11 < 0 and the 2×2 matrix S' := −Ã11·Ã_cc + Ã_c1 Ã_1c (c = {2,3}) has a negative eigenvalue
     (det S' < 0 or tr S' < 0), then Ã has ≥ 2 negative eigenvalues (Haynsworth), so θ_2 < (5/2)L ≤ (5/2)λ_1 and xi < 5/2.
     Side conditions to certify: βA1 − A0 > 0, M ≻ 0 (M11 > 0, leading 2×2 minor > 0, det M > 0).
     Why three functions: near the square λ_2 and λ_3 are nearly degenerate in the isosceles-trapezoid plane (b,c)
     (gap ≈ 9|p|², margin 2.5λ_1−λ_2 ≈ 14|p|²), so a single fixed second test function loses accuracy across a box.
 (4) Taylor models: every base quantity (6 entries of K, 6 of M, A2) is a polynomial in p (vertices affine in p,
     quadrature nodes affine in the vertices, exact Gauss rules). Represent each as an order-2 Taylor model in δ = p − p0:
     exact (ball/interval) value, gradient, Hessian at p0 + a rigorous remainder over the box (order-3 and higher
     terms). Propagate through the criterion polynomials with Taylor-model arithmetic, then bound each criterion over
     the box (quadratic part over the box + remainder). Evaluate fixed polynomials at moving nodes by a finite Taylor
     expansion at the centre node position with exactly computed derivatives (as quadcert.jl's taylor_eval! does) —
     this avoids the ×10^2–10^4 dependency blow-up documented in ../global/README.md.
Quadrature: exactness argument as in ../global/README.md (Duffy rule on the two triangles of the diagonal v1v3; the
|w|² integrand has degree 2·max(n+5, kψ−1); the inner x-integral of a degree-(n+4) polynomial needs ⌈(n+5)/2⌉ points).

## Required workflow
 1. Float64 pilot of the whole certificate (no rigour): validate against FEM (../explore/quadfem.py, P2, k=6) at ≥ 20
    points incl. near the square (|p|∞ ∈ [0.04, 0.1]), near-triangles and width/diam ≈ 0.1–0.2; then estimate, by an
    adaptive Float64 bisection that mimics the rigorous test with a remainder model calibrated on ~20 rigorous boxes,
    the number of boxes and the total rigorous time. State the prediction with units BEFORE any long run.
 2. Rigorous implementation (Julia, Arblib or IntervalArithmetic; threads), adaptive bisection over F ∩ a-priori box,
    D4 fundamental domain only, JSONL certificate log (one line per certified/discarded box: box, reason, n, bounds),
    resumable, deterministic. Choose n per box (smaller far from the square).
 3. Run it to completion if the predicted time is ≤ 12 h on 8 threads; otherwise stop after the pilot and report.
 4. Independent replay script that re-verifies every logged box from the log alone (recomputing the certificate).
 5. README_global2.md: math, rigour argument, the exactness claims, timings, statistics, and any non-rigorous shortcut
    (forbidden in the certified run; if unavoidable, label it loudly).
