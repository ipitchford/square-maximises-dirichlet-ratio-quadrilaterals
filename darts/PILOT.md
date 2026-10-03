# Darts: pilot of the rigorous cover (2026-10-01)

All timings are wall-clock seconds on this machine (Apple silicon, 16 GB), single thread at `nice -n 10`, while the global2
Julia job ran with 8 threads (load average 50–120). "Box" = one leaf of the adaptive bisection.

## 1. Certificate cost per box (measured)
| chart (cert type) | mesh level k | CR dofs / P2 dofs | solve (eig + derivative vectors) | σ₂ (verified LDLᵀ) | rigorous TM part | total |
|---|---|---|---|---|---|---|
| BD (2 pieces)            | 4 | 736 / 961     | 0.03–0.05 s | 0.03–0.1 s | 0.1 s   | **0.2–0.4 s** |
| BD                       | 5 | 3008 / 3969   | 0.2–0.4 s   | 0.2–0.7 s  | 0.15 s  | **0.6–1.3 s** |
| BD                       | 6 | 12160 / 16129 | 1–2 s       | 1–3 s      | 0.3 s   | **3–6 s** |
| SPIKE2 (4 pieces)        | 5 | 6080 / 3969   | 0.6–1 s     | 0.5–1 s    | 0.4 s   | **1.5–2.5 s** |
| SPIKE (collar, 2 pieces) | 5 | 3040 / 1953   | 0.2 s       | 0.2 s      | 0.2 s   | **~0.6 s** |

(Order-3 Taylor models in 4 normalised parameters; second-order corrected CR test vectors; first-order corrected P2 vectors.)

## 2. What limits the box size (measured)
* Zeroth-order test vectors: Temple term ≈ 2500·r² (BD units) and growing with mesh refinement (the residual of a transported
  eigenvector is mesh-singular at the reflex corner) — useless beyond r ≈ 1%. First-order vectors: ε² ≈ 3.4 at 2% boxes, 62 at 4%;
  second-order vectors: ε² ≈ 0.0085 at 2%, 0.59 at 4%, 45 at 8% (moderate dart, k = 5). See `explore`-style probes in this folder.
* After that the bottleneck is the verified Temple shift σ₂ ≤ λ₂,h on the box: 775 of 782 failures in the first MAIN pilot were
  "σ₂ ≤ ρ". The global strain transfer (min over pieces) was replaced by one verified LDLᵀ of the piecewise Loewner minorant
  Σ_q μ_q K_q(p_c) − s Σ_q j_q M_q(p_c) per box (`gsigma.sigma_box`), which certifies boxes the old transfer rejected
  (e.g. σ₂ = 26.5 vs ρ ≤ 22.4 where the old bound gave 15.7).
* Typical certified leaf: relative half-width ≈ 3% of a_y (resp. |c_y|) in each coordinate (BD chart, k = 4).
* Darts with ξ close to 1 (thin, elongated) need tiny boxes because Temple needs λ₂,h − ρ > 0: the method stalls there (two
  pilot roots stalled at relative widths 1e-5). These are exactly the configurations for thin-domain lemmas (Lemmas 5.1, 5.2 cover
  tan b_A + tan b_C ≤ 0.1 and flat hulls with base angles ≤ 5.7°; the intermediate-thin band is not covered — see README).

## 3. Region volumes (Monte Carlo, `mc_regions.py`) and pilot extrapolation
Root = 10% relative box (log-polar units for A and C; ±10% in t for SPIKE).
| region | coordinates | volume | roots | pilot: boxes / root | pilot: s / root | predicted CPU time |
|---|---|---|---|---|---|---|
| MAIN  (|A|,|C| ≤ 3, tan b_A ≥ 0.3, b_A ≤ b_C) | (log|A|, b_A, log|C|, b_C) | 0.632 | 6 320 | 48 (old σ₂) / 4–30 (new) | 18 (old) / 1–14 (new) | 15–30 h |
| SHORT (0.02 ≤ tan b_A ≤ 0.3, a_x ≤ 1.5) | same | 0.295 | 2 952 | not piloted (thin A-triangles: expect ≥ MAIN) | — | ≥ 10 h |
| SPIKE (0.02 ≤ tan b_A ≤ 0.3, a_x ≥ 1.5, any |A|) | (log t, log|C|, b_C) | 16.2 | 8 098 | 1–15 | 2–20 | 10–40 h |
| COLLAR (tan b_A ≤ 0.02, any |A|) | (log|C|, b_C) | 12.1 | 1 206 | 1–45 | 0.2–14 | ~2 h |
| FAR / DEEP (a_y ≥ R or |c_y| ≥ R) | 3 resp. 2 parameters | — | — | charts implemented, not tuned (see §5) | — | — |
Total for the compact core alone (first estimate, old σ₂): ≈ 40–80 CPU-hours.
**Revised estimate from the production runs (per-box minorant σ₂, order-3 TMs; rates measured near the critical curves, where
boxes are smallest):** MAIN-CBOX: 799 boxes / 795 s for a log-polar volume ≈ 0.024 ⇒ ≈ 3.3·10⁴ s per unit volume ⇒ MAIN ≈ 6 h,
SHORT ≈ 3–6 h (thin A-triangles, assumed 1–2× MAIN density); SPIKE: 56 boxes / ≈ 100 s for volume ≈ 0.108 ⇒ ≈ 4 h for 16.2;
COLLAR: 4 boxes / 1.5–4 s for volume ≈ 0.04 ⇒ ≈ 0.3 h for 12.1. **Compact core ≈ 14–21 CPU-hours ≈ 7–11 h wall on the 2 permitted
threads**, uncertain by a factor ≈ 2 (box density varies with the margin). This excludes the regions that the present certificate
cannot reach at all (intermediate-thin darts, long arms; README §6), so the full cover is not ready to run.
**Decision: the full cover is NOT feasible within ~10 h of compute; it was not started.** Instead the critical neighbourhoods (where
ξ → 7/3, the only places where the margin is small) were certified completely (§4).

## 4. Runs actually performed (rigorous, logged, replayed) — final versions with dyadic roots
See README.md §2 for the exact parameter sets. Logs in `certs/` (an earlier non-dyadic set in `certs_old/` is superseded).
| run | chart | roots | certified boxes | discards | failures | wall time |
|---|---|---|---|---|---|---|
| crit_collar (C ∈ CBOX, tan b_A ≤ 0.02) | SPIKE @ t = 0.02, k = 5 | 4 | 4 | 0 | 0 | 1.5 s |
| crit_spike_X1 (C ∈ CBOX, 0.02 ≤ tan b_A ≤ 0.3, a_x ≥ 1) | SPIKE2 (X₀ = 1), k = 5 | 52 | 56 | 0 | 0 | 58 s |
| main_cbox (C ∈ CBOX, tan b_A ≥ 0.3, a_y ≤ 1.5479) | BD, k = 4–5 | 75 | 606 | 69 | 0 | 196 s |
| shallow_crit (BD box around the shallow critical curve) | BD, k = 4–5 | 256 | 5258 | 377 | 0 | 1666 s (27.8 min) |
In the shallow run 83% of the time (1380 s) went to the 64 roots with the thinnest A-triangles (a_y ∈ [0.085, 0.263]): box
widths scale with a_y in the BD chart (LEMMAS §7), as predicted by the pilot.
Replay times: collar 2.6 s, spike 80 s, MAIN-CBOX 155 s, shallow 1555 s (single thread). Total compute of this session (pilots, probes, FEM
surveys, production runs, replays): ≈ 3.5 CPU-hours.

## 5. Far arms (Lemma 9.3) — point tests of the implemented FAR chart (k = 4, C = equilateral apex)
Ratio bound λ₂(dart at a_y = R)/min(μ₁(core), tail) at single points (no box), u = cot b_A:
L = 1.5, R = 6: u = −1.5: 1.74, −0.5: 1.92, 0: 2.06, 0.5: 1.78, 1.5: 1.43, 3.0: 2.38. FEM truth for these darts: ξ ≈ 1.3–2.0.
The loss at u = 3 comes from the coarse P2 mesh on the long arm (λ₂ 9% too high at k = 4) and from the strip core; fix (not done):
bound λ₂ by the dart with the arm truncated (Dirichlet) at moderate length, and refine. DEEP chart (both arms cut) certified the
symmetric test box u_A = u_C = 0.

## 6. Second session (2026-10-01 afternoon): the full cover
**Why the FEM certificate alone cannot finish.** (Measured.) The Crouzeix–Raviart lower bound on the affine reference meshes is
accurate only when every triangle of the chart is fat and resolved: for C = (−0.70, −0.40) (CBD with angles 150°/13°/17°) the
discrete λ₁,h of the collar chart is 130.7 (k = 4), 201.5 (k = 5), 244.3 (k = 6) against λ₁ ≈ 260; even CBD alone gives 160.5 (k = 4).
The certificate therefore stalls on thin or obtuse triangles, long arms (Liu's h grows with the arm) and thin spikes.
**The 1-D (Born–Oppenheimer) certificate** (`bo1d.py`, `bo2.py`, LEMMAS §11): lower bound for λ₁ of a whole family from the
longest slice components of the union (outer domain) plus per-dart chord bounds and exponential tails for unbounded arms; upper
bound for λ₂ from P2 on an inner dart (`p2fix.py`). Cost 0.05–0.3 s per box. Float survey (500 survey darts, best of 9 axes): the
1-D bound is 9% below λ₁ in the median (min 64%); λ₂/E < 2.3 for ξ ≤ 1.4 (all), 96% for ξ ≤ 1.8, 74% for ξ ∈ [1.8, 2.0), rarely
above 2.0 — so the hybrid cover tries BO2 first and the FEM certificate for the fat, high-ξ boxes.
Examples: the stalled collar box at C = (−0.70, −0.40) and the stalled spike box at C = (−0.037, −0.167) certify in 0.05 s; the
pilot MAIN root that needed 483 FEM boxes (183 s) is ONE BO box (0.26 s); deep notch 30°/30° (both arms ≥ 20): ratio 2.28.
**Far arms**: inner dart at ρ = R; λ₂ of the truncated (tapering) arm converges slowly to the strip threshold, so R = 20 was chosen
(at the hardest probed V, (b_A, b_C) = (28.6°, 48.7°): ratio 2.03 at R = 7, 1.63 at R = 20 for a point; boxes ≤ 0.03 rad).
**Machine load.** The afternoon runs shared the Mac with ~45 foreign python3.14 processes (load average 90–500); wall times below
are therefore 2–10× the CPU times.

## 7. Third session (2026-10-01 evening): why PP stalled, and the fixes (measured on known-hard roots, then the pool)
**Diagnosis of the stall** (old code, PP_0_4_0_8/9/10: S-shaped short darts, A short and just left of D, C beyond D). The box FEM
certificate failed on σ₂ (LDLᵀ inertia of the Loewner minorant) although the true ξ ≈ 1.7: the strain factors μ_q of the thin
pieces (≈ 0.6–0.86 at boxes of 0.016 rad) cut the minorant's λ₂ by 30–40% (e.g. λ₂,h = 190 → minorant 115 at
r = (0.0156, 0.0156, 0.0156, 0.0312)); lower Temple shifts (down to 0.08 of the gap) did not help. Boxes had to shrink to
(0.0078, 0.0078, 0.0078, 0.0156) (thousands per root; PP_0_4_0_8 unfinished after > 45 min, 4,417 FEM boxes).
**INC certificate** (Lemma 11.4, `inc.py`): outer dart (apex angles at B and D maximal over the box) ⊇ every dart ⊇ inner dart;
λ₁ ≥ L₁(outer dart) is a POINT FEM certificate (no strain), λ₂ ≤ P2(inner dart). Ratio λ₂(Q⁻)/λ₁(Q⁺) on hard boxes: 1.95–1.99 at
r = (0.0078, 0.0156, 0.0078, 0.0156), 2.2–2.27 at twice that, 2.8–3.2 at four times; sensitivity per coordinate (doubling from
(0.0156, 0.0312, 0.0156, 0.0312)): b_A +0.17–0.25, l_A +0.04–0.10, b_C +0.28–0.43, l_C +0.01–0.02 (the split rule uses these
contributions). The CR mesh level of the outer dart matters for thin darts with large λ: at one near-diagonal point
λ₁,h = 216 / 292 / 323 at k = 4 / 5 / 6 (true ξ ≈ 1.5), L₁ of an outer dart 96.7 / 195.4 / 260.5 (cost 0.3 / 0.5 / 5 s), so the
level is raised adaptively when the float ratio is mesh-limited.
**BO2 heuristics (the larger effect).** (i) The P2 upper bound at k = 4 overestimates λ₂ of long darts by up to 38% (inner dart with
C at 12: λ₂^{P2} = 485.5 / 412.3 / 363.8 / 352.9 at k = 4 / 5 / 6 / 7), and the old loop stopped after k = 4 when λ₂^{P2} > 2.6E;
now k = 4, 5, 6 with a refinement-aware stop. (ii) xi_bo was taken from the LAST inner dart tried (often the C-only triangle,
λ₂ 3–6× larger) and, when BO2 was skipped, the stale value of a much larger ancestor was inherited by all descendants (FP_1_0_1_6:
xi_bo = 24 inherited, BO2 never retried, 3,193 FEM boxes); now the minimum over inner darts and mesh levels is used and a skipped
BO2 is retried after two halvings. (iii) Near misses retry with 2K slab cells and a local direction search (PP_1_10_1_7, long A at
21°: E = 33.5 → 36.3, enough for λ₂^{P2} = 85.0). (iv) For long arms the BO split scale for angles is sin b/ρ (the transverse
spread ρ·Δb of the hull against the arm thickness), and the BD-chart FEM (useless beyond ρ ≈ 4.5 on these meshes) is not tried.
(v) The mirror predicate is exact for dyadic boxes touching the diagonal b_A = b_C (the union-based enclosure rounded outward and
left such boxes undecidable: LOUD failures in PP_0_4_0_9/12 tests).
**Timings on the known-hard roots** (6 concurrent processes, nice 10): PP_0_4_0_8 152.8 s (old: > 45 min unfinished),
PP_0_4_0_9 458.2 s (old: > 1.5 h unfinished), PP_0_4_0_12 598.5 s (old heuristics: LOUD failures), FP_1_0_1_6 6.3 s (old: 2,548 s,
8,796 boxes), all with FAIL 0.
**The pool** (`pool.py`, `pool_launch.sh`, `pool_status.py`): 6 workers pull roots of PP, SP_spk, SP_col, PF and FP from one fixed
pseudo-random order (claims by atomic link), each worker with its own resumable logs; the running mean cost per root is then an
unbiased estimate for the remaining roots.
**Profiling the first pool hour (coordinator's request) and the cost benchmark.** Completed pool roots (code of 18:34–19:26):
BO2 took 93–97% of the time; FP's 14,426 certificates came from two roots computed by the pre-fix code; the inherited
"skip BO2 if the parent's xi_bo > 3.5" rule left 1,070 boxes of PP_0_4_1_13 without any attempt (its certificates have
λ₂^{P2}/E = 1.43). Where the margin is thin, the 1-D bound E is 20–30% below λ₁ (true ξ vs certified λ₂^{P2}/E:
SP_col_0_0_7_8 1.90 vs 2.46; SP_spk_6_0_3_9 1.80 vs 2.45; PP long A ≈ 2.0 vs 2.46), while the box FEM certifies boxes 2–4× larger
per side there (strain-limited beyond). Further cover changes (heuristics only): float screening of 49 slicing directions and the
far-arm variants (`bo2f.py`, 2 ms per evaluation vs 30–60 ms rigorous; best direction identical within 1% on 30/30 test cases;
the rigorous bound is computed for the best three), a fast path re-using the previous direction, BO2 skipped only where INC is in
charge, FEM tried after a failed BO2 and first after a FEM success (gated by its recent record), the exact ρ-monotone no-dart
predicate (accepts all 12,884 earlier discards and 5,509 skipped cells).
`bench.py` (the first N roots of each family in the pool's pseudo-random order, an unbiased sample), seconds per root:
| family (roots) | bench1 (19:34) | bench3 (20:00) |
|---|---|---|
| PP (10) | 54.0 | 47.4 |
| SP_spk (8) | 95.5 | 41.3 |
| SP_col (6) | 125.2 | 69.2 |
| FP (6) | 50.2 | 40.0 |
| PF (6) | 6.7 | 5.4 |
All FAIL 0. Prediction at 20:07: 3,691 roots ≈ 43 CPU-h ≈ 7.3 h wall at 6 workers (±30%: heavy tails, e.g. PP_1_10_1_7 364 s);
final replays ≈ 11 CPU-h.
**Heavy roots, work splitting and two fixes (21:00–21:50).** The heavy PP roots (b_A ∈ [14°, 29°]) have ξ_true ≈ 1.75–2.05
(the coordinator's P2-FEM values and ours), so their small boxes are method losses: (i) INC on PP_1_7_5_8 at its certified
boxes (0.0039, 0.031–0.0625, 0.031–0.0625, 0.031–0.0625): outer-dart construction 9–12% on λ₁, inner-dart 10–12% on λ₂,
point-FEM discretisation 2.5–3.3% (k = 5), P2 ≈ 0 — the construction loss sums over the four coordinates because λ₁ and λ₂ are
taken at opposite corners; (ii) the box FEM, which compares λ₁ and λ₂ pointwise, certifies (0.0156, 0.031, 0.031, 0.031) in the
PP_1_6_2_8 region (ξ ≈ 2.04) where INC needs 2–4× narrower b_A, b_C: PP order is now BO2, FEM, INC (37.7% vs 22.7% of that root
in 10 min); (iii) long arms (ρ ≥ 4.5): the FARL chart of the far families on the superset family 'arm ρ ≥ ρ_min' certifies 2×
(and after fix (b) 16×) larger boxes than BO2 per side. Fixes: (a) a worker re-taking an item it had itself handed off skipped it
(tiling gap caught by replay3 in the test); (b) `tm.TM.recip`: python-flint's `ball**3` is NaN for a ball containing 0, which made
every FARL FEM attempt with a parameter-independent piece fail (σ₂ = NaN); the tail bound now uses the exact upper end of q.
Work splitting (`pool.py`, `cover5.process_item`): items run ≤ 5 min, then hand their pending sub-boxes to `queue/`; tested with
3 s / 5 s limits including `kill -9` of a worker mid-item (replay3: 148/148 and 164/164 leaves, failures 0).
**Long arms: never split ln ρ (22:40).** Every method treats an arm with ρ ≥ 4 on the box as the far arm ρ ≥ ρ_min (BO2 far
variant, FARL chart; the inner dart uses ρ_min, and λ₂ barely depends on the length of a long arm), so splitting its ln ρ
coordinate only multiplies boxes. Measured on PP_1_10_4_10 (ρ_A 12–20, ρ_C 2.7–4.5), same code otherwise, 6 min each:
12.7% of the root covered with the old split scale vs 87.5% with ln ρ_A never split. Applied to PP and to the polar arm of PF/FP.
Also: the spike/collar families had no 'not a dart' discard; boxes with C tiny near B (convex quadrilaterals, no darts) were
uncertifiable by the family-uniform spike certificate (λ₂ of the shortest spike vs inf over a ≥ 1 of λ₁ exceeds 2.5 although every
dart has ξ ≈ 1.4): LOUD failure in SP_spk_5_0_11_1, fixed by the predicate of LEMMAS 13.1(a) for spikes.

## 8. Adversarial review of 2026-10-02 (REFUTED) — repairs, gates, measurements
**Pool outcome before the review.** All 3,891 family roots completed (FAIL 0) at 05:16; the replay of that package was stopped at
05:57 because the review found production-used rounding defects (README §0, §9). Those logs (certs3/) are evidence for the
re-certification only.
**Verification of the review's claims (before acting).** The 81 slivers were reproduced exactly (79 SP + 2 SF columns with both
neighbours skipped; widths 2^-59 at t = 0.02 and 2^-56 at t = 0.10125 and 0.2278125); every witness path was reproduced
(test_containment.py prints the old logged values next to the exact ones: e.g. the FARL tails of pp:80526, fp:325, fp:321 were
logged as 70.4266, 76.4947, 73.3508 (subdivision maxima, nearest-rounded) and are now 70.58074855293309, 76.51867366619814,
73.39535836531002, each ≤ the exact value; the closed-form convex maximum is tighter as well as rigorous).
**Gates (coordinator's rule).** (1) `rigor_lint.py`: AST lint of the 37 transitive decision sources; every float()/fin()/round()/
.mid(), ball power, bare max()/min() and nextafter is either removed or annotated with a tag and a reason (292 sites; `--list`).
(2) `test_containment.py`: 60,362 strict containment checks against exact rationals / 300-bit references (no tolerance).
(3) Independent review of exactly the code pinned by `review/round2_files.sha256` (diff: `review/round2.diff`).
**Measurements with the repaired verifier (`recert.py --time FAM 40`, 840 leaves, 0 failures), seconds per leaf:** PP BO2 0.061,
INC 0.193, FEM 0.194, FARL 0.198; SP_col BO2 0.064, FEM 0.123; SP_spk BO2 0.050, FEM 0.280; FP BO2 0.075, FARL 0.200; PF BO2
0.056; FF/SF BO2 0.04–0.05; discards < 0.001. Legacy replay2 with the repaired verifier: main_cbox 675/675 (99 s), kf_collar
218/218 (88 s), crit 60/60, kf_spike shard 569/569 (425 s), shallow shard ok; all 1,042 legacy discards re-verify (exact decimal
thresholds). Pipeline tests: FF recert → merge → replay3 540/540; forced subdivision of two leaves → replay 540/540.
Prediction: re-certification 10.4 CPU-h (+ subdivisions, expected < 1 CPU-h) ≈ 1.5–2.1 h wall at 7 workers; fresh bound replay
≈ 11.7 CPU-h ≈ 1.8 h; total ≈ 3.3 h (range 3–4.5 h).

**Round 3 (review round 2: BLOCK on the replay evidence chain; numerics and Lemma 9.6 cleared).** Every SERIOUS 2-7 claim was
reproduced before acting (wrong-family replay of crit_collar as SPIKE_X1 exited 0; SHARDS / QUICK unvalidated; summary-only
replay_ok; log re-read for its hash; no FP or runtime binary binding; no aggregate receipt; legacy roots without exact ends:
4/4 crit_collar, 52/52 crit_spike_X1, 63/75 main_cbox, 256/256 shallow_crit, 378/882 kf_spike, all cover5 cells exact).
Repairs in this package's own code (README §9, round 3); SERIOUS 1 (external sealed replay wrappers) not acted on, by instruction.
Checks: rigor_lint 39 files, 451 annotated sites (169 in exempted scopes, each in rigor_allowlist.json), 0 failures;
test_containment 60,594 checks, 0 failures (now with exact Lemma 9.6 H1-H3 over FARL chart boxes); test_replay_infra 23 probes
pass (SHARDS 0/-1/abc/07/65 and QUICK 0/false/yes rejected, only proof mode can PASS, wrong family, synthetic / dropped / forged
receipts, log changed mid-replay, FE_UPWARD / FE_DOWNWARD / FE_TOWARDZERO / FZ, hash-chained receipt of a machinery run);
check_charts, accounting, lemma tables, Corollary H', CBOX cover pass.  Overhead of the new binding: about 1.3 s per replay
process (2,703 runtime files hashed at start and end); the FP guard per leaf is a few microseconds.  Time estimate unchanged.

**Round 4 (review round 3: harness / provenance; launch path before the queued run), measurements 2026-10-02.**
Memory (/usr/bin/time -l, peak RSS): one replay3 shard of certs3/pp.merged.jsonl (621,307 lines; every shard parses the whole
log) 1.50 GB, 355 s wall at shard 0/64 (parse + manifest + tiling of all 2,047 roots ≈ 170 s, then 3,958 leaves in 184 s); the
driver's validation of the PP log (log_leaf_hashes) 0.74 GB, after the shards have exited; one recert worker on three heavy PP
roots (PP_1_7_3_8, PP_1_10_3_8, PP_1_10_2_7: 2,148 leaves, BD k = 4..6, FARL k = 4, 5, INC; 0 subdivisions) 0.43 GB, plus resume
bookkeeping of a PP share (status-only since round 4; tracemalloc: 59 MB at W = 4, 73 MB at W = 3; it was ~3.2 kB per leaf =
0.25-0.33 GB with whole records) and the FARL problem cache (~0.02 GB): 0.52 GB per worker.  Sizing (coordinator's rules): recert
W = 3 (3 x 0.52 x 1.2 + 0.3 = 2.17 GB -> --memory-gb 2.5; W = 4 would need 2.8 GB); replay SHARDS = 3 (3 x 1.50 x 1.1 + 0.4 =
5.35 GB -> 5.5 GB at the certs3 line count; launch_replay.sh rescales by the certs4 PP line count).
Recert samples with the round-4 code (`recert.py --time FAM 300 11`, 4,800 leaves): 0 failures in every (family, method) except
SP_col FEM, 27/300 = 9% ('FEM: cart box': the leaf's t-range exceeds the collar chart's t1 = float(0.02), refused by the repaired
collar FEM since round 1), about 290 leaves, which recert subdivides (earlier straddle test: certified by BO2 after one or two
splits). Seconds per leaf: PP BO2 0.061, FEM 0.204 (max 4.4), FARL 0.174, INC 0.188; SP_spk BO2 0.051, FEM 0.234; SP_col BO2
0.038, FEM 0.100; FP BO2 0.046, FARL 0.170; PF BO2 0.041; discards < 0.001.  INC: recert recomputes E (fresh apexes, fresh
inc.l1_point at the recorded CR level, fresh P2 from the recorded level up) -- the 50 old INC leaves of replay shard 0/64 whose
stored E exceeded the fresh bound all re-certify (P2 level unchanged).  Leaf CPU 9.6 h (sampled families) + FF / SF 0.02 h.
Predicted wall time on P-cores: recert at W = 3 ≈ 3.2 h of leaves + splits / starts ≈ 0.1 h + end tail (largest root 10,989
leaves ≈ 0.3 h) ≈ 3.4-3.7 h (W = 4: 2.6-2.9 h); replay at SHARDS = 3 ≈ (9.6 + 1.3 legacy) / 3 ≈ 3.6 h + gates ≈ 0.2 h + per-shard
parse ≈ 0.1 h ≈ 3.9 h (SHARDS = 4: ≈ 3.0 h).  Workers placed on E-cores run about 1.5-2.5x slower.
