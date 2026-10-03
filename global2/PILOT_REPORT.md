# Pilot report — global step (Float64 pilot, calibrated on rigorous boxes)

Date 2026-10-01. Machine: 8 cores, 16 GB. Code: `pilot.jl`, `pilot_summary.jl`, `calib.jl`, `cover.jl --pilot`.

## A. Calibration of the rigorous cost (20 boxes, `runs/calib.txt`)

20 boxes that passed the Float64 pilot (all degree n = 8; the pilot never needed n > 8 on passing boxes) were run
through the full rigorous certificate (Arb, 128-bit balls, order-3 Taylor models):

* 20/20 certified (the Float64 test is the same Taylor-model computation; only rounding differs);
* time per certified box, one thread: **0.61 s** mean (0.52 s Arb + 0.09 s Float64), max 0.71 s.

## B. Whole domain (D4 fundamental domain ∩ a-priori box, cube = 0.049)

* Coarse level (rigorous exclusion tests only, refine to widths ≤ 0.3): **867 cells** survive.
* 25 cells sampled uniformly, each run through the adaptive Float64 pilot with a cap of 2 500 processed boxes:
  13 finished (mean 337 processed boxes, 143 certificate leaves), **12 hit the cap (48 %)**, 4 of them with zero
  certificates (refinement along junctions of different exclusion conditions).
* Throughput measured on 8 threads: 60 007 processed boxes in 977 s (`runs/shell025_pilot.log`), i.e.
  16.3 ms wall per processed box.

**Prediction (stated before any long run):**
* lower bound (capped cells counted at the cap): **≥ 1.19×10⁶ processed boxes, ≥ 1.03×10⁵ certificates,
  ≥ 5.4 h Float64 screening + ≥ 1.9 h Arb = ≥ 7.3 h on 8 threads**;
* realistic estimate: the capped cells do not close. In the near-flat-vertex layer (min c_i/c_{i+2} < 0.08,
  volume 0.049 = 16 % of R) the true loss of the certificate exceeds the margin at box radius 0.005
  (`scratch/t_e2box.jl`), giving ≈ 8×10⁷ boxes; junction cells refine to the minimum width. **Predicted ≥ 10⁷
  boxes, ≥ 45 h on 8 threads (≫ 12 h). The full run was therefore not started.**

## C. Sub-region: the shell 0.049 < |p|∞ ≤ 0.225 (feasible)

* Float64 pilot of the 0.25-cube (stopped at 60 007 processed boxes): its passing boxes cover 100 % of
  |p|∞ ≤ 0.225 (Monte Carlo, 4 000 points); 9 703 of them intersect that cube.
* **Prediction: ≈ 10⁴ certified boxes, ≈ 25 min on 8 threads.**
* Outcome (rigorous run, `runs/shell0225_cert.jsonl`): 8 216 certified boxes, 0 failures, 1 171 s (19.5 min).

## D. Transport cover (TRANSPORT_SPEC.md) — prediction stated before the run (2026-10-01 07:21)

Certificate: affine transport with first-order δ-corrections (fallback: plain transport), order-4 Taylor models,
pointwise Taylor-model eigenvalue bounds + Andrews–Clutterbuck bootstraps, both split diagonals, compound discard.
Region: F ∩ a-priori box minus the cube [0,1/5]³×[−1/5,1/5] (inside the certified shell), minus E2, E3 and the
DEFERRED thin band E3T = {width/diam ≤ 2/5} (see README §11 — not an established lemma).

Pilot history on 48 sampled coarse cells (the deciding numbers):
| certificate variant | coarse cells | mean certificates/cell (median) | capped | extrapolated certificates |
|---|---|---|---|---|
| plain transport, box-constant β/L0 | 5208 | 273 (78) | 6 | ≥ 1.4×10⁶ |
| + Taylor-model β(p), L0(p) | 5208 | 193 (72) | 4 | ≥ 1.0×10⁶ |
| + first-order δ-corrections, order 4, best variant | 5208 | 98 (22) | 3 | ≥ 5.1×10⁵ (98 % of them at width/diam < 0.35) |
| + thin band w/D ≤ 2/5 deferred, progress-guaranteed splits | 1737 | 26 (3) | 0 | 4.6×10⁴ |

Throughput: 4 248 processed boxes in 107.7 s wall on 8 threads (Float64 screening) = 25 ms per processed box.
Calibration (`runs/calib_t.txt`): 20/20 pilot-passing boxes certified in Arb; 0.20 s per certified box on one
thread (0.07 s Float64 choice + 0.13 s Arb).
**Prediction: ≈ 1.5×10⁵ processed boxes and ≈ 4.6×10⁴ certificates (±50 %: two heavy cells on the near-triangle
boundary dominate the mean); wall time ≈ 3 900 s screening + ≈ 1 000 s Arb ≈ 1.4 h on 8 threads (range 0.7–2.5 h).
Run started.**

**Outcome (2026-10-01).** 100 342 processed boxes, 31 426 certificates (0 FAIL/ERROR), wall ≈ 5 030 s + 163 s ≈ 1.44 h
on 8 threads — inside the predicted 0.7–2.5 h; boxes and certificates 33 % below the central prediction. The run was
interrupted by a machine crash/reboot at ≈ 08:45 (SIGTERM after 96 929 boxes) and resumed from its log at 08:52
(README §10.2). The machine was shared throughout (another agent's jobs, a video render, later the coordinator's
verification jobs); the job had between ≈ 1.5 and 4 of the 8 cores, so per-box times in the log are inflated
relative to the calibration (mean 0.67 s vs 0.20 s per certified box).
