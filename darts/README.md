# Darts (non-convex quadrilaterals): ξ = λ₂/λ₁ < 5/2 — status 2026-10-02 (Anonymous)

## 0. Status (honest)
**PROVED (2026-10-02, 18:43), computer-assisted.** Theorem 13.2 holds: every family leaf was re-certified with the repaired verifier (certs4/, 4 workers, 2.69 h, 0 FAIL; 232 SP_col FEM leaves subdivided), and a fresh, fully bound, fail-closed replay of all 14 logs (462,849 leaves; 8 gates; 2 shards per log; 5.63 h; pin 7e564f3f…) printed PASS. Receipt replays/receipts/20261002T130555_a1b1ffb8.jsonl, final digest 35d4a6fbdbbdafacb001746a5b3bcebdc4973f735b26b7f70cd6d6be90c66eee, file sha256 a2759b4fcd03361de493049622f5eb3aa2368259f8ad86d169b09bc6ae4c3a7a. Independently verified by ../verify/xreceipt.py (exact state machine, sha256 ce1d59bc…; reviewed pin; runtime anchors from the structure-only run 20261002T121100_1f0f004a) and covered by ../verify/xdart.py (160,000 samples, 0 uncovered). The dart certificates have no second, independent implementation (owner's decision). Records: ../final_run/. History of the package follows.

**History: NOT YET A PROOF (2026-10-02, morning).** An adversarial review (2026-10-02) refuted the previous package: nearest-float conversions in
rigorous decisions (FEM/P2 Schur test, FEM strain factors, atan / cot endpoints of the spike and FARL boxes, the FARL tail, the BO2
projection range), 81 uncovered t-slivers between skipped spike grid cells, an unguarded NaN path in the Taylor-model range, an
unbound replay (manifest, decision sources, environment, cached charts) and gaps in Lemmas 5.1, 5.2, 9.3, 9.4 and the missing
oblique-cut lemma (now Lemma 9.6). All of these are repaired in the sources (§9). The family logs (certs3/) were computed with the
defective verifier: every leaf must be re-certified with the repaired verifier and the whole package replayed fresh and bound
before Theorem 13.2 is claimed. **Previous replay results are void.** The repaired verifier passes the rigour lint
(`rigor_lint.py`), the containment tests (`test_containment.py`), the chart rebuild (`check_charts.py`) and the exact accounting
(`accounting.py`); the re-certification has not been run (it waits for an independent review of the patched code).
The second review (2026-10-02, round 2) cleared the numerics and Lemma 9.6 and blocked on the replay evidence chain (shard-count and
quick-mode validation, time-of-check snapshots, the FP environment, an aggregate receipt, family identity, unchecked replay reuse);
those findings are repaired in this package's own code (§9, round 3: `replay_all.py`, `replay_common.py`, `fpenv.py`, `binding.py`,
`replay_ok.py`, `rigor_lint.py`, `test_containment.py`, `test_replay_infra.py`). The third review (round 3) blocked on harness and
provenance only (no mathematical or certificate regression); its items and the launch-path defects found before the queued launch
are repaired in round 4 (§9). Launch method (decided by the user): long runs are queued through the machine's `research-run`
launcher (`launch_recert.sh`, `launch_replay.sh`); the proof replay is this package's own hardened driver (`replay_all.py`), not an
external replay wrapper. Threat model: §5.

Proof structure (LEMMAS.md §13, Theorem 13.2): after the reflection in BD (b_A ≤ b_C), every dart lies in a lemma region or in
the parameter region of one of the runs below; a run proves its region iff its fail-closed replay passes (all roots complete,
exact tiling, every discard re-decided, every certificate recomputed). The table in §3 lists every run with its replay result.

## 1. Conventions
**BD chart** (all runs, old and new). B = (0, 0), D = (1, 0) (D the reflex vertex, B the opposite vertex), A = (a_x, a_y) with
a_y > 0, C = (c_x, c_y) with c_y < 0; the dart is Q = int(tri(B,D,A) ∪ tri(B,D,C) ∪ BD), with D reflex iff
f_D = (1 − c_x)a_y + c_y(a_x − 1) < 0. Angles b_A = ∠ABD = atan2(a_y, a_x), b_C = ∠CBD = atan2(−c_y, c_x), in RADIANS, both in
(0, π); ρ_A = |A|, ρ_C = |C|; polar form A = ρ_A (cos b_A, sin b_A), C = ρ_C (cos b_C, −sin b_C); t = tan b_A.
θ_B = b_A + b_C; β = notch depth, 1 − β = (sin b_C/ρ_A + sin b_A/ρ_C)/sin θ_B (Lemma 13.1).
**Mirror**: the reflection in BD maps (a_x, a_y, c_x, c_y) to (c_x, −c_y, a_x, −a_y) and swaps b_A, b_C; all runs cover the domain
b_A ≤ b_C (boxes with b_C ≤ b_A everywhere are discarded as 'mirror'; the diagonal b_A = b_C is covered from the side
b_A ≤ b_C, boxes being closed).

## 2. Analytic results (complete proofs in LEMMAS.md)
| result | darts covered (ξ < 5/2) |
|---|---|
| Prop 5.1 | β ≤ 0.0339 |
| Lemma 3.1 (collapse collar) | dist(A, segment BD) ≤ 0.035·inradius(CBD) (and A ↔ C) |
| Lemma 5.1 | b_A, b_C acute and tan b_A + tan b_C ≤ 1/10 (ξ ≤ 2.4223) |
| Lemma 5.3 (Bessel zeros certified, `lemma_sector.py`) | θ_B ≤ 0.1772 rad (ξ < 2.4893; certified maximum 2.489242112…) |
| Lemma 5.2 | hull angles at A and C both ≤ arctan(1/10) = 5.7106° (in particular θ_B ≥ 174.29°) (ξ ≤ 2.4356) |
| Lemma 9.5 (`lemma_longarm.py`) | b_X ∈ [b_k, 180° − b_k], \|BX\| ≥ L_k (table in LEMMAS §9) |
| **Lemma 12.1 (thin stem)** | max over 0 ≤ x ≤ 1 of the vertical stem chord h_A(x) + h_C(x) ≤ 0.05 (ξ ≤ 2.4356) |
| Lemmas 4.1, 6.1, 9.1–9.4, 9.6, 11.1–11.4, 13.1 | tools: chord bound, tails, spike-head reductions, oblique cut of a long arm (FARL, 9.6), 1-D (Born–Oppenheimer) bounds with tails, outer/inner-dart inclusion certificate (11.4), polar dart condition |

External inputs: the triangle theorem ξ(T) ≤ 7/3 (Siudeja 2010; Arbon–Mannan–Psenka–Ragavan 2020; used by Prop 5.1 and Lemma 3.1)
and the Crouzeix–Raviart interpolation constant 0.1893 (Liu–Kikuchi 2010; used by FEM certificates only).

## 3. The runs (Theorem 13.2)
| run (log) | region (Theorem 13.2) | leaves (old logs) | status |
|---|---|---|---|
| `certs/crit_collar`, `certs/crit_spike_X1` | spike collar / spike, C ∈ CBOX | 4 / 56 | old FEM leaves; to be replayed with the repaired verifier |
| `certs/main_cbox` | Theorem M region | 675 | idem |
| `certs/shallow_crit` | Theorem H box | 5,635 | idem |
| `certs2/kf_collar`, `certs2/kf_spike` | C in the KF sector | 218 / 4,682 | idem |
| PP | set (7) of 13.2 | 306,017 (84,142 discards, 126,002 BO2, 51,382 INC, 44,491 FEM incl. 27,861 FARL) | to be re-certified (`recert.py`) |
| SP_col / SP_spk | sets (5) / (6), polar C | 41,885 / 62,671 | to be re-certified |
| PF / FP / FF | far arms | 8,545 / 29,339 / 540 | to be re-certified |
| SF_col / SF_spk | spike × far C | 393 / 996 | to be re-certified |

## 4. Logs: paths, coordinates, implicit constraints, statuses (for independent coverage checks)
Every log is JSONL. Record types: `leaf` (status `cert` or `discard`; `FAIL`/`defer` never occur in a log that passes replay),
`split`, `esc` (bookkeeping of the bisection), `root_done`. A box is p0 ± r (componentwise, closed); a coordinate with r_i = 0 is
inactive (its p0_i is a fixed value or a dummy 0, see below). A `cert` leaf proves ξ < 5/2 for EVERY configuration described by its
box (including non-darts, which is harmless). A `discard` leaf names the predicate that holds on the whole box.

**Old runs (first session; BD-chart Cartesian coordinates; method = FEM certificate, `cert` records carry no `method` key):**
| log | p0 / r coordinates | darts covered by a `cert` leaf | discard reasons (exact strings) |
|---|---|---|---|
| `certs/crit_collar.jsonl` | (t, 0, c_x, c_y), t = 0.02 fixed (r₀ = r₁ = 0) | every A = a(1, t′) with 0 < t′ ≤ 0.02 and ANY a > 0 (Lemmas 9.1–9.2), C in the (c_x, c_y) box | none |
| `certs/crit_spike_X1.jsonl` | (t, 0, c_x, c_y) | A = a(1, t′), t′ in the t-box, **a ≥ 1** (Lemma 9.4, X₀ = 1), C in the box | none |
| `certs/main_cbox.jsonl` | (a_x, a_y, c_x, c_y) | the box itself | `not a dart (f_D >= 0 on box)`; `tan b_A <= 0.3 (COLLAR/SPIKE certificates, C in CBOX)` |
| `certs/shallow_crit.jsonl` | (a_x, a_y, c_x, c_y) | the box itself | `not a dart (f_D >= 0 on box)`; `Prop 5.1` |

**KF runs (certs2/; same coordinates as the old collar/spike runs; FEM certificate, no `method` key):**
| log | p0 / r | covered by `cert` | discard reasons |
|---|---|---|---|
| `certs2/kf_collar.jsonl` | (t, 0, c_x, c_y), t = 0.02 fixed | A = a(1, t′), 0 < t′ ≤ 0.02, any a > 0; C in the box | `C outside sector: b_C < 38 deg`, `C outside sector: b_C > 100 deg`, `C outside sector: \|C\| < 0.5`, `C outside sector: \|C\| > 1.8` |
| `certs2/kf_spike.jsonl` | (t, 0, c_x, c_y), t ∈ [0.02, 0.3] | A = a(1, t′), t′ in the t-box, **a ≥ 1**; C in the box | same four sector reasons |

**Family runs (certs3/; cover5):** each family's proof log is `certs3/<fam>.merged.jsonl`, the concatenation (by `merge_pool.py`,
which fails if a (root, path) tree record, handoff or item_done occurs twice) of the static shard logs `<fam>.<i>.jsonl` and the
pool worker logs `<fam>.w<ID>.jsonl` (`pool.py`: dynamic load balancing, fixed pseudo-random root order, work items claimed via
`claims/<fam>/<item>`). Work splitting: an item (a root, or a sub-box of a root) runs for at most 5 min, then each pending sub-box
is recorded as `{"type": "handoff", "root", "path", "p0", "r", "depth", "k"}` and queued (`queue/<fam>/<path>.json`) for any
worker; every item ends with `{"type": "item_done", "root", "path", "ok", "counts", "handed_off"}`. A root's records may thus
be spread over several worker logs; a root is complete iff its leaves tile its box exactly (replay3 also requires every handoff to
have an item_done, all item_done ok, and the root item done; a `root_done` record is written when a root finished in one item). The
static PP shards `pp.0/1/2.jsonl` contain only completed roots: the records of their three unfinished roots (PP_0_4_0_8/9/10, all
under the old heuristics) were removed before the pool took over (originals in `pilot_runs/legacy_originals/`).
`cert` leaves have `method` = `BO2` (1-D lower bound + P2 upper bound on an inner dart, Lemma 11.3; optional key `arms` =
`farA`/`farC`/`farAC`: the named polar arm(s) replaced by the far arm ρ ≥ ρ_min for the lower bound, a superset family, LEMMAS
11.3 remark (4); the box's darts are still exactly those of p0 ± r), `INC` (family PP only,
Lemma 11.4: outer dart (`Aout`, `Cout`) containing and inner dart (`Ain`, `Cin`) contained in every dart of the box; `E` = rigorous
lower bound of λ₁ of the outer dart at CR level `k`; P2 test of the inner dart at level `kp`) or `FEM` (gcert certificate on the
Cartesian BD box `pc` ± `rc`, which contains every dart of the family box; FP/PF, and PP boxes with exactly one arm of ρ ≥ 4 on
the box and > 5 somewhere: the FARL chart, whose `pc`/`rc` are (u = cot b of the far arm, 0, Cartesian box of the other apex),
with the far arm ρ ≥ its smallest radius on the box). Angles in radians, ρ = exp(l).
| family (log) | p0 = (·,·,·,·) | darts of a box (implicit constraints) |
|---|---|---|
| PP (`certs3/pp.merged.jsonl`) | (b_A, ln ρ_A, b_C, ln ρ_C) | A = ρ_A(cos b_A, sin b_A), C = ρ_C(cos b_C, −sin b_C), all four coordinates in the box |
| SP_col (`certs3/sp_col.merged.jsonl`) | (t, 0, b_C, ln ρ_C) | A = a(1, t′) with t′ ∈ [max(t − r_t, 0), t + r_t], **any a > 0**; C polar as in PP (a FEM leaf, collar chart at t₁ = float(0.02), requires t + r_t ≤ t₁: Lemma 9.2 covers t ≤ t₁ only) |
| SP_spk (`certs3/sp_spk.merged.jsonl`) | (t, 0, b_C, ln ρ_C) | A = a(1, t′), t′ in the t-box, **a ≥ 1** (a_x ≥ 1); C polar |
| PF (`certs3/pf.merged.jsonl`) | (b_A, ln ρ_A, b_C, 0) | A polar; C = ρ(cos b_C, −sin b_C) with **ρ ≥ 20** (any length) |
| FP (`certs3/fp.merged.jsonl`) | (b_A, 0, b_C, ln ρ_C) | A = ρ(cos b_A, sin b_A) with **ρ ≥ 20**; C polar |
| FF (`certs3/ff.merged.jsonl`) | (b_A, 0, b_C, 0) | both arms with **ρ ≥ 20** |
| SF_col (`certs3/sf_col.merged.jsonl`) | (t, 0, b_C, 0) | A = a(1, t′), any a > 0; C with ρ_C ≥ 20 |
| SF_spk (`certs3/sf_spk.merged.jsonl`) | (t, 0, b_C, 0) | A = a(1, t′), a ≥ 1; C with ρ_C ≥ 20 |
Discard reasons in the family runs (exact strings; each re-decided in replay by `fam5`):
`theta_B >= pi on the box (not a dart)` · `mirror (b_C <= b_A on the box)` (older logs: `mirror (b_C < b_A on the box)`, same predicate) · `Lemma 5.2 (theta_B >= 174.29 deg)` ·
`Lemma 5.3 (theta_B <= 0.1772)` · `collar family (tan b_A < 0.02 on the box)` [delegated to SP_col / SF_col / kf_collar] ·
`spike family (tan b_A < 0.3 and a_x >= 1 on the box)` [delegated to SP_spk / SF_spk / kf_spike] ·
`inside the shallow box of Theorem H (certs/shallow_crit.jsonl)` · `C in CBOX, tan b_A >= 0.3, a_y <= 1.5479: Theorem M (certs/main_cbox.jsonl)` ·
`no dart with beta > 0.0339 (Prop 5.1 or not a dart)` · `Lemma 12.1 (vertical stem chords <= 0.05 on the box)` ·
`C inside the KF sector (certs2/kf_collar, certs2/kf_spike)` (SP families only; the KF sector is 38° ≤ b_C ≤ 100°, 0.5 ≤ ρ_C ≤ 1.8).
Lemma regions outside all grids (Theorem 13.2): ρ_A < e^{−2.5} (Lemmas 5.2/5.3 via 13.1(b)); ρ_C < e^{−4} (no dart unless
t < 0.02, then Lemma 12.1). Grids: see `run_fam.py` (GRIDS) and the manifests in `manifests/`.

## 5. Replay (fail closed) and accounting
**Threat model.** The replay machinery protects the proof receipt against ACCIDENTAL corruption (torn or partial writes, killed
processes, disk or copy errors), stale or changed files (sources, logs, manifests, chart caches, pin list, runtime files changed
before or during a run), concurrent agents editing files in this directory, and environment drift (interpreter flags, startup
hooks, PATH, numerics variables, FP environment, loaded libraries, niceness). Deliberate tampering by the owner of the machine
(who could rewrite the sources, the pin, the receipt and the checker together) is OUT OF SCOPE; the final digest of a receipt is
anchored outside this directory (HANDOFF, paper) so that a later change of a receipt is detectable.
`./replay_all.sh` (exit status 0 iff everything passes): certified tables of Lemmas 5.3 and 9.5, mesh conformity of the FEM charts,
Corollary H′, CBOX cover, `accounting.py` (KF quadtree tiling, family grids cover the claimed ranges, far/polar overlap), and
`replay2.py` (old and KF logs) / `replay3.py` (family logs): status whitelist, per-status schemas (no unknown keys), duplicate
(root, path) detection, completeness (no pending root), exact dyadic tiling, discards re-decided, certificates recomputed from
scratch; nonzero exit on any unknown, malformed, pending or failed record. Each checked record is written to `replays/*.jsonl`
with the sha256 of its log record, the family, and the binding sha256 (`binding.py`: the transitive closure of local imports of
the replay script, every cached chart file, the sha256 of the manifest bytes the replay parsed, the family, and the runtime
identity: interpreter binary and libpython, every installed file of python-flint / numpy / scipy / scikit-fem / mpmath, the NumPy
build and BLAS configuration, the FLINT version and arb context, the OS, the floating-point environment and the numerics
environment variables); the summary adds the sha256 of the sorted multiset of checked record hashes. A replay parses an immutable
in-memory snapshot of its log and manifest and records the sha256 of exactly those bytes; at the end it fails unless the files on
disk, the binding and the FP environment are unchanged, no production-only module was loaded and every local module came from this
directory (`replay_common.py`). The FP environment must be IEEE round-to-nearest with gradual underflow at start, after every leaf
and at the end (`fpenv.py`: run-time behavioural tests in Python and NumPy, and on arm64 the FPCR RMode / FZ / FIZ / AH bits).
Every manifest names its family; a replay under another family fails. Every leaf box, skipped cell and chart box must have exact
float ends (`nanguard.exact_box`); legacy / KF processed roots are represented by their exact dyadic snapping (LEMMAS §13).
**Aggregate driver** `replay_all.py` (wrapper `replay_all.sh`: `SHARDS` a canonical integer 1..64, `QUICK` unset or exactly 1,
`PIN` the reviewed pin list): verifies every pinned file and that every transitive decision source is pinned; runs the gates
(`lemma_sector`, `lemma_longarm`, `check_charts.py 8`, `rigor_lint.py`, `test_containment.py`, `hull_to_box`, `check_cbox_cover`,
`accounting`), the merges and, FRESH, exactly SHARDS replay processes per log into a new `replays/run_<ID>/` (nothing is reused);
validates EVERY receipt record of every shard (`replay_ok.validate`: keys, ok, family, binding, counts, leaves_sha256, and the
multiset of record hashes equal to the one recomputed from the bound log); re-checks the pin, the logs, the manifests and the
runtime identity at the end; and writes an append-only, hash-chained receipt `replays/receipts/<ID>.jsonl` (canonical records with
strict sequence, chained and self hashes, schema; write loop, fsync of file and directory; re-read and validated before success;
the printed FINAL DIGEST is the end record's self hash). Before any gate, merge or replay it checks, in seconds, that the
re-certification directory is complete (`certs_complete.py`: every manifest root exactly one root_done with ok = true, no FAIL leaf
-- all FAIL leaves are listed --, no torn, garbage or duplicate record, every handoff with its item_done, every claim of a worker
1..N and done) and aborts otherwise. Only a proof-mode run without any failure prints `PASS` and exits 0; quick mode prints
`STRUCTURE-ONLY (NOT A PROOF)`, test mode `MACHINERY TEST OK (NOT A PROOF)` (exit 3). The driver runs at the TARGET niceness 10
(no double nice under a launcher that already applies nice 10). `test_replay_infra.py` probes all of this (82 probes).
**Launch (queued, research-run).** `launch_recert.sh` queues `recert_all.sh 3` (`--cpu 3 --memory-gb 2.5`: measured peak 0.43 GB
per recert worker + resume bookkeeping and caches, sizing rule W x 0.52 x 1.2 + 0.3 <= 2.5 GB); recert_all.sh records N, the code
identity and the families at the first start and refuses any other N or code on resume, runs its workers isolated (-I -S -B,
bootstrap.py) under `env -i` at niceness 10, and exits with the verdict of `certs_complete.py`. `launch_replay.sh` (after the
round-4 review) sizes the replay from the certs4 PP line count (measured peak 1.50 GB per replay3 shard of the 621,307-line PP log;
SHARDS x 1.65 + 0.4 <= 6 GB: 3 shards, 5.5 GB) and queues `replay_all.sh` with the reviewed pin. `smoke_launch.sh recert|replay`:
one-minute dry runs of exactly this path under research-run (1 CPU; not yet run: executing research-run was denied to the
preparing agent by the permission system, so running them is the user's decision).
**Replay results: none valid at present** (the previous results were obtained with the defective verifier and are void; the fresh
bound replay runs after the re-certification).

## 6. Numerical validation (not part of the proof)
* `test_containment.py` (gate 2 of the 2026-10-02 repair): every bound function used by the certificates against exact rationals
  or 300-bit ball references, strict containment, no tolerance (helpers fdown/fup/amax/amin/enclose/exact_box/dec, dyadic snapping,
  Taylor-model ranges incl. zero-midpoint / huge / tiny radii, Schur upper bound, strain factors, atan / cot endpoints, FARL tail,
  BO2 projection extrema, spike tail, and the reviewer's witnesses pp:66327, pp:189234, pp:80526, sp_spk:7, fp:321, fp:325).
* `test_inc.py`, `validate_fem5.py`: float cross-checks of the INC and FEM certificates (not part of the proof).
* `recert.py --time FAM N`: recomputation of random leaves of each family / method with the repaired verifier (timing; no log).

## 7. What remains
1. The re-certification of every family leaf (`launch_recert.sh` -> `recert_all.sh 3`; failing leaves are subdivided and
   re-certified; about 9% of the SP_col FEM leaves, whose t-range exceeds the collar chart's t1 = 0.02, are expected to be
   subdivided, all other methods 0/300 in samples), then, after the round-4 review clears, a completely fresh, fully bound
   `launch_replay.sh` -> `replay_all.sh` whose receipt digest is anchored externally.
2. The old and KF logs are replayed (recomputed) with the repaired verifier; a failing old leaf would need a new certificate.
3. Unchanged from before: no part of the dart case is open in the decomposition of Theorem 13.2 once the runs replay.

## 9. Repairs after the adversarial review of 2026-10-02 (finding → fix)
* Schur test (FATAL 1): `schur.schur_upper` keeps ball endpoints through abs, max, product, division and comparison (gcert, p2up).
* Strain factors (FATAL 2): `gsigma.piece_factors` forms amc, b, λ_min in balls, μ rounded down, j up; σ₂ subtraction directed.
* Grid slivers (FATAL 3): every cover5 manifest cell has exact float ends (`run_fam.exact_cell`, outward; legacy / KF processed roots
  are represented by their exact dyadic snapping, LEMMAS §13); `remanifest.py` re-encoded the
  manifests (root snapping unchanged, every skipped predicate re-checked); `accounting.py` checks coverage on represented boxes.
* atan endpoints (FATAL 4): `fdown(atan)`, `fup(atan)` in the SP/SF arms and predicates.
* FARL (FATAL 5): `charts.cut_tail` (exact convex maximum, directed), `fam5.farl_setup` (directed cot, exact cut cells, directed R,
  exact-endpoint chart box); Lemma 9.6 states and proves the oblique cut with R > L > k₀u.
* Replay binding (SERIOUS): `binding.py`; strict manifest schema with exact p0 and r; leaf multiset hash; `check_charts.py`.
* BO2 projection extrema: `bo2.proj_range` (left down, right up); exact subdivision ends (`nanguard.cuts`); exact mid-angle.
* NaN path in `tm.TM.range`: branch when 0 is in the quadratic coefficient; every candidate checked finite.
* Every other conversion in a decision: directed through `nanguard.fdown` / `fup` only (`rigor_lint.py` enforces it); decimal
  thresholds exact (`nanguard.dec`); exact dyadic snapping (`grid.dyadic_root`); exact inverse of the reference edge matrix in `gcert.hmax_box`.
* Lemmas 5.1 / 5.2 (acute hypotheses), 9.4 (a_x = 1), 9.3 (exact hypotheses), 5.3 (tan 0.0997 > 0.1 handled; α₃₁ label).

**Round 3 (review round 2, 2026-10-02):** SHARDS / QUICK validated, only proof mode can PASS, exactly SHARDS successful and
validated outputs required per log (`replay_all.py`); immutable byte snapshots of log and manifest, end-of-run identity re-check
(`replay_common.py`); FP environment asserted and guarded per leaf (`fpenv.py`), interpreter / numerical packages / BLAS / FLINT /
OS bound (`binding.runtime_identity`); append-only hash-chained aggregate receipt binding pin list, driver, gates, merges, every
shard output and summary, runtime identity and the end-of-run unchanged state; family in every manifest, binding and receipt,
exact agreement required; no reuse in proof runs, and `replay_ok.validate` checks every receipt against the recomputed leaf
multiset; `rigor_lint.py` NumPy-taint rule (F5) and per-site allowlist of every site in an exempted function / module scope
(`rigor_allowlist.json`), production-only modules verified unloaded by every replay; `test_containment.py` exact Lemma 9.6
conditions H1–H3 over FARL chart boxes and a non-vacuity check for fp:321; `fam5.farl_setup` / `charts.farl_problem` reject
k₀ ≤ 0; legacy manifests carry their family; exactness wording qualified for legacy roots.

**Round 4 (review round 3: BLOCK on harness / provenance, no mathematical or certificate regression; threat model above).**
(1) Runtime inventory: a fresh inventory at every end check (no memoisation) of the interpreter, every loaded non-OS Mach-O image
(dyld: the Python framework binary, every extension module and bundled library), the standard library, the numerical packages and
every *.pth and top-level module of the site directories; module origins checked against it, against the pinned sources or against
the installers' RECORD hashes (late imports), no bytecode read (-B, empty pycache prefix), sys.flags / _xoptions / pycache prefix
bound, `python -O` refused (`binding.py`, `bootstrap.py`). (2) Chart caches: the exact refcache/* inventory and hashes are pinned
and verified before any unpickling, before every replay and at the end; a pinned run never writes a cache (`chart.py`); the
rebuild-and-compare gate stays (`check_charts.py`). (3) The pin file is read once and its bytes parsed and hashed (private
read-only copy for the children); every consumed log and manifest snapshot must equal its pin entry or the recorded merge output.
(4) Shard outputs created exclusively (O_EXCL), snapshotted once after the child exits; exactly those bytes are validated,
summarised and hashed; a missing summary fails the run. (5) Receipt: self hash per record, strict sequence and schema, write loop,
fsync of file and directory, re-read before success, final digest printed for external anchoring. (6) Bootstrap: absolute
package path and interpreter (never from $0), in-process niceness, children with a sanitised environment and -I -S -B
(`replay_all.sh`, `bootstrap.py`). (7) Lint aliases (from-numpy imports, ufunc reductions, np.all / any, .var / .std) and
comment-versus-string detection; Lemma 9.6 test also by analytic extrema over u-intervals (`farl_conditions_interval`); FIZ / AH /
FZ16 probes (N/A where the CPU lacks FEAT_AFP); receipt metadata cross-checked against the log record; exit 0 only for a proof
PASS. Launch path (found before the queued launch): target niceness 10 instead of an increment (research-run already applies
nice 10); default worker counts 6 (the queue's cap; the sized launch uses 3); recert resume bound to the recorded N and code (a
different N would orphan claims), atomic claims, one process per worker ID, torn tails cut at restart, status-only resume
bookkeeping (memory); pre-replay completeness gate `certs_complete.py`; queued launch scripts and dry-run script.

## 8. Files
Certificates: `gcert.py`, `gsigma.py`, `ldl.py`, `tm.py`, `poly.py`, `chart.py`, `charts.py`, `refexact.py` (FEM); `bo1d.py`, `bo2.py`,
`p2up.py`, `p2fix.py`, `bocert.py` (1-D bounds and P2); `inc.py` (INC, Lemma 11.4); `nanguard.py` (fail-closed handling of non-finite
values in every verifier reduction); `bo2f.py` (float screening of slicing directions: heuristic only, never part of a certificate). Covers: `cover2.py` (old), `cover3.py`, `cover4.py`, `cover5.py`; families
`fam5.py`; runners `run_crit*.py`, `run_main_cbox.py`, `run_box.py`, `run_spike_all.py`, `run_fam.py`, `prod_*.sh`; the pool
`pool.py` (work items, claims, queue), `pool_launch.sh`, `pool_scaler.sh`, `pool_keepalive.sh`, `pool_status.py` (volume-based
ETA), `pool_inprogress.py`, `pool_strip.py` (removes a failed root's records before it is redone), `bench.py` (cost benchmark). Predicates
`regions.py`. Replay `replay2.py`, `replay3.py`, `replay_ok.py` (strict validation of shard outputs; reuse of an older replay is
for development only and never part of a proof run), `legacy_replays.sh`, `finalize.sh`, `accounting.py`, `manifests_legacy.py`, `merge_shards.py`, `merge_pool.py`,
`summarize.py` (runs table), `audit_farl_nan.py`; repair tools `binding.py`, `grid.py`, `schur.py`, `rigor_lint.py`,
`test_containment.py`, `check_charts.py`, `remanifest.py`, `recert.py`, `recert_all.sh` (re-certification), `certs_complete.py`
(pre-replay completeness gate), `replay_all.py` (aggregate driver; `replay_all.sh` is its wrapper), `bootstrap.py` (isolated
start of every proof process), `replay_common.py`, `fpenv.py`, `test_replay_infra.py`, `rigor_allowlist.json`, `make_pin.py` (pin
list), `launch_recert.sh` / `launch_replay.sh` / `smoke_launch.sh` (queued launch through research-run); review bundles `review/`
(`round2*`, `round3*`, `round4/`; `PINNED.sha256` is the current pin); tests of the repair `recert_tests/` (timing samples, pipeline and subdivision tests,
legacy replay samples; not part of the proof);
validation `validate_fem*.py`, `test_bo1d.py`, `test_inc.py`; lemma scripts `lemma_sector.py`, `lemma_longarm.py`, `hull_to_box.py`,
`check_conformity.py`, `check_cbox_cover.py`. Logs `certs/`, `certs2/`, `certs3/`; manifests `manifests/`; replay outputs
`replays/`; pilots `pilot_runs/` (not part of the proof). Proofs `LEMMAS.md`; pilot/cost notes `PILOT.md`.
