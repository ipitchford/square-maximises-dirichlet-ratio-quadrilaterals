# Agent notes to the coordinator (dart repair, 2026-10-02)

## Round 2 note (gates 1-2 pass; patched files; estimate) — sent after the first review
- Sol's round-1 claims verified before acting: the 81 slivers reproduced exactly (79 SP + 2 SF; measured widths 2^-59, 2^-56,
  2^-56 where the review said 2^-58, 2^-56, 2^-55); every witness path reproduced.
- Patched: nanguard (fdown/fup the single directed helper; amax/amin/upper_min/lower_max; exact_box, cuts, enclose, dec/decq),
  schur (ball-only Schur), grid (exact dyadic_root), binding, gcert, gsigma, p2up, tm (0-in-a branch), bo2, bo1d, charts
  (spike_tail, cut_tail), chart (cache_file, fresh builds), inc, fam5, regions, run_fam (exact_cell), replay2, replay3,
  replay_ok, replay_all.sh, accounting (rewritten), check_cbox_cover, hull_to_box, lemma_sector, manifests_legacy, merge_pool;
  manifests re-encoded; new tools rigor_lint, test_containment, check_charts, remanifest, recert; LEMMAS 5.1/5.2/5.3/9.3/9.4,
  new 9.6; README status.
- Results then: rigor_lint 37 files / 292 sites / 0 failures; test_containment 60,362 checks / 0 failures; check_charts 18 charts
  identical; accounting ok; recert samples 840 leaves / 0 failures; legacy replays pass; pipeline tests pass.
- Estimate: re-certification 10.4 CPU-h (1.5-2.1 h wall at 7 workers); fresh replay about 11.7 CPU-h (about 1.8 h); total about
  3.3 h (range 3-4.5 h).

## Round 2 addendum
- Self-review found 232 SP_col FEM leaves whose t-range exceeded t1 = float(0.02), where the collar chart is valid (Lemma 9.2
  covers t <= t1 only). The cover was not affected (the accounting relies only on the exact cell [0, t1]), but those leaves'
  certificates did not cover their whole boxes. SP.cart_box (collar) now refuses such boxes; re-certification subdivides them
  (6 timed: 0.3-2.8 s each).

## Round 3 note (this round: review round 2 BLOCK; SERIOUS 2-7 and MINOR items repaired; SERIOUS 1 not acted on)
- Claims verified first: wrong-family replay (crit_collar as SPIKE_X1) exited 0; SHARDS/QUICK unvalidated; replay_ok read only
  the summary; replays re-read the log for its hash; the binding lacked FP / binary / OS identity; no aggregate receipt; legacy
  roots without exact ends (4/4 crit_collar, 52/52 crit_spike_X1, 63/75 main_cbox, 256/256 shallow_crit, 378/882 kf_spike).
- Repairs: replay_all.py (aggregate driver; wrapper replay_all.sh validates SHARDS 1..64 canonical and QUICK unset or 1; only
  proof mode prints PASS; exactly SHARDS validated outputs per log; nothing reused; append-only hash-chained receipt binding the
  pin list, driver, gates, merges, shard outputs and summaries, runtime identity and the end-of-run unchanged state);
  replay_common.py (immutable snapshots, end-of-run identity checks, production modules verified unloaded); fpenv.py (IEEE
  round-to-nearest with gradual underflow asserted at start, per leaf and at the end; FPCR bits on arm64); binding.py
  (family, manifest bytes hash, interpreter/libpython, every file of python-flint/numpy/scipy/scikit-fem/mpmath, NumPy/BLAS
  config, FLINT version and context, OS, FP state, numerics env vars); replay_ok.py (every receipt validated against the
  recomputed leaf multiset); manifests carry their family; rigor_lint F5 NumPy taint + per-site allowlist
  (rigor_allowlist.json); test_containment exact Lemma 9.6 H1-H3 and fp:321 non-vacuity; k0 <= 0 rejected; docs qualified.
- Checks: rigor_lint 39 files, 451 annotated sites (169 in exempted scopes, each in rigor_allowlist.json), 0 failures;
  test_containment 60,594 checks, 0 failures; test_replay_infra 23 probes pass; check_charts OK; accounting ok; lemma tables,
  Corollary H', CBOX cover pass; driver machinery run with review/PINNED.sha256: 66 pinned files, none unpinned, receipt
  hash chain valid, end state unchanged.
- Bundle: review/PINNED.sha256 (the proof-run pin), review/round3_files.sha256 (everything reviewed), review/round3_incremental.diff
  (round 2 -> round 3), review/round3_full.diff (pre-review -> round 3; files without a pre-review copy listed in its header).
- Estimate unchanged: about 3.3 h wall at 7 workers (re-certification 1.5-2.1 h, fresh bound replay about 1.8 h); the new
  binding adds about 2.6 s per replay process and the per-leaf FP guard a few microseconds.
