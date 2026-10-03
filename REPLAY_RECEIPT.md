# Replay receipt

Fresh-extraction replay of version `0.1.0-candidate` (date of publication 2026-10-03).

**What was replayed.** The component files of this archive (code, inputs, final logs, pins: everything under `local/`, `global2/`, `bandB/`
and `darts/`) were frozen as a replay input archive, extracted into a new directory, and every component was replayed there on
2–3 October 2026. Each run went through the producer's local job queue, with a clean environment and peak memory measured. The release
archive contains the same component files byte for byte; the release gate checks this against the replay input's manifest.

| Component | Command | Result | Wall time | Peak memory | Record |
|---|---|---|---|---|---|
| Theorem 2, local step | `local/replay_local_v2.sh` | 960 direction boxes, 0 failures, cube $1/20$; corrections, remainders and box records identical to the shipped ones | 10 min | 4.9 GB | `local/replay_v2/RECEIPT.json` `4210833462c2611b…` (fresh run: `fresh_replay/local_replay_v2_fresh/`) |
| Lemma $B'$, band | `bandB/replay.jl` (3 threads) | chart 1: 45 169 leaves, chart 2: 8 048 leaves, tiling exact; 44 268 certificates re-run; 0 errors | 3 min | 5.1 GB | `fresh_replay/band.out`, `fresh_replay/band_runs/` |
| Theorem 3, shell | `global2/replay.jl` (3 threads) | 8 217 leaves, 0 uncovered; 8 216 of 8 216 certificates re-verified; REPLAY OK | 30 min | 1.0 GB | `fresh_replay/shell.out` |
| Transport cover | `global2/replay_t.jl` (3 threads) | 50 183 leaves; 31 426 certificates re-verified and all exclusion and compound leaves re-checked; smallest certified margin 5.05e-06; REPLAY OK | 61 min | 2.5 GB | `global2/runs/transport_replay_fresh_RECEIPT.json` `1a14a85eec38f0bb…` |
| Theorem 5, darts | `darts/replay_all.sh` (SHARDS=2) | all 462 849 leaves of 14 logs; 8 gates, 8 merges, 14 replays; PASS (5.40 h) | 324 min | 1.6 GB | `darts/replays/receipts/20261002T221832_7c5eb4b8.jsonl` (final digest `bfe58c0b9ffe7759…`, file `9055fcbce245ca81…`) |

**Dart receipt verification.** The separately written state machine `verify/xreceipt.py` accepts the structure-only receipt
(`darts/replays/receipts/20261002T202211_d322f388.jsonl`) and the proof receipt, checking the proof receipt against the runtime anchors of the structure-only run
(driver inventory `5b2e45a77b8c5982…`, shard inventory `fcfbc91a85a1cab3…`) and the pin `58bd96bf9ee620e3…`.
Its mutation test (`verify/test_xreceipt.py`) rejects all 27 re-chained forgeries of the proof receipt.
Outputs: `fresh_replay/xreceipt_quick.out`, `fresh_replay/xreceipt_proof.out`, `fresh_replay/test_xreceipt_proof.out`.

**Smoke tests before the full runs.** These were partial replays of the shell, transport and band from the same extraction, and a
structure-only dart replay (`fresh_replay/smoke_*.out`, `fresh_replay/darts_quick.out`).

**Earlier replays.** Every component had previously been replayed exhaustively in the working tree:
- shell: `global2/runs/shell0225_replay_failclosed_full.log`;
- transport: `global2/runs/transport_replay_clean_RECEIPT.json`;
- band: `bandB/replay.log`;
- local: `local/replay_v2/RECEIPT.json`;
- darts: `darts/replays/receipts/20261002T130555_a1b1ffb8.jsonl`.

The dart receipt binds the bytes before the portability patch (`darts/review/PORTABILITY_PATCH_2026-10-02.md`).

**Scope.** These are producer-side replays on one machine (Apple M2, macOS 27.0). They are not an independent rerun by an unaffiliated
party, a reimplementation or a formal verification.

**Environment.** See `ENVIRONMENT.txt`. Python 3.13.5 with python-flint 0.9.0 (FLINT 3.6.0); Julia 1.12.7 with Arblib.jl 1.8.1.
