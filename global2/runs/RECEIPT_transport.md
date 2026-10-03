# RECEIPT — transport cover (T): frozen log and clean replay

## 1. Frozen log (recorded 2026-10-01 15:37:25 BST before the clean replay started)

| item | value |
|---|---|
| file | `global2/runs/transport_cert.jsonl` (made read-only: chmod a-w) |
| SHA-256 | `cb187ebb05553ae4a45ad273e4e1304400fdc8246999e2234c1ea835eddd7bd7` |
| bytes | 31430611 |
| records | 100342 |
| pre-reboot prefix | `transport_cert.jsonl.bak_prereboot`, SHA-256 `0527947407eba898c89b92c4f22939321c6a5fb44cde5c2dda86595efa50fa99` (line prefix of the log; `check_log_t.py`) |
| identical to the log audited by Sol at 09:23:43 | yes |

Superseded: `superseded_unbound_replay_results.jsonl` (results of the earlier, resumed replay invocations of 09:00–10:35;
keyed by path only, not hash-bound, therefore not used by the clean replay) and its log `transport_replay.log`.

## 2. Clean replay (started 2026-10-01 15:38; results bound to log, verifier and record hashes)

Command (resumable; re-run the same command after any interruption — only hash-bound results are reused):
`nice -n 10 julia --heap-size-hint=1500M -t 3 replay_t.jl runs/transport_cert.jsonl --results runs/transport_replay_clean.jsonl`
Output: `runs/transport_replay_clean.log`; results `runs/transport_replay_clean.jsonl`; invocations
`runs/transport_replay_clean_invocations.jsonl`; machine receipt `runs/transport_replay_clean_RECEIPT.json` (written only
by an exhaustive, failure-free run).

### Result (2026-10-01 17:42:25)

| item | value |
|---|---|
| verdict | **REPLAY OK** — exhaustive: all 50182 non-SHELL leaves (31426 certified, 18756 exclusion/compound) re-verified under one verifier identity; 1 SHELL leaf by exact containment |
| log SHA-256 | `cb187ebb05553ae4a45ad273e4e1304400fdc8246999e2234c1ea835eddd7bd7` (100342 records, 50183 leaves: COMPOUND 3117, E2 13933, E3T 1568, E4 138, SHELL 1, certified 31426) |
| verifier SHA-256 | `4322dcc65eff1ccd0c1d93df584637c47590009c1b7330c56d22a3b12e9976e1` (julia 1.12.7; Arblib 1.8.1; PREC 128) |
| results SHA-256 | `c9fa7f8b6dca757d1fa0c2c0d99c77f3e51928156677d53e1ee346c621422c69` (`runs/transport_replay_clean.jsonl`, 50182 bound records) |
| receipt SHA-256 | `16ba10880dfd2b532ce055a2db1c80ccf0085da8b87634b9a8af836632312345` (`runs/transport_replay_clean_RECEIPT.json`) |
| counts | COMPOUND: 3117 checked, 0 failed, 0 TM-only, 280 thread-s; E2: 13933 checked, 0 failed, 0 TM-only, 57 thread-s; E3T: 1568 checked, 0 failed, 0 TM-only, 120 thread-s; E4: 138 checked, 0 failed, 0 TM-only, 0 thread-s; certified: 31426 checked, 0 failed, 0 TM-only, 21848 thread-s |
| smallest certified margin | 5.051419e-06 (largest κ* 2.499994948580528 at `T8.b0.b0.b1.b0.d1.b0.a1.c0.d0.b0.a1.c1.d0.b0.d1.a1.b0.c0`) |
| wall time | 7481 s over 1 invocation(s): 2026-10-01 15:37:42 → 2026-10-01 17:42:23 (7481 s, 50182 leaves, 3 threads) |

Verifier files (SHA-256):

| file | SHA-256 |
|---|---|
| `assemble.jl` | `4f19c660c82f0c848a77cb26a459903350865a8d191ce8dfd63d0e4430f3381b` |
| `bounds_tm.jl` | `5ceae71f70a03cc3a1ea30e7db34a3b6e99a8556fb8b8fb41420ee0df6e10990` |
| `config.toml` | `513d0ac93f748c2426b21a212623bced25caf29168198a8a3e64d2e8c4d20bca` |
| `cover.jl` | `939d131c4f9a460fecc3da58ef4a5e969ad4c5b986a9633db36390dbbd572ea1` |
| `cover_t.jl` | `53cb3882063a6e07fd2ae6d739690b77a8034d9383f451ac8213a5dd8aa9bdd8` |
| `criterion.jl` | `0eb9de883991e0a0b8a399b1f1e820cc0439d7bc4031eab0b170c78d958401ae` |
| `fsolve.jl` | `e2cb6bf847b80c69c3205380a3860293e094f908624b47c23e9ed13250e87576` |
| `geom.jl` | `e562e4a2dc2a9b39ac4248694ef95b43279ed6c0944f6b98958f83182d6164f7` |
| `logcheck.jl` | `09a82973f4ec816e816390eee9774b2625e5033cab7ef1e15a85ec916f53176a` |
| `replay_t.jl` | `3356f5a2f421cc8c99161bc31f6aaa4ddf40d8cc4b85b3ee72c76fe293d2d62e` |
| `rules.jl` | `4342c9abc3b5a20f57e84293316c56fa687b14bf52775004e92de89f80a874e4` |
| `strictjson.jl` | `7f5f15359241e66048d397d215f0f44be8840a1c0535197f8f360b5894abf817` |
| `tm.jl` | `d5f58430725b99acb9d860cfb297a04f78c97de301193278b76c93cfffc26562` |
| `transport.jl` | `e420613a22caec710c47e78f2fb10403888b24b5bb4c0993da96564511f11b30` |
| `transport1.jl` | `7f44f7c1d8f6d4b94eee8f160588db28bda9631729355a60fb016b3ef091da59` |
| `xcheck_cf.jl` | `a4f3206249a1ce0d7e0cf28cb66cfe920dd3323f3bc7d55432c6710930b2a4c5` |

Post-run check (2026-10-01 17:48): all 16 verifier files, the log and the results file on disk still have the SHA-256 values above: yes.
