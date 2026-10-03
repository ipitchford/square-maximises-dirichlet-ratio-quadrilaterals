# Dart proof replay: run records of 2 October 2026

These are the run records of the final dart pipeline: the launch outputs of the structure-only and proof replays, the outputs of the
separately written receipt verifier (`verify/xreceipt.py`) and its mutation test, the sampled certificate recomputations, the
sampled coverage check (`verify/xdart.py` on the re-certified logs) and the replay-infrastructure test.

Runbook: `verify/post_recert_prepare.py`, then gates, then the quick (structure-only) replay under the shared job queue, then
`verify/xreceipt.py --mode quick`, then `verify/test_xreceipt.py` (mutations), then sampled certificate recomputation
(`replay3 --shard 0/64` per family), then the proof replay (shared queue), then `verify/xreceipt.py --mode proof` with the quick
run's runtime anchors, then `XDART_CERTS=certs4 verify/xdart_run.py`.

The receipts of these two replays bind the pre-portability bytes of nine verifier files (`../../review/PORTABILITY_PATCH_2026-10-02.md`).
`portability/` contains the patch script, its complete diff and the infrastructure test after the patch. The replay that verifies
the shipped bytes is the fresh-extraction replay reported in the archive's `REPLAY_RECEIPT.md`.
