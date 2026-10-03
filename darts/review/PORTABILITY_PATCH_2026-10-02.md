# Portability patch of 2 October 2026

**What changed.** Before the frozen archive was built, nine pinned files of the dart verifier and four separately written checking or
test files contained the producer's absolute package path (`/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts`) or a fixed
interpreter path. Run from any other location, the verifier would have imported modules from that path. The runtime binding check is
fail-closed: a module loaded from outside the package aborts the replay. A copy of the archive could therefore not be replayed at all.

The patch replaces each such constant by the directory of the file itself (`os.path.dirname(os.path.realpath(__file__))`). The
launcher `replay_all.sh` now takes the package directory from its own physical location, with symlinks resolved, so a symlink to it
still selects the canonical driver, and takes the interpreter from `PY`. `verify/xreceipt.py` gained `--darts`, `--home` and
`--recorded-darts`, so that a receipt can be checked against a copy of the package at another path. The complete diff is
`PORTABILITY_PATCH_2026-10-02.diff` (78 changed lines in 13 files). No predicate, certificate, tolerance, log or control-flow
decision changed.

**Consequences.**
- The pin changed: `PINNED.pre_portability_2026-10-02.sha256` (SHA-256 `7e564f3f…`, 90 files) became `PINNED.sha256` (90 files,
  the same paths; exactly the nine patched pinned files differ).
- The proof replay of 2 October 2026 (receipt `replays/receipts/20261002T130555_a1b1ffb8.jsonl`, PASS) and its structure-only run
  (`20261002T121100_1f0f004a.jsonl`) bind the pre-patch bytes. They are kept as records. They can be verified only against a copy in
  which the reverse of the diff has been applied.
- The replay that verifies the shipped bytes is a new, fresh-extraction replay of the archive under the new pin. It is reported in
  the archive's `REPLAY_RECEIPT.md`.

**Checks after the patch.** `rigor_lint.py`: 0 failures. `test_replay_infra.py`: 82 of 82 probes pass, including the probe that a
symlinked launcher next to a counterfeit driver still runs the canonical driver (output in
`run_records/2026-10-02/portability/test_replay_infra_after_portability.txt`).
