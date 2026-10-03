# Round 4: notes for the reviewer (2026-10-02, Anonymous)

Scope: the round-3 review items 1-7 (harness / provenance) and the launch-path defects found before the queued launch.
Threat model (README §5): accidental corruption, stale or changed files, concurrent agents, environment drift; deliberate
tampering by the machine owner is out of scope. Pin: `review/PINNED.sha256` (90 files, sha256 0ba05176...681b50).
Files changed since the round-3 pin, tagged: `changed_since_round3.txt`; diff against the verified round-3 versions:
`incremental.diff`. No KERNEL semantic change: gcert.py, ldl.py and hull_to_box.py changed by RIGOR comments only
(AST-identical to round 3, checked by reconstruction); chart.py changed only in the cache-loading path.

## Round-3 items
1. Runtime inventory (`binding.py`, `bootstrap.py`, `replay_common.py`): computed fresh at every call, never memoised. It covers
   every loaded non-OS Mach-O image (dyld, which includes the Python framework binary), the interpreter, the stdlib, the files of
   the numerical packages, and every *.pth and top-level module of the site directories. Module origins are checked against
   the inventory, the pinned sources, or the installers' RECORD hashes for late imports. Cached bytecode is never read (-B and
   an empty pycache prefix). sys.flags, _xoptions and the pycache prefix are bound. `python -O` is refused.
2. Chart caches (`chart.py`, `check_charts.py`, `make_pin.py`): the refcache/* inventory is pinned. It is verified before any
   unpickling, before each replay and at the end. A pinned run never writes a cache. The rebuild-and-compare gate is kept.
3. Pin and inputs (`replay_all.py`): the pin is read once and those bytes are parsed and hashed; the children get a private
   read-only copy. Each log and manifest snapshot must equal its pin entry or the recorded merge output.
4. Shard outputs: created with O_EXCL and snapshotted once after the child exits. Exactly those bytes are validated, summarised
   and hashed. A missing summary fails the run (`summaries_of`).
5. Receipt: per record, a self hash, strict sequence and schema checks, a write loop, and fsync of the file and directory. The
   receipt is re-read before success and the final digest is printed.
6. Bootstrap: the package path and interpreter are absolute constants. Niceness is set in-process. Children run under env -i
   with -I -S -B.
7. Lint aliases and string-versus-comment detection; Lemma 9.6 also checked by analytic extrema; FIZ / AH / FZ16 probes (on
   this CPU FZ and FZ16 are detected, AH and FIZ do not stick, i.e. FEAT_AFP is absent); receipt metadata is cross-checked
   against the log; exit 0 only for a proof PASS.

## Launch path (research-run execs jobs as `/usr/bin/nice -n 10 CMD`)
* Double nice: `replay_all.py` does `os.nice(max(0, 10 - os.nice(0)))`. `recert_all.sh` takes its increment from
  `/bin/ps -o nice= -p $$`. Probes: under nice -n 10, and from niceness 0, the driver and every child (and every recert worker)
  run at 10.
* Worker counts: the defaults are 6 (replay_all.sh, replay_all.py, recert_all.sh). The sized launch uses 3.
* Resume: N, the code identity and the families are recorded at the first start. A resume with any other value is refused (exit
  2), including the default N. Claims are created atomically (a file already holding the ID is linked into place). A claim that
  names no worker 1..N stops the worker. A second process with the same worker ID is refused (flock). A torn tail is cut off at
  the next start. Claims are not adopted: a root's partial records would sit in another worker's output and duplicate on merge.
  With N fixed, orphans cannot arise.
* Completeness gate (`certs_complete.py`, a few seconds): runs at the end of recert_all.sh, and in the driver before any gate,
  merge or replay (abort). It lists every FAIL leaf. The merge inputs must equal the bytes it checked.
* Memory (peak RSS from /usr/bin/time -l):
  * replay3 shard of the PP log: 1.50 GB.
  * driver validation of PP: 0.74 GB, after the shards exit.
  * recert worker on heavy PP roots: 0.43 GB.
  * PP-share resume bookkeeping: 59-73 MB. This is status-only now; with whole records it was ~0.3 GB.
  Sizing: recert W = 3 at 2.5 GB; replay SHARDS = 3 at 5.5 GB, rescaled from the certs4 PP line count by `launch_replay.sh`.
* Dry runs: NOT performed under research-run. Executing `research-run` (even `--help` / `status`) was denied to this agent by the
  permission system; no workaround was attempted, and whether the dry runs are run is the user's decision. `smoke_launch.sh
  recert|replay` holds the exact commands. Local equivalents under env -i and nice -n 10 (cwd /), see `dry_runs_local.txt`:
  * recert FF, N=1: exit 0, 15 s, niceness 10, complete.
  * replay test mode with all 8 gates, crit_collar plus FF, SHARDS=1: exit 3 (MACHINERY TEST OK, NOT A PROOF), 147 s. Driver
    and children at niceness 10; final digest a45bce19...
* Gates on the pinned code: rigor_lint (45 files, 524 annotated sites, 0 failures), test_containment (62,767 checks, 0
  failures), check_charts (all caches rebuilt IDENTICAL), accounting ok, and the lemma gates. test_replay_infra: 82 probes, 0
  failures (`test_replay_infra_output.txt`).
* Recert samples (300 per family and method, seed 11): no failures anywhere except SP_col FEM, 27/300 = 9%. These leaves fail
  with 'FEM: cart box': their t-range exceeds the collar t1 = 0.02. recert subdivides them; an earlier test certified them by
  BO2 after one or two splits.
* INC: recert recomputes E. The 50 old leaves whose stored E exceeded the fresh bound all re-certify.
