"""Fill README_global2.md placeholders [[REPLAY0]], [[REPLAY]], [[MARGIN]] from the clean replay receipt, and append the
receipt summary to runs/RECEIPT_transport.md.  Run once, after the clean replay has written its receipt."""
import json, hashlib
R = json.load(open('runs/transport_replay_clean_RECEIPT.json'))
inv = [json.loads(l) for l in open('runs/transport_replay_clean_invocations.jsonl')]
inv = [q for q in inv if q['src_sha'] == R['verifier_sha256'] and q['log_sha'] == R['log_sha256']]
pk = R['per_kind']
g = lambda n: f"{n:,}".replace(",", " ")
nex = sum(v['checked'] for k, v in pk.items() if k != 'certified')
tm = sum(v['tm_only_leaves'] for v in pk.values())
wall = sum(q['wall'] for q in inv)
inv_txt = "; ".join(f"{q['start']} → {q['end']} ({q['wall']:.0f} s, {q['checked']} leaves, {q['threads']} threads)" for q in inv)
marg = R['smallest_certified_margin']
rep0 = (f"**clean, hash-bound replay of all {g(pk['certified']['checked'] + nex)} non-SHELL leaves: REPLAY OK** "
        f"(receipt `runs/RECEIPT_transport.md`, `runs/transport_replay_clean_RECEIPT.json`; log SHA-256 "
        f"`{R['log_sha256'][:16]}…`). Smallest certified margin 5/2 − κ* = {marg:.3e}.")
rep = (f"the **clean** replay (no cached results; started 2026-10-01 15:37 from an empty results file, 3 threads at nice 10) "
       f"re-verified **every one of the {g(pk['certified']['checked'] + nex)} non-SHELL leaves** of the frozen log "
       f"(SHA-256 `{R['log_sha256']}`): {g(pk['certified']['checked'])} certificates recomputed from the log alone "
       f"(Float64 choice at the logged degree/diagonal/first-order flag, Arb reference integrals, Arb Taylor models, the "
       f"logged lower-bound candidate rebuilt in Arb, Arb criterion at exactly κ = 5/2; for the 3 262 ladder-fallback "
       f"certificates additionally a certified κ* by exact-rational bisection), and {g(nex)} exclusion/compound leaves "
       f"re-checked by the centred-form checker (TM-only fallbacks: {tm}); the SHELL leaf by exact containment. "
       f"Failures: **0**. Coverage rebuilt from the 24 roots: 0 uncovered, 0 mismatches, 0 unreached records. Every result "
       f"line is bound to the SHA-256 of its log line, of the log and of the verifier (`{R['verifier_sha256'][:16]}…`: "
       f"{len(R['verifier_files'])} source/config files + Julia/Arblib versions + precision). Wall time {wall:.0f} s "
       f"over {len(inv)} invocation(s) ({inv_txt}). Results file `runs/transport_replay_clean.jsonl`, SHA-256 "
       f"`{R['results_sha256']}`; output `runs/transport_replay_clean.log`; receipt `runs/transport_replay_clean_RECEIPT.json` "
       f"(summary in `runs/RECEIPT_transport.md`).\n\n"
       f"*Superseded:* the earlier replay invocations of 09:00–10:35 (`runs/transport_replay.log`, results "
       f"`runs/superseded_unbound_replay_results.jsonl`) also ended with REPLAY OK and 0 failures, but they resumed from "
       f"results keyed only by path while the verifier was being edited; they are not used as evidence.")
mrow = (f"{marg:.3e} (exact-rational bisection in the clean replay; largest κ* = {R['largest_bisected_kappa_star']:.15f} "
        f"at `{R['at']}`); every certificate holds at κ = 5/2 with strict inequalities; the ladder-logged certificates "
        f"have margin ≥ {R['ladder_min_margin']:.3e}")
s = open('README_global2.md').read()
for k, v in (("[[REPLAY0]]", rep0), ("[[REPLAY]]", rep), ("[[MARGIN]]", mrow)):
    assert k in s, k
    s = s.replace(k, v, 1)
s = s.replace("**Replay** (`nice -n 10 julia -t 4 replay_t.jl runs/transport_cert.jsonl`, resumable; results\n`runs/transport_cert_replay_results.jsonl`, log `runs/transport_replay.log`): ",
              "**Replay** (`nice -n 10 julia --heap-size-hint=1500M -t 3 replay_t.jl runs/transport_cert.jsonl --results runs/transport_replay_clean.jsonl`): ", 1)
open('README_global2.md', 'w').write(s)
res_sha = hashlib.sha256(open('runs/transport_replay_clean.jsonl', 'rb').read()).hexdigest()
assert res_sha == R['results_sha256'], "results file changed after the receipt"
with open('runs/RECEIPT_transport.md', 'a') as f:
    f.write(f"""
### Result ({R['finished']})

| item | value |
|---|---|
| verdict | **{R['verdict']}** — {R['scope']} |
| log SHA-256 | `{R['log_sha256']}` ({R['log_records']} records, {R['leaves']} leaves: {', '.join(f'{k} {v}' for k, v in R['leaves_by_status'].items())}) |
| verifier SHA-256 | `{R['verifier_sha256']}` ({R['verifier_extra']}) |
| results SHA-256 | `{R['results_sha256']}` (`runs/transport_replay_clean.jsonl`, {R['results_bound_records']} bound records) |
| receipt SHA-256 | `{hashlib.sha256(open('runs/transport_replay_clean_RECEIPT.json','rb').read()).hexdigest()}` (`runs/transport_replay_clean_RECEIPT.json`) |
| counts | {'; '.join(f"{k}: {v['checked']} checked, {v['failed']} failed, {v['tm_only_leaves']} TM-only, {v['thread_s']:.0f} thread-s" for k, v in pk.items())} |
| smallest certified margin | {marg:.6e} (largest κ* {R['largest_bisected_kappa_star']:.15f} at `{R['at']}`) |
| wall time | {wall:.0f} s over {len(inv)} invocation(s): {inv_txt} |

Verifier files (SHA-256):

| file | SHA-256 |
|---|---|
""" + "".join(f"| `{k}` | `{v}` |\n" for k, v in sorted(R['verifier_files'].items())))
print("README and RECEIPT updated")
