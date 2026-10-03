"""Independent post-run verification of a dart replay receipt (darts/replays/receipts/<ID>.jsonl) — exact state machine (v2).

Written separately from the replay driver (it does not import it). v1 checked properties of whatever records were present; the
round-5 adversarial review reproduced canonical, correctly chained receipts with missing shards, missing logs, invented gates,
duplicate merges or the wrong certificate directory that v1 accepted. v2 accepts exactly one record sequence:

    start, complete, gate x 8 (driver order), merge x 8 (family order), replay x 14 (6 legacy logs, then 8 families), end

and checks every record's fields exactly (types strict: an int is never a bool), the cross-links between records (pin entries,
merge outputs, completeness inputs, binding, child environment, runtime identity), and the files as they are NOW (pin, pinned
files, completeness inputs = the certs4 directory listing, merged logs, number of leaf records per log = sum of leaves checked).
Modes: --mode proof (status PASS, every summary recert = True) or --mode quick (structure-only run: status
'STRUCTURE-ONLY (NOT A PROOF)', recert = False).  Runtime anchors (--driver-inv / --shard-inv, from an EARLIER run such as the quick
run) make the runtime identity an externally anchored precondition.
Usage: xreceipt.py RECEIPT --pin-sha256 HEX [--mode proof|quick] [--driver-inv HEX] [--shard-inv HEX] [--darts DIR] [--home DIR]
  --darts: the package directory the receipt was produced in (default: ../darts beside this file, physical path);
  --home: the HOME of the run (default: $HOME);
  --recorded-darts: the package directory recorded in the receipt, if the files are now read from another copy (--darts):
  recorded absolute paths are compared with this root, files are read from --darts.          exit 0 iff XRECEIPT PASS
"""
import sys, os, json, hashlib, re, argparse

DARTS = os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "darts"))
GATES = [["lemma_sector.py"], ["lemma_longarm.py"], ["check_charts.py", "8"], ["rigor_lint.py"], ["test_containment.py"],
         ["hull_to_box.py", "56", "64", "0.0339", "0.12", "0.5", "0.875"], ["check_cbox_cover.py"], ["accounting.py"]]
LEGACY = [("certs/crit_collar.jsonl", "COLLAR", "manifests/crit_collar.json", "crit_collar"),
          ("certs/crit_spike_X1.jsonl", "SPIKE_X1", "manifests/crit_spike_X1.json", "crit_spike_X1"),
          ("certs/main_cbox.jsonl", "MAIN_CBOX", "manifests/main_cbox.json", "main_cbox"),
          ("certs/shallow_crit.jsonl", "BOX", "manifests/shallow_crit.json", "shallow_crit"),
          ("certs2/kf_collar.jsonl", "KF_COLLAR", "manifests/kf_collar.json", "kf_collar"),
          ("certs2/kf_spike.jsonl", "KF_SPIKE", "manifests/kf_spike.json", "kf_spike")]
FAMILIES = ["PP", "SP_col", "SP_spk", "PF", "FP", "FF", "SF_col", "SF_spk"]
CERTS = "certs4"
SAME_AS_DRIVER = ("executable", "python", "implementation", "packages", "flint", "numpy_config", "os", "arb_prec", "TM_ORDER",
                  "tm_order_effective", "flags")
HEX = re.compile(r"[0-9a-f]{64}")


def sha_bytes(b): return hashlib.sha256(b).hexdigest()
def sha_file(p):
    with open(p, "rb") as fh: return sha_bytes(fh.read())
def canon(o): return json.dumps(o, sort_keys=True, separators=(",", ":"), ensure_ascii=True)   # the producer's canonical form
def is_int(x): return isinstance(x, int) and not isinstance(x, bool)
def is_hex(x): return isinstance(x, str) and HEX.fullmatch(x) is not None


class V:
    def __init__(self): self.problems = []
    def req(self, cond, msg):
        if not cond: self.problems.append(msg)
        return bool(cond)


def parse_chain(raw, v):
    recs = []; prev = "0" * 64
    if not v.req(raw.endswith(b"\n"), "receipt does not end with a newline (torn final record)"): return recs
    for i, line in enumerate(raw.split(b"\n")[:-1]):
        try: r = json.loads(line)
        except Exception: v.req(False, f"record {i}: not JSON"); break
        if not v.req(isinstance(r, dict) and set(r) == {"seq", "prev_sha256", "time", "kind", "data", "self_sha256"}, f"record {i}: keys"): break
        if not v.req(canon(r).encode() == line, f"record {i}: not in canonical form"): break
        if not v.req(is_int(r["seq"]) and r["seq"] == i, f"record {i}: seq {r['seq']!r}"): break
        if not v.req(r["prev_sha256"] == prev, f"record {i}: chain broken"): break
        body = {k: x for k, x in r.items() if k != "self_sha256"}
        if not v.req(sha_bytes(canon(body).encode()) == r["self_sha256"], f"record {i}: self_sha256 mismatch"): break
        if not v.req(isinstance(r["data"], dict), f"record {i}: data is not an object"): break
        prev = r["self_sha256"]; recs.append(r)
    return recs


def count_leaves(path):
    n = 0
    with open(path, "rb") as fh:
        for line in fh:
            if b'"leaf"' in line:
                d = json.loads(line)
                if d.get("type") == "leaf": n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("receipt"); ap.add_argument("--pin-sha256", required=True); ap.add_argument("--mode", default="proof", choices=("proof", "quick"))
    ap.add_argument("--driver-inv"); ap.add_argument("--shard-inv")
    ap.add_argument("--darts", default=DARTS); ap.add_argument("--home", default=os.environ.get("HOME", ""))
    ap.add_argument("--recorded-darts", default=None)
    a = ap.parse_args()
    globals()["DARTS"] = os.path.realpath(a.darts)
    REC = a.recorded_darts or DARTS          # recorded root (string comparisons); DARTS: where files are read
    v = V()
    raw = open(a.receipt, "rb").read()
    recs = parse_chain(raw, v)
    kinds = [r["kind"] for r in recs]
    expect = ["start", "complete"] + ["gate"] * len(GATES) + ["merge"] * len(FAMILIES) + ["replay"] * (len(LEGACY) + len(FAMILIES)) + ["end"]
    if not v.req(kinds == expect, f"record sequence {kinds} != the exact expected sequence ({len(expect)} records)"):
        return report(v, None)
    D = [r["data"] for r in recs]
    st, cp, gates, merges, reps, end = D[0], D[1], D[2:10], D[10:18], D[18:32], D[32]
    # ---------------------------------------------------------------- pin (on disk, reviewed digest) and pinned files
    pin_path = os.path.join(REC, "review/PINNED.sha256")
    pin_raw = open(os.path.join(DARTS, "review/PINNED.sha256"), "rb").read()
    v.req(sha_bytes(pin_raw) == a.pin_sha256, "review/PINNED.sha256 on disk is not the reviewed pin")
    pins = {}
    for line in pin_raw.decode().split("\n"):
        if line.strip():
            h, f = line.split("  ", 1); pins[f] = h
    bad = [f for f, h in pins.items() if sha_file(os.path.join(DARTS, f)) != h]
    v.req(not bad, f"{len(bad)} pinned files differ from the pin now: {bad[:5]}")
    # ---------------------------------------------------------------- start
    N = st.get("shards")
    v.req(is_int(N) and 1 <= N <= 64, f"start.shards {N!r}")
    N = N if is_int(N) and 1 <= N <= 64 else 1
    v.req(st.get("mode") == a.mode, f"start.mode {st.get('mode')!r} != {a.mode}")
    v.req(st.get("certs") == CERTS, f"start.certs {st.get('certs')!r}")
    v.req(st.get("pin_sha256") == a.pin_sha256, "start.pin_sha256 is not the reviewed pin")
    v.req(st.get("pin_file") == pin_path, f"start.pin_file {st.get('pin_file')!r}")
    v.req(is_int(st.get("pin_count")) and st.get("pin_count") == len(pins), f"start.pin_count {st.get('pin_count')!r} != {len(pins)}")
    v.req(st.get("argv") == [os.path.join(REC, "replay_all.py"), "--pin", pin_path, "--shards", str(N), "--mode", a.mode, "--certs", CERTS],
          f"start.argv is not the plain {a.mode} invocation: {st.get('argv')}")
    v.req(st.get("driver_sha256") == pins.get("replay_all.py") and is_hex(pins.get("replay_all.py")), "start.driver_sha256 is not the pinned replay_all.py")
    v.req(st.get("wrapper_sha256") == pins.get("replay_all.sh") and is_hex(pins.get("replay_all.sh")), "start.wrapper_sha256 is not the pinned replay_all.sh")
    v.req(st.get("niceness") == 10 and is_int(st.get("niceness")), f"start.niceness {st.get('niceness')!r}")
    v.req(st.get("failures_so_far") == [], f"start.failures_so_far {st.get('failures_so_far')!r}")
    run_id = st.get("run_id")
    v.req(isinstance(run_id, str) and re.fullmatch(r"[0-9]{8}T[0-9]{6}_[0-9a-f]{8}", run_id) is not None, f"start.run_id {run_id!r}")
    run_dir = os.path.join(DARTS, "replays", f"run_{run_id}"); run_dir_rec = os.path.join(REC, "replays", f"run_{run_id}")
    v.req(st.get("pin_copy") == os.path.relpath(os.path.join(run_dir_rec, "PIN.sha256"), REC), f"start.pin_copy {st.get('pin_copy')!r}")
    want_env = {"PATH": "/usr/bin:/bin", "HOME": a.home, "TMPDIR": os.path.join(run_dir_rec, "tmp"), "OMP_NUM_THREADS": "1",
                "OPENBLAS_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1", "DARTS_PIN": os.path.join(run_dir_rec, "PIN.sha256")}
    v.req(st.get("child_env") == want_env, f"start.child_env {st.get('child_env')!r}")
    ident = st.get("runtime_identity")
    v.req(isinstance(ident, dict) and is_hex(ident.get("inventory_sha256")) and is_int(ident.get("inventory_count")), "start.runtime_identity malformed")
    ident = ident if isinstance(ident, dict) else {}
    if a.driver_inv: v.req(ident.get("inventory_sha256") == a.driver_inv, f"driver runtime inventory {ident.get('inventory_sha256')} != anchor {a.driver_inv}")
    v.req(not (ident.get("fp") or {}).get("behaviour_failures"), "driver FP behaviour failures")
    # ---------------------------------------------------------------- completeness
    v.req(cp.get("certs") == CERTS and cp.get("families") == FAMILIES, f"complete.certs/families {cp.get('certs')!r} {cp.get('families')!r}")
    recert_n = open(os.path.join(DARTS, CERTS, "RECERT_N")).read().strip()
    v.req(is_int(cp.get("N")) and str(cp.get("N")) == recert_n, f"complete.N {cp.get('N')!r} != RECERT_N {recert_n}")
    v.req(cp.get("problems") == 0 and is_int(cp.get("problems")) and cp.get("fail_leaves") == 0 and is_int(cp.get("fail_leaves")), "complete: problems / FAIL leaves")
    man_roots = {F: len(json.load(open(os.path.join(DARTS, "manifests", F.lower() + ".json")))["roots"]) for F in FAMILIES}
    v.req(cp.get("roots") == man_roots, f"complete.roots {cp.get('roots')!r} != manifest roots {man_roots}")
    inputs = cp.get("inputs") if isinstance(cp.get("inputs"), dict) else {}
    on_disk = sorted(f for f in os.listdir(os.path.join(DARTS, CERTS)) if f.endswith(".jsonl") and not f.endswith(".merged.jsonl")
                     and any(f.startswith(F.lower() + ".") for F in FAMILIES))
    v.req(sorted(inputs) == on_disk, f"complete.inputs {len(inputs)} files != the certs4 worker logs on disk now {len(on_disk)}")
    for f, h in inputs.items():
        p = os.path.join(DARTS, CERTS, f)
        v.req(is_hex(h) and os.path.isfile(p) and sha_file(p) == h, f"input {f} changed or missing")
    # ---------------------------------------------------------------- gates
    for g, d in zip(GATES, gates):
        v.req(d.get("command") == g and d.get("exit") == 0 and is_int(d.get("exit")) and is_hex(d.get("output_sha256")), f"gate {g}: {d.get('command')!r} exit {d.get('exit')!r}")
        out = os.path.join(DARTS, d.get("output") or "/nonexistent")
        v.req(os.path.dirname(out) == run_dir and os.path.isfile(out) and sha_file(out) == d.get("output_sha256"), f"gate {g}: output file missing or changed")
    # ---------------------------------------------------------------- merges
    merged_sha = {}
    for F, d in zip(FAMILIES, merges):
        L = F.lower(); mp = f"{CERTS}/{L}.merged.jsonl"
        want_in = {f: h for f, h in inputs.items() if f.startswith(L + ".")}
        v.req(d.get("family") == F and d.get("exit") == 0 and is_int(d.get("exit")) and d.get("merged") == mp, f"merge {F}: family/exit/path {d.get('family')!r} {d.get('exit')!r} {d.get('merged')!r}")
        v.req(d.get("inputs") == want_in and len(want_in) > 0, f"merge {F}: inputs differ from the completeness inputs")
        v.req(is_hex(d.get("merged_sha256")) and os.path.isfile(os.path.join(DARTS, mp)) and sha_file(os.path.join(DARTS, mp)) == d.get("merged_sha256"), f"merge {F}: merged log missing or changed")
        merged_sha[F] = d.get("merged_sha256")
    # ---------------------------------------------------------------- replays
    plan = [(log, fam, man, "replay2.py", tag, pins.get(log)) for log, fam, man, tag in LEGACY] + \
           [(f"{CERTS}/{F.lower()}.merged.jsonl", F, f"manifests/{F.lower()}.json", "replay3.py", F.lower(), merged_sha.get(F)) for F in FAMILIES]
    shard_invs = set(); total_checked = 0; nsum = 0
    for (log, fam, man, script, tag, want_log), d in zip(plan, reps):
        pre = f"replay {tag}"
        v.req((d.get("tag"), d.get("script"), d.get("family"), d.get("log"), d.get("manifest")) == (tag, script, fam, log, man), f"{pre}: identity {d.get('tag')!r} {d.get('script')!r} {d.get('family')!r} {d.get('log')!r}")
        v.req(is_hex(want_log) and d.get("log_sha256") == want_log, f"{pre}: log sha256 is not the pinned / merged one")
        v.req(d.get("manifest_sha256") == pins.get(man) and is_hex(pins.get(man)), f"{pre}: manifest sha256 is not the pinned one")
        v.req(d.get("shards") == N and is_int(d.get("shards")), f"{pre}: shards {d.get('shards')!r}")
        v.req(d.get("exit_codes") == {str(i): 0 for i in range(N)}, f"{pre}: exit codes {d.get('exit_codes')!r}")
        outs = d.get("outputs") if isinstance(d.get("outputs"), dict) else {}
        v.req(sorted(outs) == sorted(f"{tag}.{i}.jsonl" for i in range(N)) and all(is_hex(h) for h in outs.values()), f"{pre}: outputs {sorted(outs)}")
        for n_, h in outs.items():
            p = os.path.join(run_dir, n_)
            v.req(os.path.isfile(p) and sha_file(p) == h, f"{pre}: shard output {n_} missing or changed")
        v.req(d.get("validation_problems") == 0 and is_int(d.get("validation_problems")), f"{pre}: validation problems {d.get('validation_problems')!r}")
        sums = d.get("summaries")
        if not v.req(isinstance(sums, list) and len(sums) == N, f"{pre}: {len(sums) if isinstance(sums, list) else sums!r} summaries, need {N}"): continue
        vh = d.get("expected_binding")
        v.req(is_hex(vh), f"{pre}: expected_binding {vh!r}")
        checked = 0
        for i, sm in enumerate(sums):
            sp = f"{pre} shard {i}/{N}"
            if not v.req(isinstance(sm, dict), f"{sp}: summary is not an object"): continue
            v.req(sm.get("shard") == f"{i}/{N}" and sm.get("log") == log and sm.get("family") == fam, f"{sp}: shard/log/family {sm.get('shard')!r} {sm.get('log')!r} {sm.get('family')!r}")
            v.req(sm.get("partial") is False, f"{sp}: partial {sm.get('partial')!r}")
            v.req(sm.get("recert") is (a.mode == "proof"), f"{sp}: recert {sm.get('recert')!r}")
            v.req(sm.get("failures") == 0 and is_int(sm.get("failures")), f"{sp}: failures {sm.get('failures')!r}")
            v.req(is_int(sm.get("leaves_checked")) and is_int(sm.get("ok")) and sm.get("ok") == sm.get("leaves_checked"), f"{sp}: ok {sm.get('ok')!r} != leaves_checked {sm.get('leaves_checked')!r}")
            v.req(sm.get("niceness") == 10 and is_int(sm.get("niceness")), f"{sp}: niceness {sm.get('niceness')!r}")
            v.req(sm.get("binding_sha256") == vh and sm.get("log_sha256") == want_log and sm.get("manifest_sha256") == pins.get(man), f"{sp}: binding / log / manifest hash")
            bs = sm.get("binding_sources") if isinstance(sm.get("binding_sources"), dict) else {}
            v.req(bool(bs) and all(pins.get(s) == h and is_hex(h) for s, h in bs.items()), f"{sp}: binding sources not all pinned")
            v.req({"refcache/" + f: h for f, h in (sm.get("binding_refcache") or {}).items()} == {f: h for f, h in pins.items() if f.startswith("refcache/")}, f"{sp}: chart caches != pin")
            env = sm.get("environment") if isinstance(sm.get("environment"), dict) else {}
            v.req(env.get("env") == {k: x for k, x in want_env.items() if k.startswith(("OMP_", "OPENBLAS_", "VECLIB_", "DARTS_"))}, f"{sp}: child environment")
            v.req(env.get("pycache_prefix") == os.path.join(run_dir_rec, "pycache"), f"{sp}: pycache prefix {env.get('pycache_prefix')!r}")
            v.req(not (env.get("fp") or {}).get("behaviour_failures"), f"{sp}: FP behaviour failures")
            v.req(is_hex(env.get("inventory_sha256")) and is_int(env.get("inventory_count")), f"{sp}: runtime inventory malformed")
            shard_invs.add(env.get("inventory_sha256"))
            for k in SAME_AS_DRIVER:
                if not v.req(k in env and k in ident and json.dumps(env.get(k), sort_keys=True) == json.dumps(ident.get(k), sort_keys=True), f"{sp}: runtime field {k} differs from the driver's"): break
            checked += sm.get("leaves_checked") if is_int(sm.get("leaves_checked")) else 0
        n_leaves = count_leaves(os.path.join(DARTS, log)) if os.path.isfile(os.path.join(DARTS, log)) else -1
        v.req(checked == n_leaves, f"{pre}: shards checked {checked} leaves, the log has {n_leaves}")
        total_checked += checked; nsum += len(sums)
    v.req(len(shard_invs) == 1, f"{len(shard_invs)} distinct runtime inventories across all shards of all logs")
    if a.shard_inv: v.req(shard_invs == {a.shard_inv}, f"shard runtime inventory {sorted(map(str, shard_invs))} != anchor {a.shard_inv}")
    # ---------------------------------------------------------------- end
    want_status = "PASS" if a.mode == "proof" else "STRUCTURE-ONLY (NOT A PROOF)"
    v.req(end.get("status") == want_status, f"end.status {end.get('status')!r} != {want_status!r}")
    v.req(end.get("failures") == [] and end.get("end_problems") == [] and end.get("end_state_unchanged") is True, f"end: failures {end.get('failures')!r} problems {end.get('end_problems')!r}")
    stats = dict(records=len(recs), shards_per_log=N, summaries=nsum, leaves_checked=total_checked, final_digest=recs[-1]["self_sha256"],
                 receipt_file_sha256=sha_bytes(raw), driver_inventory=ident.get("inventory_sha256"),
                 shard_inventory=sorted(map(str, shard_invs))[0] if shard_invs else None, verifier_sha256=sha_file(os.path.realpath(__file__)))
    return report(v, stats)


def report(v, stats):
    print("XRECEIPT", json.dumps(stats))
    for p in v.problems[:200]: print("PROBLEM:", p)
    if len(v.problems) > 200: print(f"PROBLEM: ... {len(v.problems) - 200} more")
    print("XRECEIPT", "PASS" if not v.problems else f"FAIL ({len(v.problems)} problems)")
    return 0 if not v.problems else 1


if __name__ == "__main__":
    sys.exit(main())
