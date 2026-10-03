"""Mutation test of verify/xreceipt.py: a receipt that passes must pass, and every semantic mutation of it — re-chained with
correct canonical form, prev hashes and self hashes, so that only the state machine can catch it — must fail.
Usage: test_xreceipt.py PASSING_RECEIPT --pin-sha256 HEX --mode quick|proof [--darts DIR] [--recorded-darts DIR] [--home DIR]
       (the last three are passed to xreceipt.py)                exit 0 iff the base passes and every mutant fails"""
import sys, os, json, copy, hashlib, subprocess, argparse, tempfile

HERE = os.path.dirname(os.path.realpath(__file__))
PY = os.environ.get("PY") or sys.executable


def canon(o): return json.dumps(o, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
def sha(b): return hashlib.sha256(b).hexdigest()


def rechain(recs):
    out = []; prev = "0" * 64
    for i, r in enumerate(recs):
        r = {k: v for k, v in r.items() if k != "self_sha256"}
        r["seq"] = i; r["prev_sha256"] = prev
        r["self_sha256"] = sha(canon(r).encode()); prev = r["self_sha256"]
        out.append(r)
    return b"".join(canon(r).encode() + b"\n" for r in out)


import io, contextlib, functools
sys.path.insert(0, HERE)
import xreceipt as X                                     # the frozen verifier, executed in-process (same bytes as the CLI)
X.count_leaves = functools.lru_cache(maxsize=None)(X.count_leaves)   # inputs are immutable during the test: memoise the
X.sha_file = functools.lru_cache(maxsize=None)(X.sha_file)           # expensive reads (counts and hashes of ~600 MB of logs)


def run(raw, a, extra=()):
    with tempfile.NamedTemporaryFile("wb", suffix=".jsonl", delete=False) as fh:
        fh.write(raw); p = fh.name
    X.sha_file.cache_clear() if False else None
    argv = ["xreceipt.py", p, "--pin-sha256", a.pin_sha256, "--mode", a.mode] + list(extra) + list(a.passthrough)
    buf = io.StringIO(); old = sys.argv
    try:
        sys.argv = argv
        with contextlib.redirect_stdout(buf): rc = X.main()
    except SystemExit as e: rc = e.code
    finally:
        sys.argv = old; os.unlink(p)
    return rc, [l for l in buf.getvalue().splitlines() if l.startswith("PROBLEM")]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("receipt"); ap.add_argument("--pin-sha256", required=True)
    ap.add_argument("--mode", required=True, choices=("quick", "proof"))
    ap.add_argument("--darts"); ap.add_argument("--recorded-darts"); ap.add_argument("--home"); a = ap.parse_args()
    a.passthrough = [x for k in ("darts", "recorded_darts", "home") if getattr(a, k) for x in ("--" + k.replace("_", "-"), getattr(a, k))]
    base = [json.loads(l) for l in open(a.receipt, "rb").read().split(b"\n")[:-1]]
    ok = True
    rc, probs = run(rechain(base), a)
    print(f"{'ok  ' if rc == 0 else 'FAIL'} base receipt passes (exit {rc}) {probs[:3]}"); ok &= rc == 0
    R = lambda: copy.deepcopy(base)
    firstrep = next(i for i, r in enumerate(base) if r["kind"] == "replay")
    famrep = max(i for i, r in enumerate(base) if r["kind"] == "replay")
    m1 = next(i for i, r in enumerate(base) if r["kind"] == "merge")
    gate0 = next(i for i, r in enumerate(base) if r["kind"] == "gate")
    def mut(name, f, extra=()):
        nonlocal ok
        recs = R()
        try: f(recs)
        except IndexError:      # a mutation that needs >= 2 shards per log does not apply to a 1-shard receipt
            print(f"n/a  {name}  (needs at least 2 shards per log; this receipt has {base[0]['data'].get('shards')})"); return
        rc, probs = run(rechain(recs), a, extra)
        good = rc != 0
        ok &= good
        print(f"{'ok  ' if good else 'FAIL'} rejects: {name}  (exit {rc}; {probs[0][9:120] if probs else 'no problem reported'})")
    def setd(i, k, v): return lambda rs: rs[i]["data"].__setitem__(k, v)
    mut("one summary dropped from a replay", lambda rs: rs[famrep]["data"]["summaries"].pop())
    mut("all summaries of a log dropped", setd(firstrep, "summaries", []))
    mut("a replay record dropped", lambda rs: rs.pop(firstrep))
    mut("an invented extra gate", lambda rs: rs.insert(gate0, copy.deepcopy(rs[gate0])))
    mut("a gate command altered", lambda rs: rs[gate0]["data"].__setitem__("command", ["rigor_lint.py"]))
    mut("certs3 instead of certs4", setd(0, "certs", "certs3"))
    mut("duplicate merge (second = first)", lambda rs: rs.__setitem__(m1 + 1, copy.deepcopy(rs[m1])))
    mut("non-hex runtime inventory", lambda rs: rs[famrep]["data"]["summaries"][0]["environment"].__setitem__("inventory_sha256", "x"))
    mut("duplicate shard label", lambda rs: rs[famrep]["data"]["summaries"][1].__setitem__("shard", rs[famrep]["data"]["summaries"][0]["shard"]))
    mut("end status forged to the other mode's", lambda rs: rs[-1]["data"].__setitem__("status", "PASS" if a.mode == "quick" else "STRUCTURE-ONLY (NOT A PROOF)"))
    mut("child env DARTS_PIN changed", lambda rs: rs[0]["data"]["child_env"].__setitem__("DARTS_PIN", "/tmp/PIN"))
    mut("argv with --only", lambda rs: rs[0]["data"]["argv"].extend(["--only", "pp"]))
    mut("failures = False (bool, not int 0)", lambda rs: rs[famrep]["data"]["summaries"][0].__setitem__("failures", False))
    mut("leaves_checked and ok both +1", lambda rs: [rs[famrep]["data"]["summaries"][0].__setitem__(k, rs[famrep]["data"]["summaries"][0][k] + 1) for k in ("leaves_checked", "ok")])
    mut("replay log sha256 changed", setd(famrep, "log_sha256", "0" * 64))
    mut("two replay records swapped", lambda rs: rs.__setitem__(slice(firstrep, firstrep + 2), [rs[firstrep + 1], rs[firstrep]]))
    mut("start pin sha256 changed", setd(0, "pin_sha256", "f" * 64))
    mut("end record removed", lambda rs: rs.pop())
    mut("an abort record before end", lambda rs: rs.insert(len(rs) - 1, {"kind": "abort", "data": {"reason": "x"}, "time": rs[-1]["time"]}))
    mut("a shard at niceness 15", lambda rs: rs[famrep]["data"]["summaries"][0].__setitem__("niceness", 15))
    mut("a shard on another interpreter path", lambda rs: rs[famrep]["data"]["summaries"][0]["environment"].__setitem__("executable", "/opt/homebrew/bin/python3"))
    mut("an extra completeness input", lambda rs: rs[1]["data"]["inputs"].__setitem__("pp.w99.jsonl", "0" * 64))
    mut("shards = 2 in start (summaries for 3)", setd(0, "shards", 2) if base[0]["data"]["shards"] != 2 else setd(0, "shards", 3))
    mut("recert flag flipped in one summary", lambda rs: rs[famrep]["data"]["summaries"][0].__setitem__("recert", a.mode != "proof"))
    mut("validation_problems = 1", setd(famrep, "validation_problems", 1))
    inv = next(s["environment"]["inventory_sha256"] for r in base if r["kind"] == "replay" for s in r["data"]["summaries"])
    mut("shard runtime anchor mismatch", lambda rs: None, extra=("--shard-inv", "0" * 64))
    mut("driver runtime anchor mismatch", lambda rs: None, extra=("--driver-inv", "0" * 64))
    rc, probs = run(rechain(base), a, ("--shard-inv", inv, "--driver-inv", base[0]["data"]["runtime_identity"]["inventory_sha256"]))
    print(f"{'ok  ' if rc == 0 else 'FAIL'} base passes with its own anchors (exit {rc})"); ok &= rc == 0
    print("TEST_XRECEIPT", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
