"""Coordinator runbook, step 1 of the dart proof replay: prepare after the re-certification (fail-stop; prints the launch sizing).
Checks: the recert job has exited with status 0 and certs_complete reported complete; no recert / research-run process remains;
the re-certification code id still equals RECERT_CODE; the chart-cache directory equals the pinned inventory (no new caches);
the pin file has the reviewed digest and every pinned file its pinned hash.  Then freezes (chflags uchg) every certs4 worker
output and bookkeeping file and the directories that hold inputs (list appended to IMMUTABLE_DURING_RUNS.txt /
IMMUTABLE_DIRS_DURING_RUNS.txt), and computes SHARDS / --memory-gb from the final PP line count (rule: SHARDS = largest <= 4 with
SHARDS x 1.65 GB x L/621,307 + 0.4 <= 6.0 GB; memory rounded up to 0.5 GB).  Gates are run separately (my_gates.sh).
Usage: post_recert_prepare.py --pin-sha256 HEX        exit 0 iff ready"""
import os, sys, stat, json, glob, hashlib, subprocess, argparse, math

ROOT = "/Users/admin/claude-sandbox/quad-ppw-2026-09-30"
D = ROOT + "/darts"; C = D + "/certs4"
PY = "/Users/admin/.venvs/analysis-py313/bin/python"


def sha(p):
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def die(msg): print("NOT READY:", msg); sys.exit(1)


def freeze(paths, listfile):
    done = []
    for p in sorted(paths):
        st = os.lstat(p)
        if stat.S_ISLNK(st.st_mode): die(f"symlink {p}")
        if not st.st_flags & stat.UF_IMMUTABLE:
            os.chflags(p, st.st_flags | stat.UF_IMMUTABLE)
        done.append(p)
    with open(listfile, "a") as fh: fh.write("\n".join(done) + "\n")
    return len(done)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--pin-sha256", required=True); ap.add_argument("--dry-run", action="store_true"); a = ap.parse_args()
    if not a.dry_run:
        lo = open(C + "/launch.out").read()
        if "research-run exit 0" not in lo: die("launch.out has no 'research-run exit 0' (recert not finished or failed):\n" + lo[-800:])
        if '"complete": true' not in lo: die("certs_complete did not report complete")
        ps = subprocess.run(["/bin/ps", "-Ao", "command"], capture_output=True, text=True).stdout
        if "recert.py" in ps or "recert_all.sh" in ps: die("recert processes still running")
        rr = subprocess.run(["/Users/admin/.local/bin/research-run", "status"], capture_output=True, text=True).stdout
        if "darts-recert" in rr: die("research-run still lists darts-recert")
    code = subprocess.run(["/usr/bin/env", "-i", "PATH=/usr/bin:/bin", "HOME=/Users/admin", PY, "-I", "-S", "-B", D + "/bootstrap.py", D + "/recert.py", "--code-id"],
                          capture_output=True, text=True, cwd=D).stdout.strip()
    if code != open(C + "/RECERT_CODE").read().strip(): die(f"code id {code} != RECERT_CODE")
    pin_raw = open(D + "/review/PINNED.sha256", "rb").read()
    if hashlib.sha256(pin_raw).hexdigest() != a.pin_sha256: die("pin digest is not the reviewed one")
    pins = {}
    for l in pin_raw.decode().split("\n"):
        if l.strip():
            h, f = l.split("  ", 1); pins[f] = h
    bad = [f for f, h in pins.items() if sha(os.path.join(D, f)) != h]
    if bad: die(f"pinned files changed: {bad[:5]}")
    have = sorted(os.listdir(D + "/refcache")); want = sorted(f[len("refcache/"):] for f in pins if f.startswith("refcache/"))
    if have != want: die(f"chart-cache inventory differs from the pin: extra {sorted(set(have) - set(want))}, missing {sorted(set(want) - set(have))}")
    errs = [f for f in glob.glob(C + "/recert.w*.err") if os.path.getsize(f)]
    if errs: die(f"worker error logs not empty: {errs}")
    L = sum(sum(1 for _ in open(p, "rb")) for p in glob.glob(C + "/pp.w*.jsonl"))
    if a.dry_run:
        print(json.dumps(dict(dry_run=True, code_id=code, pins=len(pins), caches=len(have), pp_lines_so_far=L))); print("DRY-RUN OK"); return
    # freeze certs4 outputs and bookkeeping (writers have exited), then the input directories
    files = [p for p in glob.glob(C + "/*") if os.path.isfile(p)] + glob.glob(C + "/claims/*/*") + glob.glob(C + "/locks/*") + \
            [ROOT + "/verify/" + f for f in ("xreceipt.py", "test_xreceipt.py", "post_recert_prepare.py", "xdart.py", "xdart_run.py")]
    nf = freeze(files, ROOT + "/IMMUTABLE_DURING_RUNS.txt")
    dirs = [D, D + "/manifests", D + "/refcache", D + "/certs", D + "/certs2", D + "/certs3", D + "/review", D + "/review/round4",
            C + "/old", C + "/claims", C + "/locks", ROOT + "/verify"] + glob.glob(C + "/old/*") + glob.glob(C + "/claims/*")
    dirs = [p for p in dirs if os.path.isdir(p)]
    nd = freeze(dirs, ROOT + "/IMMUTABLE_DIRS_DURING_RUNS.txt")
    # the replay driver must still be able to write: certs4/ (merged logs) and replays/ stay writable
    for w in (C, D + "/replays"):
        if os.lstat(w).st_flags & stat.UF_IMMUTABLE: die(f"{w} must stay writable")
    shards = max([s for s in (1, 2, 3, 4) if s * 1.65 * L / 621307 + 0.4 <= 6.0] or [1])
    mem = math.ceil((shards * 1.65 * L / 621307 + 0.4) * 2) / 2
    print(json.dumps(dict(ready=True, code_id=code, pin=a.pin_sha256, frozen_files=nf, frozen_dirs=nd, pp_lines=L, shards=shards, memory_gb=mem)))
    print("READY")


if __name__ == "__main__":
    main()
