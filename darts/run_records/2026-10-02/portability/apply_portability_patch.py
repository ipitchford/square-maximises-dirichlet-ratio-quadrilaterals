"""Portability patch of 2 October 2026: replace the producer's absolute package paths in the dart verifier, its launcher and the
separately written receipt / coverage checkers by paths derived from each file's own location, and make the interpreter of the
launchers configurable.  No predicate, certificate, tolerance or control flow changes.  Every replacement must match exactly once
(otherwise abort with nothing written).  Backups: <file>.bak_20261002_portability.  The unified diff is written beside this script."""
import os, sys, shutil, difflib
ROOT = "/Users/admin/claude-sandbox/quad-ppw-2026-09-30"
OLD = "/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts"
HERE_NEW = "HERE = os.path.dirname(os.path.realpath(__file__))"
P = {
 "darts/replay3.py":  [(f"HERE = '{OLD}'", HERE_NEW)],
 "darts/replay_ok.py": [(f"HERE = '{OLD}'", HERE_NEW)],
 "darts/replay2.py":  [(f"HERE = '{OLD}'", HERE_NEW)],
 "darts/merge_pool.py": [(f"HERE = '{OLD}'", HERE_NEW)],
 "darts/check_charts.py": [(f"sys.path.insert(0, '{OLD}')", "sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))")],
 "darts/manifests_legacy.py": [("import sys, json, itertools, numpy as np\n", "import sys, os, json, itertools, numpy as np\n"),
                               (f"sys.path.insert(0, '{OLD}')", "sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))")],
 "darts/run_fam.py": [("import sys, json, time, math, itertools, numpy as np\n", "import sys, os, json, time, math, itertools, numpy as np\n"),
                      (f"sys.path.insert(0, '{OLD}')", "sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))"),
                      (f"open(f'{OLD}/manifests/{{name.lower()}}.json')",
                       "open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'manifests', f'{name.lower()}.json'))")],
 "darts/accounting.py": [("import sys, json, math\n", "import sys, os, json, math\n"),
                         (f"sys.path.insert(0, '{OLD}')", "sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))"),
                         (f"D = '{OLD}/'", "D = os.path.dirname(os.path.realpath(__file__)) + '/'")],
 "darts/replay_all.sh": [
   ("# Launcher of the aggregate fail-closed proof replay replay_all.py (see its docstring).  The package directory and the interpreter\n"
    "# are ABSOLUTE constants (never derived from $0, so a symlink or a copy of this script cannot select another driver); the driver\n",
    "# Launcher of the aggregate fail-closed proof replay replay_all.py (see its docstring).  The package directory is the physical\n"
    "# directory of this launcher (portability patch of 2026-10-02; before it, an absolute constant of the producer's machine); the\n"
    "# interpreter is PY (an absolute path to an executable; default: the producer's interpreter).  The driver\n"),
   ("#   PIN     the reviewed pin list (default <package>/review/PINNED.sha256)   CERTS_DIR  re-certified family logs (default certs4)\n",
    "#   PIN     the reviewed pin list (default <package>/review/PINNED.sha256)   CERTS_DIR  re-certified family logs (default certs4)\n"
    "#   PY      absolute path of the Python interpreter (default /Users/admin/.venvs/analysis-py313/bin/python)\n"),
   (f"D={OLD}\nPY=/Users/admin/.venvs/analysis-py313/bin/python\n",
    "D=$(cd \"$(/usr/bin/dirname \"$0\")\" && /bin/pwd -P) || exit 2\n"
    "PY=\"${PY:-/Users/admin/.venvs/analysis-py313/bin/python}\"\n"
    "case \"$PY\" in /*) ;; *) echo \"PY must be an absolute path, got '$PY'\"; exit 2 ;; esac\n"
    "[ -x \"$PY\" ] || { echo \"PY is not an executable file: $PY\"; exit 2; }\n")],
 "darts/test_replay_infra.py": [("PY = '/Users/admin/.venvs/analysis-py313/bin/python'", "PY = os.environ.get('PY') or sys.executable")],
 "verify/xreceipt.py": [
   ("Usage: xreceipt.py RECEIPT --pin-sha256 HEX [--mode proof|quick] [--driver-inv HEX] [--shard-inv HEX]   exit 0 iff XRECEIPT PASS\n",
    "Usage: xreceipt.py RECEIPT --pin-sha256 HEX [--mode proof|quick] [--driver-inv HEX] [--shard-inv HEX] [--darts DIR] [--home DIR]\n"
    "  --darts: the package directory the receipt was produced in (default: ../darts beside this file, physical path);\n"
    "  --home: the HOME of the run (default: $HOME).                                            exit 0 iff XRECEIPT PASS\n"),
   (f'DARTS = "{OLD}"\n', 'DARTS = os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "darts"))\n'),
   ('    ap.add_argument("--driver-inv"); ap.add_argument("--shard-inv")\n    a = ap.parse_args()\n',
    '    ap.add_argument("--driver-inv"); ap.add_argument("--shard-inv")\n'
    '    ap.add_argument("--darts", default=DARTS); ap.add_argument("--home", default=os.environ.get("HOME", ""))\n'
    '    a = ap.parse_args()\n    globals()["DARTS"] = os.path.realpath(a.darts)\n'),
   ('"HOME": "/Users/admin", "TMPDIR"', '"HOME": a.home, "TMPDIR"')],
 "verify/test_xreceipt.py": [('PY = "/Users/admin/.venvs/analysis-py313/bin/python"', 'PY = os.environ.get("PY") or sys.executable')],
 "verify/xdart.py": [
   (f'    ix = LeafIndex("{OLD}/certs2/kf_spike.jsonl")', '    ix = LeafIndex(DARTS + "certs2/kf_spike.jsonl")'),
   (f'DARTS = "{OLD}/"\n', 'DARTS = os.environ.get("XDART_DARTS") or os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "darts"))\nDARTS = DARTS.rstrip("/") + "/"\n')],
}
new = {}
for rel, reps in P.items():
    path = os.path.join(ROOT, rel); s = open(path).read(); t = s
    for a, b in reps:
        n = t.count(a)
        if n != 1: sys.exit(f"ABORT (nothing written): {rel}: expected 1 match, found {n}: {a[:70]!r}")
        t = t.replace(a, b)
    new[rel] = (s, t)
# xdart.py: DARTS is used by the --test-index block above its definition; move the definition check
s, t = new["verify/xdart.py"]
if t.index('LeafIndex(DARTS + "certs2/kf_spike.jsonl")') < t.index("DARTS = os.environ.get"):
    # the --test-index block runs at import time before DARTS is defined: define DARTS early instead
    t = t.replace('DARTS = os.environ.get("XDART_DARTS") or os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "darts"))\nDARTS = DARTS.rstrip("/") + "/"\n', "", 1)
    head_anchor = 'if __name__ == "__main__" and "--test-index" in __import__("sys").argv:\n'
    assert t.count(head_anchor) == 1
    t = t.replace(head_anchor, 'import os\nDARTS = os.environ.get("XDART_DARTS") or os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "darts"))\n'
                  'DARTS = DARTS.rstrip("/") + "/"\n\n' + head_anchor, 1)
    new["verify/xdart.py"] = (s, t)
diff = []
for rel, (s, t) in new.items():
    shutil.copy2(os.path.join(ROOT, rel), os.path.join(ROOT, rel) + ".bak_20261002_portability")
    open(os.path.join(ROOT, rel), "w").write(t)
    diff += difflib.unified_diff(s.splitlines(True), t.splitlines(True), fromfile="a/" + rel, tofile="b/" + rel, n=1)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "PORTABILITY_PATCH_2026-10-02.diff")
open(out, "w").write("".join(diff))
print(f"patched {len(new)} files; diff {out} ({sum(1 for l in diff if l.startswith(('+','-')) and not l.startswith(('+++','---')))} changed lines)")
