"""Negative controls for verify/xlocal.py.  Author: Anonymous.
The unmodified record must pass, and each of five corrupted copies must be rejected:
a flipped sign of one exact bound, a deleted direction box, a duplicated box, a shift level raised above its admissible value,
and a missing s-interval.  The copies are written to a temporary directory; the input file is never modified.
Usage: xlocal_negative_controls.py [path/to/schur_cover2.json]      exit 0 iff the base passes and every corruption is rejected"""
import sys, os, json, copy, subprocess, tempfile

HERE = os.path.dirname(os.path.realpath(__file__))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "local", "replay_v2", "schur_cover2.json")
XL = os.path.join(HERE, "xlocal.py")

def run(path):
    r = subprocess.run([sys.executable, XL, path], capture_output=True, text=True)
    return r.returncode, [l.strip() for l in r.stdout.splitlines() if "XLOCAL" in l or l.startswith("   ")][-2:]

d = json.load(open(SRC))
def flip(x):
    u = x["boxes"][5]["used"][3]; u[2] = u[2].lstrip("-"); u[3] = -u[3]
def drop(x): del x["boxes"][17]
def dup(x): x["boxes"].append(copy.deepcopy(x["boxes"][3]))
def raise_level(x): x["boxes"][0]["used"][9][0] = "49"
def drop_interval(x): x["boxes"][7]["used"] = x["boxes"][7]["used"][:9]
ok = True
with tempfile.TemporaryDirectory() as tmp:
    rc, msg = run(SRC)
    print(f"{'ok  ' if rc == 0 else 'FAIL'} unmodified record passes (exit {rc})"); ok &= rc == 0
    for name, f in (("flipped sign", flip), ("deleted box", drop), ("duplicated box", dup), ("raised level", raise_level),
                    ("missing s-interval", drop_interval)):
        x = copy.deepcopy(d); f(x); p = os.path.join(tmp, name.replace(" ", "_") + ".json"); json.dump(x, open(p, "w"))
        rc, msg = run(p)
        print(f"{'ok  ' if rc != 0 else 'FAIL'} {name}: rejected={rc != 0} | {' / '.join(msg)[:200]}"); ok &= rc != 0
print("XLOCAL NEGATIVE CONTROLS", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
