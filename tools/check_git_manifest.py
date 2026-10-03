"""Check a Git checkout of this repository against the archive manifest.  Author: Anonymous.  Licence: MIT.
The release archive (GitHub release asset and Zenodo record) holds every file listed in MANIFEST.sha256. The Git repository holds all
of them except the bulk files listed in LARGE_FILES.md (certificate logs, chart caches, replay run directories). This script checks that
every file present matches the manifest, that every absent file is listed in LARGE_FILES.md, and that the checkout has no unlisted
files.  Usage: python tools/check_git_manifest.py [--archive]   (--archive: no file may be absent)   exit 0 iff consistent"""
import hashlib, os, sys

root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
os.chdir(root)
manifest = {}
for line in open("MANIFEST.sha256", encoding="utf-8"):
    if line.strip():
        h, f = line.rstrip("\n").split("  ", 1); manifest[f] = h
large = set()
for line in open("LARGE_FILES.md", encoding="utf-8"):
    if line.startswith("| `"): large.add(line.split("`")[1])
present = set()
for dp, dn, fn in os.walk("."):
    dn[:] = [d for d in dn if d != ".git"]
    for f in fn:
        p = os.path.relpath(os.path.join(dp, f), ".")
        if p != "MANIFEST.sha256": present.add(p)
problems = []
for f, h in sorted(manifest.items()):
    if f in present:
        with open(f, "rb") as fh:
            if hashlib.sha256(fh.read()).hexdigest() != h: problems.append(f"changed: {f}")
    elif f not in large or "--archive" in sys.argv: problems.append(f"missing: {f}")
problems += [f"not in the manifest: {f}" for f in sorted(present - set(manifest))]
problems += [f"LARGE_FILES.md lists a file not in the manifest: {f}" for f in sorted(large - set(manifest))]
n_abs = len([f for f in manifest if f not in present])
print(f"manifest: {len(manifest)} files; present {len(present & set(manifest))}; absent (archive only) {n_abs}; problems {len(problems)}")
for p in problems[:50]: print("PROBLEM:", p)
sys.exit(1 if problems else 0)
