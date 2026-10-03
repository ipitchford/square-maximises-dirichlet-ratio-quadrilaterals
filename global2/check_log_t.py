"""Structural check of a transport cover log (cover_t.jl JSONL), independent of the Julia code.

usage: python check_log_t.py <log.jsonl> [<prefix_backup.jsonl>]

* every line is valid JSON with a path and a status;
* duplicate paths: reported, with a flag if their statuses conflict;
* the tree is rebuilt from the 24 roots of root_boxes_t (a-priori box split at a, b, c = S and d = ±S, S = Float64(1/5),
  bounds rounded outward exactly as cover.jl) with the children recomputed in IEEE double exactly as children();
  every node reached must have a record; a split record must have dim in 1..4 and both children are then visited;
  a leaf record must carry lo/hi equal (bitwise) to the recomputed box; every record must be reached exactly once
  (no orphans); no FAIL / ERROR / pilot_pass leaves;
* optional: the backup log must be a line prefix of the log (crash/resume) and no path recorded in the backup may be
  recorded again after it.
"""
import json, sys, math, collections
from fractions import Fraction

def f64_up(q):
    x = float(q)
    return x if Fraction(x) >= q else math.nextafter(x, math.inf)

logf = sys.argv[1]
recs = collections.defaultdict(list); nbad = 0; nlines = 0
for k, ln in enumerate(open(logf), 1):
    nlines += 1
    try:
        r = json.loads(ln)
        assert isinstance(r.get('path'), str) and isinstance(r.get('status'), str)
    except Exception as e:
        nbad += 1; print(f"line {k}: unparsable ({e}): {ln[:120]!r}"); continue
    recs[r['path']].append((k, r))
dups = {p: v for p, v in recs.items() if len(v) > 1}
conflict = {p: v for p, v in dups.items() if len({(r['status'], r.get('dim')) for _, r in v}) > 1}
print(f"lines {nlines}, unparsable {nbad}, distinct paths {len(recs)}, duplicate paths {len(dups)} (conflicting {len(conflict)})")
for p, v in list(conflict.items())[:10]:
    print("  CONFLICT", p, [(k, r['status']) for k, r in v])

S = 0.2                                    # Float64(1//5) (round to nearest)
assert Fraction(S) != Fraction(1, 5) and abs(Fraction(S) - Fraction(1, 5)) < Fraction(1, 10**16)
A = f64_up(Fraction(905, 1000)); B = f64_up(Fraction(436, 100))
lo0 = [0.0, 0.0, 0.0, -A]; hi0 = [A, B, B, A]
def cut(l, h, pts):
    xs = [l] + [x for x in pts if l < x < h] + [h]
    return [(xs[i], xs[i + 1]) for i in range(len(xs) - 1)]
ivs = [cut(lo0[i], hi0[i], [S]) for i in range(3)] + [cut(lo0[3], hi0[3], [-S, S])]
roots = []; k = 0
for ia in ivs[0]:
    for ib in ivs[1]:
        for ic in ivs[2]:
            for id_ in ivs[3]:
                k += 1
                roots.append((f"T{k}", [ia[0], ib[0], ic[0], id_[0]], [ia[1], ib[1], ic[1], id_[1]]))
assert len(roots) == 24
DIM = "abcd"
visited = collections.Counter(); leaves = collections.Counter(); missing = []; mismatch = []; badleaf = []; vol = collections.Counter()
stack = list(roots); depthmax = 0
while stack:
    path, lo, hi = stack.pop()
    v = recs.get(path)
    if not v:
        missing.append(path); continue
    visited[path] += 1
    r = v[-1][1]
    if r['status'] == 'split':
        d = int(r['dim']) - 1
        assert 0 <= d < 4, (path, r)
        mid = (lo[d] + hi[d]) / 2
        if not (lo[d] < mid < hi[d]):
            mismatch.append(path + " (degenerate split)"); continue
        h1 = list(hi); h1[d] = mid; l2 = list(lo); l2[d] = mid
        stack.append((path + "." + DIM[d] + "0", list(lo), h1)); stack.append((path + "." + DIM[d] + "1", l2, list(hi)))
        depthmax = max(depthmax, path.count('.') + 1)
    else:
        if r.get('lo') != lo or r.get('hi') != hi:
            mismatch.append(path)
        leaves[r['status']] += 1
        w = 1.0
        for i in range(4): w *= hi[i] - lo[i]
        vol[r['status']] += w
        if r['status'] in ('FAIL', 'ERROR', 'pilot_pass'):
            badleaf.append((path, r['status']))
orphans = [p for p in recs if visited[p] == 0]
twice = [p for p, c in visited.items() if c > 1]
print(f"tree: {len(roots)} roots, {sum(leaves.values())} leaves, max depth {depthmax}; missing records {len(missing)}, "
      f"lo/hi mismatches {len(mismatch)}, orphan records {len(orphans)}, nodes reached twice {len(twice)}, FAIL/ERROR/pilot leaves {len(badleaf)}")
print("leaves by status:", dict(leaves.most_common()))
print("leaf volume by status:", {s: float(f"{x:.6e}") for s, x in vol.items()}, " total %.6e" % sum(vol.values()))
tot = 1.0
for i in range(4): tot *= hi0[i] - lo0[i]
print("a-priori box volume %.6e; |leaf volume sum - box volume| = %.3e" % (tot, abs(sum(vol.values()) - tot)))
for p in (missing[:5] + mismatch[:5] + orphans[:5]):
    print("  PROBLEM", p)
ok = nbad == 0 and not conflict and not missing and not mismatch and not orphans and not twice and not badleaf
if len(sys.argv) > 2:
    bak = sys.argv[2]
    bl = open(bak).read().splitlines(); ml = open(logf).read().splitlines()
    pref = ml[:len(bl)] == bl
    bpaths = {json.loads(x)['path'] for x in bl}
    again = [json.loads(x)['path'] for x in ml[len(bl):] if json.loads(x)['path'] in bpaths]
    print(f"backup: {len(bl)} lines, line prefix of the log: {pref}; records after the resume: {len(ml) - len(bl)}, "
          f"of which re-recording a backup path: {len(again)}")
    ok = ok and pref and not again
print("LOG STRUCTURE OK" if ok else "LOG STRUCTURE PROBLEMS")
