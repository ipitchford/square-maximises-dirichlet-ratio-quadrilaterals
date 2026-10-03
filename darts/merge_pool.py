"""Merge all logs of a family (static shard logs certs3/<fam>.<i>.jsonl and pool worker logs certs3/<fam>.w<ID>.jsonl) into
certs3/<fam>.merged.jsonl.  Fails (exit 1) if a (root, path) tree record (leaf / split / esc), handoff or item_done occurs twice
across the files, if a root has more than one root_done, or if a file ends with a torn line (a root may be spread over several
worker logs when it was split into work items).  Completeness, duplicates (root, path), schemas and tiling are checked by replay3.py on the merge.
Usage: merge_pool.py FAMILY [--dir DIR]   (DIR: certs3 (default) or the re-certification directory certs4)"""
import sys, os, glob, json
HERE = os.path.dirname(os.path.realpath(__file__))
fam = sys.argv[1].lower()
DIR = sys.argv[sys.argv.index('--dir') + 1] if '--dir' in sys.argv else 'certs3'
DIR = DIR if os.path.isabs(DIR) else os.path.join(HERE, DIR)
files = sorted(f for f in glob.glob(os.path.join(DIR, fam + '.*.jsonl')) if not f.endswith('.merged.jsonl'))
where = {}; ndone = {}; bad = []; nlines = 0
out_lines = []
for f in files:
    raw = open(f, 'rb').read()
    if raw and not raw.endswith(b'\n'): bad.append(f'{f}: torn last line')
    for i, line in enumerate(raw.decode().split('\n')):
        if not line.strip(): continue
        try: d = json.loads(line)
        except Exception: bad.append(f'{f}:{i + 1}: not JSON'); continue
        rid = d.get('root')
        if rid is None: bad.append(f'{f}:{i + 1}: record without root'); continue
        # a root may be spread over several worker logs (work items); a (root, path, record class) must occur once overall
        t = d.get('type'); cls = 'tree' if t in ('leaf', 'split', 'esc') else t
        key = (rid, d.get('path'), cls)
        if t != 'root_done':
            if key in where: bad.append(f'{f}:{i + 1}: duplicate {cls} record {rid} {d.get("path")} (also in {where[key]})')
            where[key] = f
        if t == 'root_done': ndone[rid] = ndone.get(rid, 0) + 1
        out_lines.append(line); nlines += 1
for rid, c in ndone.items():
    if c != 1: bad.append(f'root {rid}: {c} root_done records')
if bad:
    for b in bad[:20]: print('MERGE FAIL', b)
    sys.exit(1)
with open(os.path.join(DIR, fam + '.merged.jsonl'), 'w') as fo:
    for line in out_lines: fo.write(line + '\n')
print(json.dumps(dict(family=fam, files=[os.path.basename(f) for f in files], lines=nlines, roots_done=len(ndone))))
