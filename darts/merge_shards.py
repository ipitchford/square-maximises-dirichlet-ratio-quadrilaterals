"""Concatenate the shard logs of a cover5 family run (certs3/<fam>.<i>.jsonl, i = 0..N-1) into certs3/<fam>.merged.jsonl after
checking that each shard finished (its .out ends with 'EXIT 0' and a DONE line) and that no root appears in two shards.
Usage: merge_shards.py FAMILY NSHARDS"""
import sys, json, os
fam, n = sys.argv[1].lower(), int(sys.argv[2])
D = '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts/certs3'
roots = {}; lines = []
ok = True
for i in range(n):
    out = open(f'{D}/{fam}.{i}.out').read().strip().splitlines()
    if not out or out[-1] != 'EXIT 0' or not any(l.startswith('DONE') for l in out): print(f'shard {i} not finished'); ok = False
    for l in open(f'{D}/{fam}.{i}.jsonl'):
        lines.append(l if l.endswith('\n') else l + '\n')
        try: d = json.loads(l)
        except Exception: continue
        r = d.get('root')
        if r is not None:
            if roots.setdefault(r, i) != i: print('root in two shards:', r); ok = False
if not ok: sys.exit(1)
with open(f'{D}/{fam}.merged.jsonl', 'w') as f: f.writelines(lines)
print(fam, 'merged', len(lines), 'lines', len(roots), 'roots')
