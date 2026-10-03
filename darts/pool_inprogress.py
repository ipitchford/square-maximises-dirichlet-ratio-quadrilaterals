"""In-progress pool work: for every unfinished root (some volume uncovered), the elapsed time since its first record... approximated
by the claim file's age, certificates so far (all logs), uncovered volume fraction, and its open work items (claimed / queued).
Usage: pool_inprogress.py"""
import sys, os, json, glob, time, numpy as np
HERE = '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts'
sys.path.insert(0, HERE)
import pool
now = time.time()
for fam in ('PP', 'SP_spk', 'SP_col', 'PF', 'FP'):
    man = {rid: (p0, r) for rid, p0, r in json.load(open(os.path.join(HERE, 'manifests', fam.lower() + '.json')))['roots']}
    roots_done, items_done = pool._scan(fam)
    cd = os.path.join(HERE, 'claims', fam)
    if not os.path.isdir(cd): continue
    claimed = [x for x in os.listdir(cd) if not x.startswith('.')]
    open_roots = sorted({rid for rid in man if rid in claimed and rid not in roots_done})
    if not open_roots: continue
    vol = {rid: 0.0 for rid in open_roots}; certs = {rid: 0 for rid in open_roots}
    for f in glob.glob(os.path.join(HERE, 'certs3', f'{fam.lower()}.*.jsonl')):
        if f.endswith('.merged.jsonl'): continue
        for line in open(f):
            if '"leaf"' not in line: continue
            d = json.loads(line)
            if d.get('type') == 'leaf' and d['root'] in vol:
                vol[d['root']] += float(np.prod([2*v for v in d['r'] if v > 0])); certs[d['root']] += d['status'] == 'cert'
    queue = pool.read_queue(fam)
    for rid in open_roots:
        V = float(np.prod([2*v for v in man[rid][1] if v > 0]))
        unf = 1 - vol[rid]/V
        if unf <= 1e-12: continue
        its = [q for q in queue if q['root'] == rid and (rid, q['path']) not in items_done]
        ncl = sum(1 for q in its if os.path.exists(pool.claim_path(fam, q['path'])))
        age = (now - os.path.getmtime(pool.claim_path(fam, rid)))/60
        print(f'{fam:6s} {rid:16s} claimed {age:6.1f} min ago  certs {certs[rid]:5d}  uncovered volume {unf:.3f}  '
              f'open items {len(its)} ({ncl} claimed)')
