"""Remove every record of the given roots from the pool worker logs (certs3/<fam>.w<ID>.jsonl), their queued work items and all
their claims, so that the roots are processed again from scratch (used after a cover-heuristic change for roots that hit a LOUD failure).
ONLY run with all pool workers stopped.  The removed lines are kept in pilot_runs/stripped/<root>.jsonl.
Usage: pool_strip.py FAM ROOT [ROOT ...]      or      pool_strip.py FAM --failed   (all roots with a FAIL record)"""
import sys, os, glob, json, subprocess
HERE = '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts'
fam = sys.argv[1]
if subprocess.run(['pgrep', '-f', 'pool.py [0-9]'], capture_output=True).stdout.strip():
    sys.exit('pool workers are running: stop them first')
files = sorted(glob.glob(os.path.join(HERE, 'certs3', fam.lower() + '.w*.jsonl')))
if sys.argv[2] == '--failed':
    roots = set()
    for f in files:
        for line in open(f):
            if '"FAIL"' in line:
                d = json.loads(line)
                if d.get('status') == 'FAIL' or (d.get('type') == 'root_done' and d.get('ok') is False): roots.add(d['root'])
else:
    roots = set(sys.argv[2:])
os.makedirs(os.path.join(HERE, 'pilot_runs', 'stripped'), exist_ok=True)
for f in files:
    keep, drop = [], []
    for line in open(f):
        if not line.strip(): continue
        d = json.loads(line)
        (drop if d.get('root') in roots else keep).append(line if line.endswith('\n') else line + '\n')
    if drop:
        with open(f + '.tmp', 'w') as fo: fo.writelines(keep)
        os.replace(f + '.tmp', f)
        for line in drop:
            rid = json.loads(line)['root']
            with open(os.path.join(HERE, 'pilot_runs', 'stripped', rid + '.jsonl'), 'a') as fo: fo.write(line)
        print(f'{os.path.basename(f)}: removed {len(drop)} lines')
for rid in sorted(roots):
    c = os.path.join(HERE, 'claims', fam, rid)
    if os.path.exists(c): os.unlink(c)
    print('released', rid)
# work items of these roots (queue files carry their root id: root ids may be string prefixes of other roots' paths)
qd = os.path.join(HERE, 'queue', fam)
if os.path.isdir(qd):
    nq = 0
    for fn in os.listdir(qd):
        if not fn.endswith('.json'): continue
        try: it = json.load(open(os.path.join(qd, fn)))
        except Exception: continue
        if it.get('root') in roots:
            os.unlink(os.path.join(qd, fn)); nq += 1
            c = os.path.join(HERE, 'claims', fam, it['path'])
            if os.path.exists(c): os.unlink(c)
    print(f'removed {nq} queue items (and their claims)')
