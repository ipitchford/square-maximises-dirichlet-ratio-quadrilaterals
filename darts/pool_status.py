"""Progress and ETA of the worker pool (pool.py), by VOLUME: per family, the parameter volume of all manifest roots (active
coordinates), the volume already covered by leaves (all logs certs3/<fam>.*.jsonl), and the recent rate = leaf volume produced
per worker-second by the items finished since EPOCH (pool_runs/w*.out lines: vol_done, dt).  ETA = sum over families of
remaining volume / rate, divided by the number of live workers.  (Cost per volume is not uniform: heavy regions cost more, so the
estimate tracks the mix of the recent work; it is printed per family together with the queue state.)
Usage: pool_status.py [--since=EPOCH]   (default: pool_runs/since_epoch)"""
import sys, os, json, glob, time, subprocess, numpy as np
HERE = '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts'
sys.path.insert(0, HERE)
arg = next((x.split('=')[1] for x in sys.argv[1:] if x.startswith('--since=')), None)
since = float(arg) if arg else float(open(os.path.join(HERE, 'pool_runs', 'since_epoch')).read())
now = time.time()
FAMS = ('PP', 'SP_spk', 'SP_col', 'PF', 'FP')
live = subprocess.run(['pgrep', '-f', 'pool.py [0-9]'], capture_output=True, text=True).stdout.split()
rate_v = {f: 0.0 for f in FAMS}; rate_t = {f: 0.0 for f in FAMS}; nitems = 0
for f in glob.glob(os.path.join(HERE, 'pool_runs', 'w*.out')):
    for l in open(f):
        if not l.startswith('{'): continue
        d = json.loads(l)
        if 'vol_done' not in d or d['t'] - d['dt'] < since: continue
        rate_v[d['fam']] += d['vol_done']; rate_t[d['fam']] += d['dt']; nitems += 1
tot_ws = 0.0; rows = []
for fam in FAMS:
    man = json.load(open(os.path.join(HERE, 'manifests', fam.lower() + '.json')))
    V = sum(float(np.prod([2*v for v in r if v > 0])) for _, _, r in man['roots'])
    C = 0.0
    for f in glob.glob(os.path.join(HERE, 'certs3', fam.lower() + '.*.jsonl')):
        if f.endswith('.merged.jsonl'): continue
        for line in open(f):
            if '"leaf"' not in line: continue
            d = json.loads(line)
            if d.get('type') == 'leaf': C += float(np.prod([2*v for v in d['r'] if v > 0]))
    q = os.path.join(HERE, 'queue', fam); nq = len([x for x in os.listdir(q) if x.endswith('.json')]) if os.path.isdir(q) else 0
    rem = max(V - C, 0.0)
    rt = rate_v[fam]/rate_t[fam] if rate_t[fam] > 0 and rate_v[fam] > 0 else None
    ws = rem/rt if rt else None
    if ws is not None: tot_ws += ws
    rows.append((fam, V, C, rem, rt, ws, nq, rate_t[fam]))
print(time.strftime('%Y-%m-%d %H:%M:%S'), f'pool status (rates from items finished since {time.strftime("%H:%M", time.localtime(since))};'
      f' {nitems} items; {len(live)} live workers)')
for fam, V, C, rem, rt, ws, nq, wt in rows:
    print(f'  {fam:7s} covered {C/V:6.1%} of volume, remaining {rem/V:6.1%};  measured on {wt/3600:5.2f} worker-h: '
          + (f'{ws/3600:6.1f} worker-h remaining' if ws is not None else 'no rate yet') + f';  queue files {nq}')
if live:
    missing = [r[0] for r in rows if r[5] is None]
    print(f'  ETA (families with a rate): {tot_ws/3600:.1f} worker-h = {tot_ws/3600/len(live):.1f} h wall at {len(live)} workers'
          + (f'; no rate yet for {missing}' if missing else ''))
