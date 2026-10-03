"""Dynamic worker pool for the cover5 family runs: load-balanced, crash-safe, roots in a fixed pseudo-random order, heavy roots
split into independent work items.
Usage:  pool.py WORKER_ID FAM [FAM ...]          (all workers should be given the same family list, in the same order)

* Every worker keeps one log per family, certs3/<fam>.w<ID>.jsonl (cover5 protocol: resumable at leaf level).
* Work items: the manifest roots (item id = root id) and sub-boxes handed off by other items (item id = the sub-box path;
  queue/<fam>/<path>.json holds root, path, box, depth, mesh level and split hint).  An item is processed for at most ITEM_S
  seconds; after that every pending sub-box of it is handed to the queue (log record 'handoff' first, then the queue file), so
  that all workers share heavy roots and progress is measurable.  Queue items are taken before new roots.
* Items are claimed atomically: the claim file claims/<fam>/<item id> is created by os.link from a temporary file that already
  holds the owner id (link fails if the claim exists), so a claim always names its owner.
* On start a worker re-creates missing queue files for its own handoff records, resumes the items it has claimed but not finished
  (their partial records are in its own log), then takes queue items and new roots (fixed pseudo-random order, seed SEED).
  A root with a root_done record in ANY log certs3/<fam>.*.jsonl is skipped; an item with an item_done record likewise.
* Root lists come from the FROZEN manifests (manifests/<fam>.json), exactly as run_fam.py."""
import sys, os, json, time, glob, numpy as np
HERE = '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts'
sys.path.insert(0, HERE)
WORK = os.environ.get('POOL_WORK', HERE)       # logs, claims, queue, pool_runs (a scratch directory for tests only)
SEED = 20261001
ITEM_S = float(os.environ.get('POOL_ITEM_S', 300.0))
ONLY = os.environ.get('POOL_ONLY')             # tests only: comma-separated root ids to process

def _scan(fam):
    """(done roots with ok root_done, done items {(root, path)}) over all logs of the family"""
    roots, items = set(), set()
    for f in glob.glob(os.path.join(WORK, 'certs3', fam.lower() + '.*.jsonl')):
        if f.endswith('.merged.jsonl'): continue
        with open(f) as fh:
            for line in fh:
                if '"root_done"' in line or '"item_done"' in line:
                    try: d = json.loads(line)
                    except Exception: continue
                    if d.get('type') == 'root_done' and d.get('ok') is True: roots.add(d['root'])
                    elif d.get('type') == 'item_done': items.add((d['root'], d['path']))
    return roots, items

def done_roots(fam): return _scan(fam)[0]

def claim_path(fam, iid): return os.path.join(WORK, 'claims', fam, iid)

def claim(fam, iid, wid):
    d = os.path.join(WORK, 'claims', fam); os.makedirs(d, exist_ok=True)
    tmp = os.path.join(d, f'.tmp.{iid}.w{wid}.{os.getpid()}')
    with open(tmp, 'w') as f: f.write(str(wid)); f.flush(); os.fsync(f.fileno())
    try:
        os.link(tmp, claim_path(fam, iid)); ok = True
    except FileExistsError:
        ok = False
    finally:
        os.unlink(tmp)
    return ok

def owner(fam, iid):
    try: return open(claim_path(fam, iid)).read().strip()
    except FileNotFoundError: return None

def order(fams):
    items = []
    for fam in fams:
        man = json.load(open(os.path.join(HERE, 'manifests', fam.lower() + '.json')))
        items += [(fam, rid, p0, r) for rid, p0, r in man['roots']]
    perm = np.random.default_rng(SEED).permutation(len(items))
    return [items[i] for i in perm]

def qdir(fam): return os.path.join(WORK, 'queue', fam)

def write_queue(fam, root, path, p0, r, depth, k, hint):
    """durable queue entry for a handed-off sub-box (atomic rename)"""
    d = qdir(fam); os.makedirs(d, exist_ok=True)
    h = {kk: (vv.tolist() if hasattr(vv, 'tolist') else vv) for kk, vv in (hint or {}).items()}
    tmp = os.path.join(d, f'.tmp.{path}.{os.getpid()}')
    with open(tmp, 'w') as f:
        json.dump(dict(fam=fam, root=root, path=path, p0=list(p0), r=list(r), depth=depth, k=k, hint=h), f)
        f.flush(); os.fsync(f.fileno())
    os.replace(tmp, os.path.join(d, path + '.json'))

def read_queue(fam):
    d = qdir(fam)
    if not os.path.isdir(d): return []
    out = []
    for fn in sorted(os.listdir(d)):
        if fn.startswith('.') or not fn.endswith('.json'): continue
        try: out.append(json.load(open(os.path.join(d, fn))))
        except Exception: continue
    return out

def single_instance(wid):
    """refuse to run if another live process holds this worker id (two writers on one log would corrupt it)"""
    lk = os.path.join(WORK, 'pool_runs', f'w{wid}.lock'); os.makedirs(os.path.dirname(lk), exist_ok=True)
    if os.path.exists(lk):
        try:
            pid = int(open(lk).read().strip()); os.kill(pid, 0)
            sys.exit(f'worker {wid} already running (pid {pid})')
        except (ValueError, ProcessLookupError, PermissionError):
            pass
    with open(lk, 'w') as f: f.write(str(os.getpid()))

def main():
    wid = int(sys.argv[1]); fams = sys.argv[2:]
    single_instance(wid)
    import fam5, run_fam
    from grid import dyadic_root
    from cover5 import Cover5
    covers = {}
    def cover(fam):
        if fam not in covers:
            os.makedirs(os.path.join(WORK, 'certs3'), exist_ok=True)
            covers[fam] = Cover5(os.path.join(WORK, 'certs3', f'{fam.lower()}.w{wid}.jsonl'), fam5.FAMILIES[fam](),
                                 fem_problem=run_fam.fem_problem(fam))
        return covers[fam]
    roots = order(fams)
    rootmap = {(fam, rid): (p0, r) for fam, rid, p0, r in roots}
    if ONLY: roots = [x for x in roots if x[1] in ONLY.split(',')]
    t0 = time.time(); n = 0
    def run(it, how):
        nonlocal n
        fam = it['fam']; cv = cover(fam); c0 = dict(cv.stats); tr = time.time()
        hf = lambda rt, pth, q, rr, dd, kk, hh: write_queue(fam, rt, pth, q, rr, dd, kk, hh)
        res = cv.process_item(it['root'], it['path'], it['p0'], it['r'], it['depth'], it['k'], it['hint'],
                              deadline=time.time() + ITEM_S, handoff_fn=hf)
        d = {kk: (round(cv.stats[kk] - c0[kk], 1 if kk != 'vol' else 12) if isinstance(cv.stats[kk], float) else cv.stats[kk] - c0[kk])
             for kk in ('cert', 'discard', 'FAIL', 'bo', 'fem', 'inc', 't_bo', 't_inc', 't_fem', 'vol')}
        vol = float(np.prod([2*v for v in it['r'] if v > 0]))
        n += 1
        d['vol_done'] = d.pop('vol')
        print(json.dumps(dict(w=wid, fam=fam, root=it['root'], item=it['path'], how=how, ok=res['ok'], handed=res['handed'],
                              vol=vol, dt=round(time.time() - tr, 2), **d, wall=round(time.time() - t0, 1),
                              t=round(time.time(), 1))), flush=True)
    def root_item(fam, rid):
        p0, r = dyadic_root(*rootmap[(fam, rid)])
        return dict(fam=fam, root=rid, path=rid, p0=list(p0), r=list(r), depth=0, k=4, hint=None)
    # 0. re-create missing queue files for this worker's own handoff records (crash between the log record and the queue file)
    for fam in fams:
        for h in cover(fam).handoff_recs:
            if not os.path.exists(os.path.join(qdir(fam), h['path'] + '.json')):
                write_queue(fam, h['root'], h['path'], h['p0'], h['r'], h['depth'], h['k'], None)
    while True:
        st = {fam: _scan(fam) for fam in fams}
        def done(it): return ((it['path'] == it['root'] and it['root'] in st[it['fam']][0])
                              or (it['root'], it['path']) in st[it['fam']][1])
        # 1. own claimed, unfinished items (roots and queue items)
        mine = [it for fam in fams for it in read_queue(fam) if owner(fam, it['path']) == str(wid) and not done(it)]
        mine += [root_item(fam, rid) for fam, rid, _, _ in roots
                 if owner(fam, rid) == str(wid) and not done(dict(fam=fam, root=rid, path=rid))]
        if mine:
            run(mine[0], 'resume'); continue
        # 2. queue items (sub-boxes of heavy roots) before new roots
        got = None
        for fam in fams:
            for it in read_queue(fam):
                if done(it) or os.path.exists(claim_path(fam, it['path'])): continue
                if claim(fam, it['path'], wid): got = it; break
            if got: break
        if got is None:
            for fam, rid, _, _ in roots:
                if rid in st[fam][0] or (rid, rid) in st[fam][1] or os.path.exists(claim_path(fam, rid)): continue
                if claim(fam, rid, wid): got = root_item(fam, rid); break
        if got is None: break
        run(got, 'item' if got['path'] != got['root'] else 'new')
    for cv in covers.values(): cv.close()
    print('DONE', json.dumps(dict(w=wid, items=n, wall=round(time.time() - t0, 1))), flush=True)

if __name__ == '__main__':
    main()
