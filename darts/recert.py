"""Re-certification of the cover5 family logs with the corrected verifier (adversarial review 2026-10-02, repair 6).

Input : certs3/<fam>.merged.jsonl (the certified tree of every root) and manifests/<fam>.json (exact grid encodings).
Output: certs4/<fam>.w<ID>.jsonl (one file per worker), later merged by merge_pool.py --dir certs4 and replayed by replay3.py.
Launched by recert_all.sh (worker count N recorded in certs4/RECERT_N at the first start; any other N, a worker ID outside 1..N, a
claim naming no worker 1..N or a second process with the same worker ID is REFUSED, so no claim can be orphaned by a restart).
For every root (claimed dynamically: certs4/claims/<fam>/<root>, created atomically by link(2) of a file that already holds the
worker ID; a restarted worker resumes its own claims):
  * the structural records of the old tree (split, esc) are copied, so the root keeps its exact tiling;
  * every old leaf is recomputed from scratch with the corrected code (recert_leaf):
      discard : the recorded reason re-checked with the current predicate (else any other family reason);
      BO2     : the recorded slicing axis, cell count, arm mode and inner-dart index; fresh rigorous 1-D bound E, fresh inner dart,
                fresh P2 test at the recorded level (then the next two levels);
      INC     : the recorded CR / P2 levels; fresh outer / inner apexes (inc.inc_darts), inclusion check, fresh L1 lower bound at
                the outer dart, fresh P2 test of the inner dart (then the next levels);
      FEM     : the recorded mesh level; fresh chart setup (BD chart on the exact Cartesian enclosure, or FARL), solve, sigma2,
                certificate;
    a success writes a new leaf record at the same path; a failure hands the leaf box to Cover5.process_item (the full cover
    search: certificates at the box first, then bisection), whose records go to the same output;
  * then a fresh root_done record (ok iff no FAIL leaf), box = the dyadic snapping of the manifest box.
Resumable: the output is read on restart (a torn last line, left by a kill during a write, is cut off first); finished leaves
(new leaf record, or item_done of the leaf path) are skipped, an unfinished subdivision resumes inside Cover5 (it reads the same
output), copied structural records are not duplicated.  Completeness is checked by certs_complete.py before any replay.
Usage: recert.py --split FAMILY                 (once per family: per-root files of the old records, certs4/old/<fam>/)
       recert.py FAMILY WORKER_ID                (production worker; exits when no unclaimed root is left)
       recert.py --time FAMILY N [SEED]          (test: recompute N random leaves per method, report times and failures)
       recert.py --code-id                       (sha256 of the re-certification code: recert_all.sh records it at the first start)"""
import sys, os, json, time, random, collections, re, fcntl, hashlib
HERE = '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts'
sys.path.insert(0, HERE)
import numpy as np
import fam5, run_fam
from grid import dyadic_root
from nanguard import exact_box

class Ctx:
    def __init__(self, name):
        self.name = name; self.fam = fam5.FAMILIES[name](); self.gcs = {}
    def gc_for(self, k):
        from gcert import GCert
        if k not in self.gcs: self.gcs[k] = GCert(run_fam.fem_problem(self.name)(k))
        return self.gcs[k]
    def gc_setup(self, pb):
        from gcert import GCert
        key = ('setup', id(pb))
        if key not in self.gcs: self.gcs[key] = (pb, GCert(pb))
        return self.gcs[key][1]

def recert_leaf(cx, L):
    """recompute one leaf record; returns (new record or None, info)"""
    fam = cx.fam
    p0 = np.array(L['p0'], float); r = np.array(L['r'], float)
    exact_box(p0, r)
    base = dict(type='leaf', path=L['path'], root=L['root'], status=L['status'], p0=L['p0'], r=L['r'])
    if L['status'] == 'discard':
        if fam.check(L['reason'], p0, r): return dict(base, reason=L['reason']), 'discard'
        why = fam.discard(p0, r)
        if why: return dict(base, reason=why), 'discard (new reason)'
        return None, 'discard fails'
    m = L['method']; base['method'] = m
    if m == 'BO2':
        from bo2 import bo_family, arm_variant
        from p2fix import p2_cert
        armA, armC = arm_variant(*fam.arms(p0, r), L.get('arms', 'polar'))
        E, _ = bo_family(armA, armC, L['theta'], L['K'])
        if not E > 0: return None, 'BO2: E <= 0'
        inn = fam.inners(p0, r)
        if not (0 <= L['inner'] < len(inn)): return None, 'BO2: inner index'
        Ai, Ci = inn[L['inner']]
        for k in range(L['k'], 7):
            res = p2_cert(Ai, Ci, E, k)
            if res.get('ok'):
                rec = dict(base, k=k, E=E, theta=L['theta'], K=L['K'], inner=L['inner'], Ain=list(Ai) if Ai is not None else None,
                           Cin=list(Ci), schur=res['schur'], lam2p2=res['lam2p2'])
                if 'arms' in L: rec['arms'] = L['arms']
                return rec, 'BO2'
        return None, 'BO2: P2 ' + str(res.get('why'))
    if m == 'INC':
        import inc
        from p2fix import p2_cert
        dd = inc.inc_darts(p0, r)
        if dd is None: return None, 'INC: no apexes'
        Ao, Co, Ai, Ci = dd
        if not inc.verify_inclusion(p0, r, Ao, Co, Ai, Ci): return None, 'INC: inclusion'
        E, info = inc.l1_point(Ao, Co, L['k'])
        if E is None: return None, 'INC: L1 ' + str(info.get('why'))
        for kp in range(L['kp'], 7):
            res = p2_cert(Ai, Ci, E, kp)
            if res.get('ok'):
                return dict(base, k=L['k'], kp=kp, Aout=list(Ao), Cout=list(Co), Ain=list(Ai), Cin=list(Ci), E=E, sigma=info['sigma'],
                            schur=res['schur'], lam2p2=res['lam2p2'], lam1h=info['lam1h'], h=info['h'], temple=info['temple'],
                            ratio=L.get('ratio', 0.0)), 'INC'
        return None, 'INC: P2 ' + str(res.get('why'))
    if m == 'FEM':
        from gsigma import sigma_box
        k = L['k']
        if hasattr(fam, 'fem_setup'):
            st = fam.fem_setup(p0, r, k)
            if st is None: return None, 'FEM: setup'
            pb, pc, rc = st; gc = cx.gc_setup(pb)
        else:
            cb = fam.cart_box(p0, r)
            if cb is None: return None, 'FEM: cart box'
            pc, rc = cb; gc = cx.gc_for(k)
        pc = np.array(pc, float); rc = np.array(rc, float)
        sol = gc.solve(pc, rc)
        sc, sinfo = sigma_box(gc, sol, pc, rc)
        if sc is None: return None, 'FEM: sigma ' + str(sinfo.get('why'))
        res = gc.certify(pc, rc, sc, sol=sol)
        if not res.get('ok'): return None, 'FEM: ' + str(res.get('why'))
        return dict(base, k=k, pc=pc.tolist(), rc=rc.tolist(), sigma=[pc.tolist(), sc, k, sinfo.get('frac')], L1=list(res['L1']),
                    schur=res['schur'], h=res['h'], temple=res['temple'], tails=res.get('tails')), 'FEM'
    return None, 'unknown method'

def split_old(name):
    """one-time: certs4/old/<fam>/<root>.jsonl with the leaf / split / esc records of each root (in log order), so that a worker
    reads only the roots it claims (the merged PP log is 170 MB)"""
    final = os.path.join(os.environ.get('RECERT_DIR', os.path.join(HERE, 'certs4')), 'old', name)
    if os.path.exists(os.path.join(final, '_COMPLETE')): return final
    if os.path.exists(final): raise SystemExit(f'{final} exists without _COMPLETE: an older interrupted split; move it away first')
    os.makedirs(os.path.dirname(final), exist_ok=True)
    d = final + f'.partial.{os.getpid()}'                   # built here, renamed into place when complete (an interrupted split
    os.makedirs(d)                                          # leaves only an ignored .partial directory, never a half-filled one)
    by = collections.defaultdict(list)
    with open(os.path.join(HERE, 'certs3', f'{name.lower()}.merged.jsonl')) as f:
        for line in f:
            t = json.loads(line).get('type')
            if t in ('leaf', 'split', 'esc'): by[json.loads(line)['root']].append(line)
            if sum(len(v) for v in by.values()) > 200000:
                for rid, lines in by.items():
                    with open(os.path.join(d, rid + '.jsonl'), 'a') as fo: fo.writelines(lines)
                by.clear()
    for rid, lines in by.items():
        with open(os.path.join(d, rid + '.jsonl'), 'a') as fo: fo.writelines(lines)
    open(os.path.join(d, '_COMPLETE'), 'w').write('ok\n')
    os.rename(d, final)
    return final
def old_records(d, rid):
    fn = os.path.join(d, rid + '.jsonl')
    if not os.path.exists(fn): return []
    return [json.loads(l) for l in open(fn)]

def recorded_n(out_dir):
    """the worker count recorded by recert_all.sh at the first start (out_dir/RECERT_N)"""
    try:
        with open(os.path.join(out_dir, 'RECERT_N')) as fh: s = fh.read().strip()
    except OSError: raise SystemExit(f'{out_dir}/RECERT_N missing: start the re-certification with recert_all.sh')
    if not re.fullmatch(r'[1-9][0-9]?', s): raise SystemExit(f'{out_dir}/RECERT_N malformed: {s!r}')
    return int(s)

def claim(cdir, rid, wid):
    """atomic claim: the claim file appears (link(2)) already holding the worker ID; False if the root is claimed"""
    tmp = os.path.join(cdir, f'.tmp.{wid}.{os.getpid()}')
    with open(tmp, 'w') as fh: fh.write(str(wid)); fh.flush(); os.fsync(fh.fileno())
    try: os.link(tmp, os.path.join(cdir, rid)); return True
    except FileExistsError: return False
    finally: os.unlink(tmp)

def cut_torn_tail(path):
    """a kill during a write leaves a torn last line: cut it off (the record was never complete; its work is redone)"""
    if not os.path.exists(path) or os.path.getsize(path) == 0: return 0
    with open(path, 'rb+') as fh:
        data = fh.read()
        if data.endswith(b'\n'): return 0
        cut = data.rfind(b'\n') + 1
        fh.truncate(cut); fh.flush(); os.fsync(fh.fileno())
    print(json.dumps(dict(event='torn_tail_removed', file=os.path.basename(path), bytes=len(data) - cut)), flush=True)   # RIGOR: info: progress
    return len(data) - cut

def code_id():
    """sha256 over the re-certification code: recert.py's local import closure, recert_all.sh and certs_complete.py"""
    from binding import closure, sha_file
    files = sorted(set(closure(['recert.py', 'certs_complete.py'])) | {'recert_all.sh'})
    return hashlib.sha256(''.join(f'{f} {sha_file(os.path.join(HERE, f))}\n' for f in files).encode()).hexdigest()

def worker(name, wid):
    from cover5 import Cover5
    import fpenv
    fpenv.check()                                            # certificates are computed under IEEE round-to-nearest only
    out_dir = os.environ.get('RECERT_DIR', os.path.join(HERE, 'certs4')); os.makedirs(out_dir, exist_ok=True)
    N = recorded_n(out_dir)
    if not (re.fullmatch(r'[1-9][0-9]?', str(wid)) and int(wid) <= N): raise SystemExit(f'worker ID {wid!r} is not in 1..{N} (RECERT_N)')
    wid = str(int(wid))
    # one process per worker ID: a second one would interleave the same output file (advisory lock, released at exit)
    os.makedirs(os.path.join(out_dir, 'locks'), exist_ok=True)
    lk = os.open(os.path.join(out_dir, 'locks', f'w{wid}.lock'), os.O_CREAT | os.O_RDWR, 0o644)
    try: fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError: raise SystemExit(f'worker {wid} is already running (locks/w{wid}.lock)')
    cdir = os.path.join(out_dir, 'claims', name); os.makedirs(cdir, exist_ok=True)
    out = os.path.join(out_dir, f'{name.lower()}.w{wid}.jsonl')
    man = json.load(open(os.path.join(HERE, 'manifests', f'{name.lower()}.json')))
    mroots = {rid: (p0, r) for rid, p0, r in man['roots']}
    # every claim must name a manifest root and a worker 1..N: a claim of any other owner would never be finished
    claims = {}; bad = []
    for c in os.listdir(cdir):
        if c.startswith('.'): continue
        with open(os.path.join(cdir, c)) as fh: owner = fh.read().strip()
        if c not in mroots or not re.fullmatch(r'[1-9][0-9]?', owner) or int(owner) > N: bad.append(f'{c}: owner {owner!r}')
        claims[c] = owner
    if bad: raise SystemExit(f'{len(bad)} invalid claims in {cdir} (N = {N}), e.g. {bad[:5]}: refusing to start')
    odir = os.path.join(out_dir, 'old', name)
    if not os.path.exists(os.path.join(odir, '_COMPLETE')): raise SystemExit(f'run recert.py --split {name} first')
    print(json.dumps(dict(event='start', family=name, worker=int(wid), N=N, pid=os.getpid(), niceness=os.nice(0),   # RIGOR: info: progress
                          isolated=sys.flags.isolated, no_site=sys.flags.no_site, dont_write_bytecode=sys.flags.dont_write_bytecode)), flush=True)
    cut_torn_tail(out)
    cx = Ctx(name)
    rng = random.Random(20261002 + int(wid)); order = sorted(mroots); rng.shuffle(order)
    cv = Cover5(out, cx.fam, fem_problem=run_fam.fem_problem(name))
    # what is already in the output (restart)
    have = set(); leaves_new = {}; items = set(); rdone = set()     # leaves_new: (root, path) -> status only (memory: ~0.1 kB, not ~3 kB, per leaf)
    with open(out) as f:
        for line in f:
            try: d = json.loads(line)
            except Exception: continue
            t = d.get('type')
            if t in ('split', 'esc', 'leaf'): have.add((t, d['root'], d['path']))
            if t == 'leaf': leaves_new[(d['root'], d['path'])] = d['status']
            if t == 'item_done': items.add((d['root'], d['path']))
            if t == 'root_done': rdone.add(d['root'])
    mine = {x for x, o in claims.items() if o == wid}
    queue = [x for x in order if x in mine and x not in rdone] + [x for x in order if x not in mine]
    t_all = time.time(); nroots = 0; stats = collections.Counter()
    for rid in queue:
        if rid not in mine and not claim(cdir, rid, wid): continue
        if rid in rdone: continue
        t0 = time.time(); counts = collections.Counter()
        recs = old_records(odir, rid)
        for d in recs:                                     # structural records first (exact tiling)
            if d['type'] in ('split', 'esc') and (d['type'], rid, d['path']) not in have:
                cv.emit(d); have.add((d['type'], rid, d['path']))
        for d in recs:
            if d['type'] != 'leaf': continue
            key = (rid, d['path'])
            if key in leaves_new:
                st = leaves_new[key]; counts[st if st in ('cert', 'FAIL') else 'discard'] += 1; continue   # (a FAIL leaf keeps root_done ok = false)
            if key in items: continue
            fpenv.check()
            try: new, info = recert_leaf(cx, d)
            except Exception as e: new, info = None, 'exception: ' + repr(e)[:200]
            stats[info.split(':')[0]] += 1
            if new is not None:
                cv.emit(new); leaves_new[key] = new['status']; counts['cert' if new['status'] == 'cert' else 'discard'] += 1
            else:
                print(json.dumps(dict(event='subdivide', root=rid, path=d['path'], old=d.get('method', d['status']), why=info)), flush=True)
                depth = len(d['path']) - len(rid)
                k = d.get('k', cv.kmin) if d['status'] == 'cert' else cv.kmin
                res = cv.process_item(rid, d['path'], d['p0'], d['r'], depth, min(int(k), cv.kmax), None)
                items.add(key)
                for kk, v in (res.get('counts') or {}).items(): counts[kk] += v
                stats['subdivided'] += 1
        if rid not in cv.done:
            sp0, sr = dyadic_root(*mroots[rid])
            cv.emit(dict(type='root_done', root=rid, ok=counts['FAIL'] == 0, time=time.time() - t0, counts=dict(counts),
                         box=dict(p0=[float(x) for x in sp0], r=[float(x) for x in sr])))   # RIGOR: exact: float() of exact floats
        nroots += 1
        print(json.dumps(dict(root=rid, dt=round(time.time() - t0, 2), counts=dict(counts), wall=round(time.time() - t_all, 1))), flush=True)   # RIGOR: info: progress
    cv.close()
    print('DONE', json.dumps(dict(roots=nroots, **stats)), round(time.time() - t_all, 1), flush=True)   # RIGOR: info: progress

def time_sample(name, n, seed=1):
    cx = Ctx(name)
    by = collections.defaultdict(list)
    with open(os.path.join(HERE, 'certs3', f'{name.lower()}.merged.jsonl')) as f:
        for line in f:
            if '"leaf"' not in line: continue
            d = json.loads(line)
            if d['type'] != 'leaf': continue
            kind = d['status'] if d['status'] == 'discard' else d['method']
            if kind == 'FEM' and name in ('PP', 'PF', 'FP') and d['pc'][1] == 0.0: kind = 'FEM-FARL'
            by[kind].append(d)
    rng = random.Random(seed); out = {}
    for kind, lst in sorted(by.items()):
        pick = rng.sample(lst, min(n, len(lst))); ts = []; fails = []
        for d in pick:
            t0 = time.time()
            try: new, info = recert_leaf(cx, d)
            except Exception as e: new, info = None, 'exception: ' + repr(e)[:200]
            ts.append(time.time() - t0)
            if new is None: fails.append((d['path'], info))
        out[kind] = dict(n=len(pick), total=len(lst), mean_s=float(np.mean(ts)), max_s=float(np.max(ts)), fails=len(fails))
        print(json.dumps(dict(family=name, kind=kind, n=len(pick), total=len(lst), mean_s=round(float(np.mean(ts)), 4),   # RIGOR: info: timing
                              p90_s=round(float(np.percentile(ts, 90)), 4), max_s=round(float(np.max(ts)), 4), fails=len(fails),   # RIGOR: info: timing
                              fail_examples=fails[:3])), flush=True)
    return out

if __name__ == '__main__':
    if sys.argv[1] == '--split':
        print(split_old(sys.argv[2]))
    elif sys.argv[1] == '--time':
        time_sample(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]) if len(sys.argv) > 4 else 1)
    elif sys.argv[1] == '--code-id':
        print(code_id())
    else:
        worker(sys.argv[1], sys.argv[2])
