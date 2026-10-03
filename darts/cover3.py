"""Adaptive rigorous box cover, crash-safe and resumable at LEAF level (successor of cover2.Cover2; same certificate).

Every box has a path id: root id, then one character per step ('0'/'1' = lower/upper child of a bisection, 'k' = the same box
retried at a higher mesh level).  The JSONL log holds one record per decision:
  {"type": "leaf",  "path", "root", "status": cert|discard|defer|FAIL, "p0", "r", ...}
  {"type": "split", "path", "root", "dirn", "knew"}        (box bisected in coordinate dirn; children start at mesh level knew)
  {"type": "esc",   "path", "root", "knew"}                (box retried at mesh level knew)
  {"type": "root_done", "root", "ok", "box", "time", "counts"}
On restart the log is read back: finished leaves are skipped and recorded splits/escalations are replayed without recomputation,
so a crash costs at most the box being processed.  A torn last line (crash during a write) is ignored and a newline is written
before appending.  The log is fsync'ed every FSYNC_S seconds.

Statuses:  cert    -- rigorous certificate (gcert.GCert.certify with sigma2 from gsigma.sigma_box)
           discard -- rigorous region predicate (discard_fn returns the reason; replay re-checks it)
           defer   -- handed to a downstream method (defer_fn returns the reason); the box must be covered by another log
                      (checked by accounting.py) -- it is NOT proved here
           FAIL    -- refined below the floor (maxdepth) without success: a LOUD failure (also printed to stderr)."""
# RIGOR-MODULE: heuristic: legacy production search driver; never executed by a replay or the accounting
import json, time, os, sys, numpy as np, warnings
warnings.filterwarnings('ignore')
from gcert import GCert
from gsigma import sigma_box
from grid import dyadic_root

FSYNC_S = 30.0

def read_log(logpath):
    """(root, path) -> ('leaf', status) | ('split', dirn, knew) | ('esc', knew);  root -> root_done record"""
    prev = {}; done = {}
    if not os.path.exists(logpath): return prev, done
    with open(logpath) as f:
        for line in f:
            try: d = json.loads(line)
            except Exception: continue                      # torn line
            t = d.get('type')
            # keys are (root, path): root ids may be prefixes of each other (e.g. X_8 / X_81), so path strings alone collide
            if t == 'leaf' and 'path' in d: prev[(d['root'], d['path'])] = ('leaf', d['status'])
            elif t == 'split': prev[(d['root'], d['path'])] = ('split', int(d['dirn']), int(d['knew']))
            elif t == 'esc': prev[(d['root'], d['path'])] = ('esc', int(d['knew']))
            elif t == 'root_done': done[d['root']] = d
    return prev, done

class Cover3:
    def __init__(self, logpath, problem_fn, discard_fn=None, defer_fn=None, kmin=4, kmax=6, split_scale_fn=None,
                 margin_escalate=0.04, maxdepth=30, verbose=True):
        self.logpath = logpath; self.problem_fn = problem_fn; self.discard_fn = discard_fn; self.defer_fn = defer_fn
        self.kmin = kmin; self.kmax = kmax; self.split_scale_fn = split_scale_fn; self.margin_escalate = margin_escalate
        self.maxdepth = maxdepth; self.gc = {}
        self.stats = dict(cert=0, discard=0, defer=0, FAIL=0, ldl=0, resumed=0)
        self.prev, self.done = read_log(logpath)
        if os.path.exists(logpath) and os.path.getsize(logpath) > 0:
            with open(logpath, 'rb') as f:
                f.seek(-1, 2); last = f.read(1)
            if last != b'\n':
                with open(logpath, 'a') as f: f.write('\n')
        self.log = open(logpath, 'a'); self.t_sync = time.time()
        self.dbg = open(logpath + '.dbg', 'a') if verbose else None
    def getgc(self, k):
        if k not in self.gc: self.gc[k] = GCert(self.problem_fn(k))
        return self.gc[k]
    def emit(self, d):
        self.log.write(json.dumps(d) + '\n'); self.log.flush()
        if time.time() - self.t_sync > FSYNC_S:
            os.fsync(self.log.fileno()); self.t_sync = time.time()
    def close(self):
        self.log.flush(); os.fsync(self.log.fileno()); self.log.close()
        if self.dbg: self.dbg.close()
    def process_root(self, root_id, p0, r):
        if root_id in self.done: return bool(self.done[root_id].get('ok'))
        p0, r = dyadic_root(p0, r)
        stack = [(root_id, np.array(p0, float), np.array(r, float), 0, None, self.kmin)]
        counts = dict(cert=0, discard=0, defer=0, FAIL=0); t0 = time.time()
        while stack:
            path, p, rr, depth, sig, k = stack.pop()
            pr = self.prev.get((root_id, path))
            if pr is not None:
                self.stats['resumed'] += 1
                if pr[0] == 'leaf':
                    if pr[1] in counts: counts[pr[1]] += 1
                    continue
                if pr[0] == 'split':
                    dirn, knew = pr[1], pr[2]
                    for s, tag in ((-1, '0'), (1, '1')):
                        q = p.copy(); q[dirn] += s*rr[dirn]/2; r2 = rr.copy(); r2[dirn] /= 2
                        stack.append((path + tag, q, r2, depth + 1, None, knew))
                    continue
                if pr[0] == 'esc':
                    stack.append((path + 'k', p, rr, depth + 1, None, pr[1])); continue
            reason = self.discard_fn(p, rr) if self.discard_fn else None
            if reason:
                self.emit(dict(type='leaf', path=path, root=root_id, status='discard', reason=reason, p0=p.tolist(), r=rr.tolist()))
                self.stats['discard'] += 1; counts['discard'] += 1; continue
            res, sig_used = self.try_cert(p, rr, sig, k)
            if res.get('ok'):
                self.emit(self.cert_record(path, root_id, p, rr, k, sig_used, res))
                self.stats['cert'] += 1; counts['cert'] += 1; continue
            dreason = self.defer_fn(p, rr, res, depth) if self.defer_fn else None
            if dreason:
                self.emit(dict(type='leaf', path=path, root=root_id, status='defer', reason=dreason, p0=p.tolist(), r=rr.tolist(),
                               xi_h=res.get('xi_h')))
                self.stats['defer'] += 1; counts['defer'] += 1; continue
            if depth >= self.maxdepth:
                self.emit(dict(type='leaf', path=path, root=root_id, status='FAIL', p0=p.tolist(), r=rr.tolist(), res=str(res)[:400]))
                self.stats['FAIL'] += 1; counts['FAIL'] += 1
                print(f'LOUD FAILURE: box refined to the floor without a certificate: root {root_id} path {path} p0 {p.tolist()} '
                      f'r {rr.tolist()} why {res.get("why")}', file=sys.stderr, flush=True)
                continue
            knew = k
            if 'lam1h' in res and k < self.kmax:
                h = res.get('h') or 0.0
                liu = res['lam1h']/(1 + (0.1893*h)**2*res['lam1h']) if h else res['lam1h']
                if (2.5*liu - res['lam2p2'])/res['lam2p2'] < self.margin_escalate: knew = k + 1
            if self.dbg:
                self.dbg.write(json.dumps(dict(path=path, p0=p.tolist(), r=rr.tolist(), k=k, knew=knew, why=res.get('why'),
                                               schur=res.get('schur'), temple=res.get('temple'), xi_h=res.get('xi_h'))) + '\n')
                self.dbg.flush()
            if knew != k and depth > 0:
                self.emit(dict(type='esc', path=path, root=root_id, knew=knew))
                stack.append((path + 'k', p, rr, depth + 1, None, knew)); continue
            dirn = self.split_dir(p, rr)
            self.emit(dict(type='split', path=path, root=root_id, dirn=int(dirn), knew=int(knew)))
            for s, tag in ((-1, '0'), (1, '1')):
                q = p.copy(); q[dirn] += s*rr[dirn]/2; r2 = rr.copy(); r2[dirn] /= 2
                stack.append((path + tag, q, r2, depth + 1, sig_used if knew == k else None, knew))
        ok = counts['FAIL'] == 0
        rec = dict(type='root_done', root=root_id, ok=ok, time=time.time() - t0, counts=counts,
                   box=dict(p0=[float(x) for x in p0], r=[float(x) for x in r]))
        self.emit(rec); self.done[root_id] = rec
        return ok
    def cert_record(self, path, root_id, p, rr, k, sig_used, res):
        if res.get('method') == 'BO':
            return dict(type='leaf', path=path, root=root_id, status='cert', method='BO', p0=p.tolist(), r=rr.tolist(), k=k,
                        E=res['E'], theta=res.get('theta'), K=res.get('K'), schur=res['schur'], lam2p2=res.get('lam2p2'))
        return dict(type='leaf', path=path, root=root_id, status='cert', method='FEM', p0=p.tolist(), r=rr.tolist(), k=k,
                    sigma=sig_used, L1=res['L1'], schur=res['schur'], h=res['h'], temple=res['temple'], tails=res.get('tails'))
    def try_cert(self, p, rr, sig, k):
        gc = self.getgc(k)
        try: sol = gc.solve(p, rr)
        except Exception as e: return dict(ok=False, why='solve:' + repr(e)[:200]), sig
        lam1 = float(sol['lamh'][0]); lam2p2 = float(sol['lamp2'][1])
        base = dict(lam1h=lam1, lam2p2=lam2p2, xi_h=lam2p2/lam1, lam2h=float(sol['lamh'][1]) if len(sol['lamh']) > 1 else None)
        sc, info = sigma_box(gc, sol, p, rr); self.stats['ldl'] += 1
        if sc is None:
            d = dict(ok=False, why='sigma_box:' + str(info.get('why'))); d.update(base); return d, None
        try: res = gc.certify(p, rr, sc, sol=sol)
        except Exception as e: res = dict(ok=False, why='certify:' + repr(e)[:200])
        res.update(base)
        return res, (p.tolist(), float(sc), k, info.get('frac'))
    def split_dir(self, p, rr):
        scale = self.split_scale_fn(p) if self.split_scale_fn else np.where(np.abs(p) > 1e-9, np.abs(p), 1.0)
        rel = np.where(rr > 0, rr/scale, -1.0)
        return int(np.argmax(rel))
