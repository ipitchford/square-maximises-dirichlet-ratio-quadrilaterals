"""Adaptive rigorous box cover for a gcert.Problem family (one problem per mesh level k).
Leaves are CERTIFIED (gcert.GCert.certify with sigma2 from gsigma) or DISCARDED by a rigorous predicate (discard_fn returns a
reason string or None).  Unused parameter directions have r_i = 0 and are never split.
Log: JSONL (one line per leaf + one 'root_done' line per root); roots already completed are skipped (resumable)."""
# RIGOR-MODULE: heuristic: production search driver (Cover2); never executed by a replay or the accounting (they import dyadic_root from grid.py)
import json, time, os, numpy as np, warnings
warnings.filterwarnings('ignore')
from gcert import GCert
from gsigma import sigma_box

from grid import dyadic_root

class Cover2:
    def __init__(self, logpath, problem_fn, discard_fn=None, kmin=4, kmax=6, split_scale_fn=None, margin_escalate=0.04,
                 maxdepth=40, verbose=True):
        self.logpath = logpath; self.problem_fn = problem_fn; self.discard_fn = discard_fn
        self.kmin = kmin; self.kmax = kmax; self.split_scale_fn = split_scale_fn; self.margin_escalate = margin_escalate
        self.maxdepth = maxdepth; self.gc = {}
        self.stats = dict(cert=0, discard=0, fail_leaf=0, ldl=0)
        self.done_roots = set()
        if os.path.exists(logpath):
            for line in open(logpath):
                try: d = json.loads(line)
                except Exception: continue
                if d.get('type') == 'root_done': self.done_roots.add(d['root'])
        self.log = open(logpath, 'a'); self.dbg = open(logpath + '.dbg', 'a') if verbose else None
    def getgc(self, k):
        if k not in self.gc: self.gc[k] = GCert(self.problem_fn(k))
        return self.gc[k]
    def emit(self, d): self.log.write(json.dumps(d) + '\n'); self.log.flush()
    def process_root(self, root_id, p0, r):
        if root_id in self.done_roots: return True
        p0, r = dyadic_root(p0, r)
        stack = [(np.array(p0, float), np.array(r, float), 0, None, self.kmin)]
        ok_all = True; t0 = time.time()
        while stack:
            p, rr, depth, sig, k = stack.pop()
            reason = self.discard_fn(p, rr) if self.discard_fn else None
            if reason:
                self.emit(dict(type='leaf', root=root_id, status='discard', reason=reason, p0=p.tolist(), r=rr.tolist()))
                self.stats['discard'] += 1; continue
            res, sig_used = self.try_cert(p, rr, sig, k)
            if res.get('ok'):
                self.emit(dict(type='leaf', root=root_id, status='cert', p0=p.tolist(), r=rr.tolist(), k=k, sigma=sig_used,
                               L1=res['L1'], schur=res['schur'], h=res['h'], temple=res['temple'], tails=res.get('tails')))
                self.stats['cert'] += 1; continue
            if depth >= self.maxdepth:
                self.emit(dict(type='leaf', root=root_id, status='FAIL', p0=p.tolist(), r=rr.tolist(), res=str(res)))
                self.stats['fail_leaf'] += 1; ok_all = False; continue
            knew = k
            if 'lam1h' in res and k < self.kmax:
                h = res.get('h') or 0.0
                liu = res['lam1h']/(1 + (0.1893*h)**2*res['lam1h']) if h else res['lam1h']
                if (2.5*liu - res['lam2p2'])/res['lam2p2'] < self.margin_escalate: knew = k + 1
            if self.dbg:
                self.dbg.write(json.dumps(dict(p0=p.tolist(), r=rr.tolist(), k=k, knew=knew, why=res.get('why'),
                                               schur=res.get('schur'), temple=res.get('temple'))) + '\n'); self.dbg.flush()
            if knew != k and depth > 0:
                stack.append((p, rr, depth + 1, None, knew)); continue
            dirn = self.split_dir(p, rr, res)
            for s in (-1, 1):
                q = p.copy(); q[dirn] += s*rr[dirn]/2; r2 = rr.copy(); r2[dirn] /= 2
                stack.append((q, r2, depth + 1, sig_used if knew == k else None, knew))
        self.emit(dict(type='root_done', root=root_id, ok=ok_all, time=time.time() - t0,
                       box=dict(p0=[float(x) for x in p0], r=[float(x) for x in r])))
        return ok_all
    def try_cert(self, p, rr, sig, k):
        gc = self.getgc(k)
        try: sol = gc.solve(p, rr)
        except Exception as e: return dict(ok=False, why='solve:' + repr(e)), sig
        lam1 = sol['lamh'][0]
        sc, info = sigma_box(gc, sol, p, rr); self.stats['ldl'] += 1
        if sc is None:
            return dict(ok=False, why='sigma_box:' + str(info.get('why')), lam1h=float(lam1), lam2p2=float(sol['lamp2'][1])), None
        try: res = gc.certify(p, rr, sc, sol=sol)
        except Exception as e: res = dict(ok=False, why='certify:' + repr(e))
        res['lam2p2'] = float(sol['lamp2'][1]); res['lam1h'] = float(lam1)
        return res, (p.tolist(), float(sc), k, info.get('frac'))

    def split_dir(self, p, rr, res):
        scale = self.split_scale_fn(p) if self.split_scale_fn else np.where(np.abs(p) > 1e-9, np.abs(p), 1.0)
        rel = np.where(rr > 0, rr/scale, -1.0)
        return int(np.argmax(rel))
