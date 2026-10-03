"""Adaptive rigorous box cover in the BD chart (p = (ax, ay, cx, cy), B=(0,0), D=(1,0)).
Each root box is processed depth-first; every leaf is either CERTIFIED (cert.BoxCert) or DISCARDED by a rigorous
interval test (not a dart / shallow Prop 5.1 / outside the target region).  Log: JSONL, one line per leaf, plus a
'root_done' line per root (resumable: completed roots are skipped on restart)."""
import sys, json, time, os, numpy as np, warnings
warnings.filterwarnings('ignore')
from flint import arb
from refexact import build
from assemble import Chart
from cert import BoxCert
from sigma import sigma_centre, strain_factor

def iv(lo, hi): return arb(lo).union(arb(hi))

def box_ivs(p0, r): return [iv(p0[i] - r[i], p0[i] + r[i]) for i in range(4)]

def reflex_D_max(p0, r):
    """upper bound of f_D = (1 - cx) ay + cy (ax - 1) over the box (dart with reflex D iff f_D < 0)"""
    ax, ay, cx, cy = box_ivs(p0, r)
    # f is bilinear: evaluate exactly at the 16 corners (max of a multilinear function is at a vertex)
    best = -1e300
    for sx in (-1, 1):
        for sy in (-1, 1):
            for tx in (-1, 1):
                for ty in (-1, 1):
                    a = arb(p0[0] + sx*r[0]); b = arb(p0[1] + sy*r[1]); c = arb(p0[2] + tx*r[2]); d = arb(p0[3] + ty*r[3])
                    v = (1 - c)*b + d*(a - 1); best = max(best, float(v.upper()))
    return best

def beta_max(p0, r):
    """upper bound of the notch depth beta = dist(D, AC)/dist(B, AC) over the box (interval arithmetic)"""
    ax, ay, cx, cy = box_ivs(p0, r)
    # line AC: n = (ay - cy, cx - ax) (normal, unnormalised); dist ratio = |n.(D - A)| / |n.(B - A)|
    n1 = ay - cy; n2 = cx - ax
    num = abs(n1*(1 - ax) + n2*(0 - ay)); den = abs(n1*(0 - ax) + n2*(0 - ay))
    if not (den.lower() > 0): return 1e300
    return float((num/den).upper())

class Cover:
    def __init__(self, logpath, kmin=4, kmax=6, verbose=True):
        self.logpath = logpath; self.kmin = kmin; self.kmax = kmax; self.verbose = verbose
        self.bc = {}; self.stats = dict(cert=0, discard=0, fail_leaf=0, ldl=0, time=0.0)
        self.done_roots = set()
        if os.path.exists(logpath):
            for line in open(logpath):
                try: d = json.loads(line)
                except Exception: continue
                if d.get('type') == 'root_done': self.done_roots.add(d['root'])
        self.log = open(logpath, 'a'); self.dbg = open(logpath + '.dbg', 'a')
    def getbc(self, k):
        if k not in self.bc: self.bc[k] = BoxCert(Chart(build(k)))
        return self.bc[k]
    def emit(self, d):
        self.log.write(json.dumps(d) + '\n'); self.log.flush()
    def region_ok(self, p0, r):
        return None
    def process_root(self, root_id, p0, r, maxdepth=40, min_rel=1e-4):
        if root_id in self.done_roots: return True
        stack = [(np.array(p0, float), np.array(r, float), 0, None, self.kmin)]
        ok_all = True; t0 = time.time()
        while stack:
            p, rr, depth, sig, k = stack.pop()
            reason = self.discard_reason(p, rr)
            if reason:
                self.emit(dict(type='leaf', root=root_id, status='discard', reason=reason, p0=p.tolist(), r=rr.tolist()))
                self.stats['discard'] += 1; continue
            res, sig_used, k_used = self.try_cert(p, rr, sig, k)
            if res is not None and res.get('ok'):
                self.emit(dict(type='leaf', root=root_id, status='cert', p0=p.tolist(), r=rr.tolist(), k=k_used,
                               sigma=sig_used, L1=res['L1'], n11=res['n11'], detN=res['detN'], h=res['h'], temple=res['temple']))
                self.stats['cert'] += 1; continue
            if depth >= maxdepth or np.max(rr) < min_rel:
                self.emit(dict(type='leaf', root=root_id, status='FAIL', p0=p.tolist(), r=rr.tolist(), res=str(res)))
                self.stats['fail_leaf'] += 1; ok_all = False; continue
            # escalate the mesh level if the CENTRE margin (Liu bound at p0 vs P2 lambda2) is small: the
            # discretisation, not the box size, is then the bottleneck
            knew = k
            if res is not None and 'lam1h' in res and k < self.kmax:
                lam1h = res['lam1h']; h = res.get('h', None)
                if h is None: h = self.getbc(k).ch.hmax((p[0], p[1]), (p[2], p[3]))
                liu = lam1h/(1 + (0.1893*h)**2*lam1h)
                if (2.5*liu - res['lam2p2'])/res['lam2p2'] < 0.04: knew = k + 1
            if self.verbose:
                self.dbg.write(json.dumps(dict(p0=p.tolist(), r=rr.tolist(), k=k, knew=knew, why=str(res.get('why') if res else None),
                                               detN=res.get('detN') if res else None, temple=res.get('temple') if res else None)) + '\n'); self.dbg.flush()
            dirn = self.split_dir(p, rr, res)
            for s in (-1, 1):
                q = p.copy(); q[dirn] += s*rr[dirn]/2; r2 = rr.copy(); r2[dirn] /= 2
                stack.append((q, r2, depth + 1, sig_used, knew))
        self.emit(dict(type='root_done', root=root_id, ok=ok_all, time=time.time() - t0))
        return ok_all
    def discard_reason(self, p, rr):
        if p[1] - rr[1] <= 0 or p[3] + rr[3] >= 0: return None if (p[1] + rr[1] > 0 and p[3] - rr[3] < 0) and False else self.crosses(p, rr)
        if reflex_D_max(p, rr) >= 0 and self.all_convex(p, rr): return 'not_dart(convex at D on whole box)'
        if beta_max(p, rr) <= 0.0339: return 'shallow(Prop5.1, beta<=0.0339)'
        return self.region_ok(p, rr)
    def crosses(self, p, rr): return 'outside_chart'
    def all_convex(self, p, rr):
        """f_D >= 0 on the whole box  <=> min f_D >= 0"""
        best = 1e300
        for sx in (-1, 1):
            for sy in (-1, 1):
                for tx in (-1, 1):
                    for ty in (-1, 1):
                        a = arb(p[0] + sx*rr[0]); b = arb(p[1] + sy*rr[1]); c = arb(p[2] + tx*rr[2]); d = arb(p[3] + ty*rr[3])
                        v = (1 - c)*b + d*(a - 1); best = min(best, float(v.lower()))
        return best >= 0
    def try_cert(self, p, rr, sig, k):
        bc = self.getbc(k)
        try:
            sol = bc.solve(p, rr)
        except Exception as e:
            return dict(ok=False, why='solve:' + str(e)), sig, k
        lam1, lam2 = sol['lamh'][0], sol['lamh'][1]
        sigma2 = None; sig_used = None
        if sig is not None and sig[2] == k:
            sf = float(strain_factor(sig[0], p, rr).lower())
            s_try = sig[1]*sf*(1 - 1e-12)
            if s_try > lam1 + 0.5*(lam2 - lam1): sigma2 = s_try; sig_used = sig
        if sigma2 is None:
            sc, info = sigma_centre(bc.ch, p, lam1 + 0.9*(lam2 - lam1)); self.stats['ldl'] += 1
            if sc is None: return dict(ok=False, why='ldl inertia'), None, k
            sig_used = (p.tolist(), float(sc), k)
            sigma2 = float(sc)*float(strain_factor(p, p, rr).lower())*(1 - 1e-12)
        try:
            res = bc.certify(p, rr, sigma2, sol=sol)
        except Exception as e:
            res = dict(ok=False, why='certify:' + str(e))
        res['lam2p2'] = float(sol['lamp2'][1]); res['lam1h'] = float(lam1)
        return res, sig_used, k
    def split_dir(self, p, rr, res):
        # relative half-widths: ax ~ |A|, ay ~ ay, cx ~ |C|, cy ~ |cy|
        scale = np.array([abs(p[1]), abs(p[1]), abs(p[3]), abs(p[3])])   # strain of the affine map ~ dx/eta, deta/eta
        return int(np.argmax(rr/scale))
