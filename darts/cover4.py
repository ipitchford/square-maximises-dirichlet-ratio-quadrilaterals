"""Hybrid adaptive rigorous cover: every box is tried first with the BO certificate (bocert.py: 1-D lower bound for lambda_1,
P2 upper bound for lambda_2) and, if that fails, with the FEM certificate (gcert: Crouzeix-Raviart/Liu + Temple for lambda_1,
P2 for lambda_2); boxes that fail both are bisected (or retried at a finer mesh).  Crash-safe and resumable at leaf level
(cover3.Cover3: path ids, split/esc records, fsync).  Leaf records carry 'method' = 'BO' (with the slicing direction 'theta',
'K') or 'FEM' (as in cover3), so replay.py can re-run exactly the logged certificate."""
import json, time, numpy as np
from cover3 import Cover3
from bocert import bo_cert_bd, bo_cert_spike
from p2up import solve_p2

class Cover4(Cover3):
    def __init__(self, logpath, problem_fn, bo_kind=None, bo_K=120, fem=True, **kw):
        super().__init__(logpath, problem_fn, **kw)
        self.bo_kind = bo_kind; self.bo_K = bo_K; self.fem = fem
        self.stats.update(bo=0, fem=0)
    def try_cert(self, p, rr, sig, k):
        gc = self.getgc(k)
        if self.bo_kind is not None:
            try:
                sol2 = solve_p2(gc, p, rr)
                if self.bo_kind == 'bd': res = bo_cert_bd(gc, p, rr, self.bo_K, sol=sol2)
                else: res = bo_cert_spike(gc, p, rr, self.bo_K, sol=sol2)
            except Exception as e:
                res = dict(ok=False, why='bo:' + repr(e)[:200])
            if res.get('ok'):
                self.stats['bo'] += 1
                return res, None
            bo_why = res.get('why'); xi_bo = res.get('xi_bo')
        else:
            bo_why = None; xi_bo = None
        if not self.fem:
            return dict(ok=False, why=f'bo: {bo_why}', xi_bo=xi_bo), None
        res, sig_used = super().try_cert(p, rr, sig, k)
        res['bo_why'] = bo_why; res['xi_bo'] = xi_bo
        if res.get('ok'):
            res['method'] = 'FEM'; self.stats['fem'] += 1
        return res, sig_used
    def emit(self, d):
        # leaf records of certified boxes: keep the method-specific data needed by replay
        super().emit(d)

    def process_root(self, root_id, p0, r):
        return super().process_root(root_id, p0, r)
