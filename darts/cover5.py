"""Adaptive rigorous cover over a parameter FAMILY (fam5.py), crash-safe and resumable at leaf level (same log protocol as
cover3: path ids, split/esc records, root_done; keys (root, path); torn-line repair; periodic fsync).

Per box:  1. BO2 certificate: E = max over slicing axes of the rigorous 1-D lower bound bo2.bo_family for the family's arms
             (lambda_1 >= E for every dart of the box); then for an inner dart Q_in (contained in every dart of the box)
             p2fix.p2_cert certifies lambda_2(Q_in) < (5/2) E; hence xi < 5/2 on the box.
          2. INC certificate (family PP, inc.py, LEMMAS Lemma 11.4): an outer dart Q+ containing and an inner dart Q- contained
             in every dart of the box; lambda_1 >= L1(Q+) (gcert at the point Q+), lambda_2 <= lambda_2(Q-) < (5/2) L1 (p2fix).
          3. FEM certificate (families with a Cartesian enclosure): gcert on the Cartesian BD box containing the family box.
          4. otherwise bisect (largest scaled half-width) -- a box at depth >= maxdepth is a LOUD failure.
Leaf records: status cert (method BO2 | FEM) or discard (reason = a predicate of the family, re-checked in replay)."""
# RIGOR-MODULE: heuristic: production search driver (methods, splitting, hints, work items); never executed by a replay or the accounting
import json, time, os, sys, math, numpy as np, warnings
warnings.filterwarnings('ignore')
from grid import dyadic_root
from cover3 import read_log, FSYNC_S
from bo2 import bo_family, arm_variant
from p2fix import p2_cert
CENTER_RULE_W = 0.1
INC_CLOSE = 2.9          # INC float ratio below which boxes are split for INC (contribution rule) and big-box FEM is skipped

class Cover5:
    def __init__(self, logpath, fam, fem_problem=None, K=160, kmin=4, kmax=6, maxdepth=30, verbose=True):
        self.logpath = logpath; self.fam = fam; self.fem_problem = fem_problem; self.K = K
        self.kmin = kmin; self.kmax = kmax; self.maxdepth = maxdepth; self.gc = {}
        self.stats = dict(cert=0, discard=0, FAIL=0, bo=0, fem=0, inc=0, ldl=0, resumed=0, t_bo=0.0, t_inc=0.0, t_fem=0.0, vol=0.0)
        self.fem_streak = 0; self.since_fem = 0; self.last_method = None      # FEM gating state (heuristic)
        self.prev, self.done = read_log(logpath)
        # work splitting (pool.py): sub-boxes handed to the queue, items finished (keys (root, path))
        self.handed = set(); self.items_done = {}; self.handoff_recs = []
        if os.path.exists(logpath):
            with open(logpath) as fh:
                for line in fh:
                    if '"handoff"' not in line and '"item_done"' not in line: continue
                    try: d = json.loads(line)
                    except Exception: continue
                    if d.get('type') == 'handoff': self.handed.add((d['root'], d['path'])); self.handoff_recs.append(d)
                    elif d.get('type') == 'item_done': self.items_done[(d['root'], d['path'])] = d
        if os.path.exists(logpath) and os.path.getsize(logpath) > 0:
            with open(logpath, 'rb') as f:
                f.seek(-1, 2); last = f.read(1)
            if last != b'\n':
                with open(logpath, 'a') as f: f.write('\n')
        self.log = open(logpath, 'a'); self.t_sync = time.time()
        self.dbg = open(logpath + '.dbg', 'a') if verbose else None
    def emit(self, d):
        if d.get('type') == 'leaf': self.stats['vol'] += float(np.prod([2*x for x in d['r'] if x > 0]))
        self.log.write(json.dumps(d) + '\n'); self.log.flush()
        if time.time() - self.t_sync > FSYNC_S:
            os.fsync(self.log.fileno()); self.t_sync = time.time()
    def close(self):
        self.log.flush(); os.fsync(self.log.fileno()); self.log.close()
        if self.dbg: self.dbg.close()
    def getgc(self, k):
        from gcert import GCert
        if k not in self.gc: self.gc[k] = GCert(self.fem_problem(k))
        return self.gc[k]
    # ------------------------------------------------------------------------------------------ certificates
    def try_bo2(self, p, rr):
        """BO2 certificate (Lemma 11.3).  Direction / arm-variant SELECTION uses the fast float screening bo2f.bo_family_f (a
        heuristic); the bound E of the record is always the rigorous bo2.bo_family at the recorded (theta, K, arms)."""
        from bo2f import bo_family_f
        fam = self.fam
        try:
            armA0, armC0 = fam.arms(p, rr)
            inners = fam.inners(p, rr)
        except Exception as e:
            return dict(ok=False, why='bo2 setup: ' + repr(e)[:150])
        if not inners: return dict(ok=False, why='no inner dart')
        if inners[0][0] is not None: inners = inners[:1]       # the C-only triangle lies in the (A, C) inner dart: never better
        long_ = any(a.kind == 'polar' and a.r2 >= 4.0 for a in (armA0, armC0) if a is not None)
        # float target: P2 lambda_2 of the inner dart (level 5 for long darts, where level 4 is far too pessimistic)
        try:
            from p2fix import _gc
            from p2up import solve_p2
            Ai, Ci = inners[0]
            kt = self.kmin + 1 if long_ else self.kmin
            if Ai is None: g0 = _gc('c', kt); q0 = np.array([1.0, 1.0, Ci[0], Ci[1]])
            else: g0 = _gc('bd', kt); q0 = np.array([Ai[0], Ai[1], Ci[0], Ci[1]])
            lam2f = float(solve_p2(g0, q0, np.zeros(4), order=0)['lamp2'][1])
        except Exception:
            lam2f = float('inf')
        # the family's directions plus a uniform grid of 32 (cheap in the float screening)
        dirs = list(fam.dirs(p, rr)) + [-math.pi/2 + i*math.pi/32 for i in range(32)]
        modes = ['polar']
        if armA0 is not None and armA0.kind == 'polar' and armA0.r2 >= 4.0: modes.append('farA')
        if armC0 is not None and armC0.kind == 'polar' and armC0.r2 >= 4.0: modes.append('farC')
        if len(modes) == 3: modes.append('farAC')
        # 0. fast path: the previous box's direction and arm variant (usually still good for a neighbouring box)
        best = (0.0, None, 'polar', self.K)
        lt, lm = getattr(self, 'last_theta', None), getattr(self, 'last_mode', 'polar')
        if lt is not None and lm in modes:
            aA, aC = arm_variant(armA0, armC0, lm)
            try: E0, _ = bo_family(aA, aC, lt, self.K)
            except Exception: E0 = 0.0
            best = (E0, lt, lm, self.K)
        # 1. float screening of all (mode, direction) pairs (skipped when the fast path already reaches the target)
        cand = []
        for mode in (modes if not best[0]*2.5 > lam2f*1.03 else []):
            aA, aC = arm_variant(armA0, armC0, mode)
            for th in dirs:
                try: ef = bo_family_f(aA, aC, th, self.K)
                except Exception: ef = 0.0
                cand.append((ef, th, mode))
        cand.sort(key=lambda c: -c[0])
        # 2. rigorous bound for the best few candidates
        for ef, th, mode in cand[:3]:
            if not ef > 0 or ef < 0.8*best[0]: continue
            aA, aC = arm_variant(armA0, armC0, mode)
            try: E, _ = bo_family(aA, aC, th, self.K)
            except Exception: E = 0.0
            if E > best[0]: best = (E, th, mode, self.K)
            if best[0]*2.5 > lam2f*1.03: break
        E, th, amode, Kused = best
        # 3. near miss: finer slab grid (2K) and a local direction search (screened in float, then rigorous)
        if E > 0 and np.isfinite(lam2f) and lam2f/1.4 < 2.5*E < 1.03*lam2f:
            aA, aC = arm_variant(armA0, armC0, amode)
            loc = []
            for dth in (0.0, -0.05, 0.05, -0.1, 0.1, -0.2, 0.2, -0.3, 0.3):
                try: loc.append((bo_family_f(aA, aC, th + dth, 2*self.K), th + dth))
                except Exception: pass
            loc.sort(key=lambda c: -c[0])
            for ef, th2 in loc[:2]:
                try: E2, _ = bo_family(aA, aC, th2, 2*self.K)
                except Exception: E2 = 0.0
                if E2 > E: E, th, Kused = E2, th2, 2*self.K
        self.last_theta = th; self.last_mode = amode
        if not E > 0: return dict(ok=False, why='E = 0')
        last = None; best_l2 = float('inf'); l2prev = None
        for ii, (Ai, Ci) in enumerate(inners):
            for k in (self.kmin, self.kmin + 1, self.kmin + 2):
                try: res = p2_cert(Ai, Ci, E, k)
                except Exception as e: res = dict(ok=False, why='p2: ' + repr(e)[:150], lam2p2=float('inf'))
                last = res; l2 = res.get('lam2p2', float('inf')); best_l2 = min(best_l2, l2)
                if res.get('ok'):
                    return dict(ok=True, method='BO2', E=E, theta=th, K=Kused, inner=ii, Ain=Ai, Cin=Ci, k=k,
                                schur=res['schur'], lam2p2=res['lam2p2'], arms=amode)
                # a finer P2 mesh lowers the bound by up to ~15% per level (long darts): stop when even that cannot suffice,
                # or when the bound has converged (< 2% change from the previous level)
                if l2 > 2.5*E/0.85**(self.kmin + 2 - k): break
                if l2prev is not None and l2 > 0.98*l2prev: break
                l2prev = l2
        return dict(ok=False, why='bo2: ' + str(last.get('why') if last else None), E=E,
                    xi_bo=(best_l2/E) if np.isfinite(best_l2) else None)
    def try_inc(self, p, rr, k0):
        from inc import inc_cert
        try: return inc_cert(p, rr, k0=k0, kmax=self.kmax)
        except Exception as e: return dict(ok=False, why='inc: ' + repr(e)[:150], ratio=float('inf'), ik=k0)
    def inc_contrib(self, p, rr, R):
        """split heuristic (no effect on rigour): estimated contribution of each half-width to the INC float ratio R"""
        from inc import inc_float
        c = np.zeros(4)
        for i in range(4):
            if not rr[i] > 0: c[i] = -1.0; continue
            r2 = rr.copy(); r2[i] /= 2
            try: Ri = inc_float(p, r2, self.kmin)[3]
            except Exception: Ri = float('inf')
            c[i] = 2*(R - Ri) if np.isfinite(Ri) else 0.0
        return np.maximum(c, np.where(rr > 0, 0.0, -1.0))
    def try_fem(self, p, rr, k):
        from gsigma import sigma_box
        if hasattr(self.fam, 'fem_setup'):
            try: st = self.fam.fem_setup(p, rr, k)
            except Exception as e: st = None
            if st is None: return dict(ok=False, why='no FEM chart (setup)'), None
            from gcert import GCert
            pb, pc, rc = st; pc, rc = np.array(pc), np.array(rc)
            key = ('farl', id(pb))
            if key not in self.gc: self.gc[key] = GCert(pb)
            gc = self.gc[key]
        else:
            cb = self.fam.cart_box(p, rr) if self.fem_problem is not None else None
            if cb is None: return dict(ok=False, why='no FEM chart'), None
            pc, rc = np.array(cb[0]), np.array(cb[1])
            gc = self.getgc(k)
        try: sol = gc.solve(pc, rc)
        except Exception as e: return dict(ok=False, why='fem solve: ' + repr(e)[:150]), None
        lam1 = float(sol['lamh'][0]); lam2p2 = float(sol['lamp2'][1])
        try: sc, info = sigma_box(gc, sol, pc, rc)
        except Exception as e: sc, info = None, dict(why='sigma: ' + repr(e)[:100])
        self.stats['ldl'] += 1
        if sc is None: return dict(ok=False, why='sigma_box:' + str(info.get('why')), lam1h=lam1, lam2p2=lam2p2), None
        try: res = gc.certify(pc, rc, sc, sol=sol)
        except Exception as e: res = dict(ok=False, why='certify: ' + repr(e)[:150])
        res.update(lam1h=lam1, lam2p2=lam2p2, pc=pc.tolist(), rc=rc.tolist())
        return res, (pc.tolist(), float(sc), k, info.get('frac'))
    # ------------------------------------------------------------------------------------------ cover
    def process_root(self, root_id, p0, r):
        """the whole root in one go (no work splitting): item + root_done"""
        if root_id in self.done: return bool(self.done[root_id].get('ok'))
        p0, r = dyadic_root(p0, r)
        res = self.process_item(root_id, root_id, p0, r, 0, self.kmin, None)
        return res['ok']
    def process_item(self, root_id, path0, p0, r0, depth0, k0, hint0, deadline=None, handoff_fn=None):
        """cover of the sub-box path0 of root root_id.  With a deadline: once it has passed, every pending sub-box is handed to the
        queue as an independent work item (log record 'handoff' first, then handoff_fn writes the queue entry).  Ends with an
        'item_done' record (and a 'root_done' record when the item is the whole root and nothing was handed off).
        Returns dict(ok, handed, counts)."""
        if path0 == root_id and root_id in self.done:
            d = self.done[root_id]; return dict(ok=bool(d.get('ok')), handed=0, counts=d.get('counts'))
        if (root_id, path0) in self.items_done:
            d = self.items_done[(root_id, path0)]; return dict(ok=bool(d.get('ok')), handed=d.get('handed_off', 0), counts=d.get('counts'))
        stack = [(path0, np.array(p0, float), np.array(r0, float), depth0, k0, hint0)]
        counts = dict(cert=0, discard=0, FAIL=0); t0 = time.time(); nhand = 0
        while stack:
            if deadline is not None and handoff_fn is not None and time.time() > deadline and len(stack) >= 2:
                for (hp, hq, hr, hd, hk, hh) in stack:
                    if hp != path0 and (root_id, hp) in self.handed: continue
                    self.emit(dict(type='handoff', root=root_id, path=hp, p0=hq.tolist(), r=hr.tolist(), depth=int(hd), k=int(hk)))
                    os.fsync(self.log.fileno())
                    handoff_fn(root_id, hp, hq.tolist(), hr.tolist(), int(hd), int(hk), hh)
                    self.handed.add((root_id, hp)); nhand += 1
                stack = []
                break
            path, p, rr, depth, k, hint = stack.pop()
            if path != path0 and (root_id, path) in self.handed: continue   # handed off from this item: another work item now
            hint = hint or {}
            pxi = hint.get('xi'); page = hint.get('age', 0); pcon = hint.get('c'); pik = hint.get('ik', self.kmin)
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
                        stack.append((path + tag, q, r2, depth + 1, knew, None))
                    continue
                if pr[0] == 'esc':
                    stack.append((path + 'k', p, rr, depth + 1, pr[1], None)); continue
            reason = self.fam.discard(p, rr)
            if reason:
                self.emit(dict(type='leaf', path=path, root=root_id, status='discard', reason=reason, p0=p.tolist(), r=rr.tolist()))
                self.stats['discard'] += 1; counts['discard'] += 1; continue
            # the box straddles a discard boundary if its CENTRE is discardable: bisect without certificate attempts
            # (only for boxes wider than CENTER_RULE_W in the scaled coordinates, so straddling boxes are eventually attempted)
            if depth < self.maxdepth and float(np.max(rr/self.fam.scale(p, 'bo', rr))) > CENTER_RULE_W and self.fam.discard(p, np.zeros(4)) is not None:
                sc = self.fam.scale(p, 'bo', rr); rel = np.where(rr > 0, rr/sc, -1.0); dirn = int(np.argmax(rel))
                self.emit(dict(type='split', path=path, root=root_id, dirn=dirn, knew=int(k)))
                for s, tag in ((-1, '0'), (1, '1')):
                    q = p.copy(); q[dirn] += s*rr[dirn]/2; r2 = rr.copy(); r2[dirn] /= 2
                    stack.append((path + tag, q, r2, depth + 1, k, hint))
                continue
            # efficiency heuristics (no effect on rigour):
            #  * BO2 is skipped only below a parent where BO2 was far (xi_bo > 3.5) AND INC was close (INC in charge), and is then
            #    re-attempted after two halvings (the 1-D bound improves quickly as boxes shrink);
            #  * FEM is gated by its recent record in this cover: tried while it keeps succeeding (and FIRST when the last
            #    certificate was a FEM one: in fat regions it certifies boxes 2-4x larger than BO2), skipped after 3 consecutive
            #    failures except for a periodic retry every 8 boxes.
            big = float(np.max(rr)) > 0.004
            fem_ok = self.fam.fem and self.fam.fem_plausible(p, rr)
            fres, sig, fem_tried = dict(ok=False, why='FEM skipped' if self.fam.fem else 'no FEM'), None, False
            def fem_attempt():
                t_ = time.time(); fr, sg = self.try_fem(p, rr, k); self.stats['t_fem'] += time.time() - t_
                self.fem_streak = 0 if fr.get('ok') else self.fem_streak + 1; self.since_fem = 0
                return fr, sg
            def emit_fem(fr, sg):
                self.emit(dict(type='leaf', path=path, root=root_id, status='cert', method='FEM', p0=p.tolist(), r=rr.tolist(),
                               k=k, pc=fr['pc'], rc=fr['rc'], sigma=sg, L1=fr['L1'], schur=fr['schur'], h=fr['h'],
                               temple=fr['temple'], tails=fr.get('tails')))
                self.stats['cert'] += 1; self.stats['fem'] += 1; counts['cert'] += 1; self.last_method = 'FEM'
            # method order: BO2, INC, FEM -- except that the method which certified the previous box goes first (regions are
            # spatially coherent: e.g. FEM certifies 2-4x larger boxes than BO2 in fat regions where the 1-D bound is ~20% low)
            cheap = getattr(self.fam, 'fem_cheap', False)       # BD-chart FEM (PP): ~0.1-0.2 s per attempt, never gated
            # PP: FEM before INC (the box FEM exploits the correlation of lambda_1 and lambda_2 over the box: in the bulk, xi ~ 2,
            # it certifies (0.016, 0.031, 0.031, 0.031) where INC needs (0.004-0.008, 0.031, 0.016, 0.031))
            order = ['BO2', 'FEM', 'INC'] if cheap else ['BO2', 'INC', 'FEM']
            if self.last_method in ('INC', 'FEM'): order.remove(self.last_method); order.insert(0, self.last_method)
            res = dict(ok=False, why='bo2 not tried'); nxi, nage = pxi, page + 1
            inc_ratio = None; ik = pik; bo_why = None; xi_bo = None; certified = False
            for meth in order:
                if meth == 'BO2':
                    if big and pxi is not None and pxi > 3.5 and hint.get('incc') and page < 2:
                        res = dict(ok=False, why='bo2 skipped (parent xi_bo > 3.5)', xi_bo=pxi); nxi, nage = pxi, page + 1
                    else:
                        t_ = time.time(); res = self.try_bo2(p, rr); nxi, nage = res.get('xi_bo'), 0; self.stats['t_bo'] += time.time() - t_
                    if res.get('ok'):
                        self.emit(dict(type='leaf', path=path, root=root_id, status='cert', method='BO2', p0=p.tolist(), r=rr.tolist(),
                                       k=res['k'], E=res['E'], theta=res['theta'], K=res['K'], inner=res['inner'],
                                       Ain=list(res['Ain']) if res['Ain'] is not None else None, Cin=list(res['Cin']),
                                       schur=res['schur'], lam2p2=res['lam2p2'], **({'arms': res['arms']} if res['arms'] != 'polar' else {})))
                        self.stats['cert'] += 1; self.stats['bo'] += 1; counts['cert'] += 1; self.last_method = 'BO2'
                        certified = True; break
                    bo_why = res.get('why'); xi_bo = res.get('xi_bo')
                elif meth == 'INC':
                    if not getattr(self.fam, 'inc', False): continue
                    t_ = time.time(); ires = self.try_inc(p, rr, pik); self.stats['t_inc'] += time.time() - t_
                    if ires.get('ok'):
                        self.emit(dict(type='leaf', path=path, root=root_id, status='cert', method='INC', p0=p.tolist(), r=rr.tolist(),
                                       k=ires['k'], kp=ires['kp'], Aout=ires['Aout'], Cout=ires['Cout'], Ain=ires['Ain'], Cin=ires['Cin'],
                                       E=ires['E'], sigma=ires['sigma'], schur=ires['schur'], lam2p2=ires['lam2p2'], lam1h=ires['lam1h'],
                                       h=ires['h'], temple=ires['temple'], ratio=ires['ratio']))
                        self.stats['cert'] += 1; self.stats['inc'] += 1; counts['cert'] += 1; self.last_method = 'INC'
                        certified = True; break
                    inc_ratio = ires.get('ratio'); ik = ires.get('ik', pik); bo_why = str(bo_why) + ' | ' + str(ires.get('why'))
                    if inc_ratio is not None and not np.isfinite(inc_ratio): inc_ratio = None
                else:
                    inc_close_ = inc_ratio is not None and inc_ratio < INC_CLOSE
                    skip_fem_ = (not cheap) and big and inc_close_ and inc_ratio >= 2.4   # INC in charge: FEM needs smaller boxes
                    first = (order[0] == 'FEM')
                    if fem_ok and not skip_fem_ and (first or cheap or self.fem_streak < 3 or self.since_fem >= 8):
                        fres, sig = fem_attempt(); fem_tried = True
                        if fres.get('ok'): emit_fem(fres, sig); certified = True; break
                    elif fem_ok:
                        self.since_fem += 1
            if certified: continue
            if xi_bo is None: xi_bo = res.get('xi_bo')
            inc_close = inc_ratio is not None and inc_ratio < INC_CLOSE                # INC within reach: split for INC
            if depth >= self.maxdepth:
                self.emit(dict(type='leaf', path=path, root=root_id, status='FAIL', p0=p.tolist(), r=rr.tolist(),
                               res=(str(bo_why) + ' | ' + str(fres.get('why')))[:400]))
                self.stats['FAIL'] += 1; counts['FAIL'] += 1
                print(f'LOUD FAILURE: {self.fam.name} root {root_id} path {path} p0 {p.tolist()} r {rr.tolist()} '
                      f'bo: {bo_why} fem: {fres.get("why")}', file=sys.stderr, flush=True)
                continue
            knew = k
            if 'lam1h' in fres and k < self.kmax and not str(fres.get('why', '')).startswith('sigma_box'):   # strain: no k help
                h = fres.get('h') or 0.0
                liu = fres['lam1h']/(1 + (0.1893*h)**2*fres['lam1h']) if h else fres['lam1h']
                if (2.5*liu - fres['lam2p2'])/fres['lam2p2'] < 0.04: knew = k + 1
            if self.dbg:
                self.dbg.write(json.dumps(dict(path=path, p0=p.tolist(), r=rr.tolist(), k=k, knew=knew, bo=bo_why, xi_bo=xi_bo,
                                               inc=inc_ratio, fem=fres.get('why'))) + '\n'); self.dbg.flush()
            if knew != k and depth > 0:
                self.emit(dict(type='esc', path=path, root=root_id, knew=knew))
                stack.append((path + 'k', p, rr, depth + 1, knew, dict(xi=nxi, age=nage, ik=ik))); continue
            # split scale: 'fem' (equal widths: the FEM strain) when the 1-D bound was far off, else 'bo' (relative angles)
            mode = 'fem' if (fem_ok and (xi_bo is None or xi_bo > 3.0)) else 'bo'
            ccon = None
            if inc_close and (xi_bo is None or xi_bo > 3.0):
                # INC split: halve the half-width with the largest estimated contribution to the INC ratio (estimated once on entry,
                # then inherited: a child's contribution vector is its parent's with the split entry halved)
                # (inherited; re-estimating at every split was tested and was slower: PP_1_6_2_8 18% vs 38% of the root in 10 min)
                ccon = np.array(pcon, float) if pcon is not None else self.inc_contrib(p, rr, inc_ratio)
                dirn = int(np.argmax(np.where(rr > 0, ccon, -2.0)))
                sc0 = self.fam.scale(p, 'fem', rr); rel0 = np.where(rr > 0, rr/sc0, np.nan)
                if not ccon[dirn] > 0 or np.nanmax(rel0)/(rr[dirn]/sc0[dirn]) > 32: ccon = None   # no signal / sliver: fall back
                else: ccon = ccon.copy(); ccon[dirn] /= 2
            if ccon is None:
                sc = self.fam.scale(p, mode, rr)
                rel = np.where(rr > 0, rr/sc, -1.0); dirn = int(np.argmax(rel))
            self.emit(dict(type='split', path=path, root=root_id, dirn=dirn, knew=int(knew)))
            for s, tag in ((-1, '0'), (1, '1')):
                q = p.copy(); q[dirn] += s*rr[dirn]/2; r2 = rr.copy(); r2[dirn] /= 2
                stack.append((path + tag, q, r2, depth + 1, knew, dict(xi=nxi, age=nage, ik=ik, incc=bool(inc_close),
                                                                    c=(ccon.tolist() if ccon is not None else None))))
        ok = counts['FAIL'] == 0
        item = dict(type='item_done', root=root_id, path=path0, ok=ok, counts=counts, handed_off=nhand)
        self.emit(item); self.items_done[(root_id, path0)] = item
        if path0 == root_id and nhand == 0 and not any(rt == root_id for rt, _ in self.handed):
            rec = dict(type='root_done', root=root_id, ok=ok, time=time.time() - t0, counts=counts,
                       box=dict(p0=[float(x) for x in p0], r=[float(x) for x in r0]))
            self.emit(rec); self.done[root_id] = rec
        return dict(ok=ok, handed=nhand, counts=counts)
