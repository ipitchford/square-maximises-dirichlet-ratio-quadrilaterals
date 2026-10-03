"""Spike / collapse face for C in K1 = {angle CBD >= 8 deg, angle CDB >= 8 deg, |C| <= 6} (README, Theorem S2).
   COLLAR : tan b_A <= 0.02, every A on such a ray (any length, also a_x < 1)  -- chart SPIKE at t = 0.02 (Lemmas 9.1, 9.2)
   SPIKE  : 0.02 <= tan b_A <= 0.3 and a_x >= 1 (any length)                  -- chart SPIKE2, X0 = 1 (Lemma 9.4)
C roots: a dyadic quadtree of the rectangle [-6, 6] x [-6, 0], refined until the cell half-width is <= REL * (distance of the
cell to B and to D) or the cell lies outside K1 (rigorous predicate; such cells are listed in the root manifest and re-checked by
replay).  SPIKE roots: C cells x geometric t slabs (ratio 1.5).
Boxes are certified, discarded (outside K1), or deferred (float xi_h < XI_DEFER: a thin configuration handed to the 1-D
Born-Oppenheimer cover).  Resumable (cover3).  Usage:
   run_spike_all.py [K1|KF] COLLAR|SPIKE LOG SHARD NSHARDS  (production; region K1 by default)
   run_spike_all.py COLLAR|SPIKE LOG pilot N SEED           (random roots, statistics only)
   run_spike_all.py COLLAR|SPIKE manifest OUT.json          (root list + outside cells, for replay)"""
import sys, json, time, math, numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts')
from regions import k1_discard, sector_discard

XI_DEFER = 1.25
T1, T2, TRATIO = 0.02, 0.3, 1.5
# C regions: K1 = {angle CBD >= 8, angle CDB >= 8, |C| <= 6};  KF = {38 <= b_C <= 100 deg, 0.5 <= |C| <= 1.8} (the high-xi zone)
def _k1x_disc():
    from regions import P_C_in_sector
    d1 = k1_discard(8, 6.0)
    def d(p0, r):
        why = d1(p0, r)
        if why: return why
        if P_C_in_sector(p0, r, 38, 100, 0.5, 1.8): return 'C inside sector [38, 100] deg x [0.5, 1.8] (covered by run KF)'
        return None
    return d
REGIONS = {'K1': dict(disc=lambda: k1_discard(8, 6.0), rect=(-6.0, 6.0, -6.0, 0.0), REL=0.12, kmin=4, prefix=''),
           'K1X': dict(disc=_k1x_disc, rect=(-6.0, 6.0, -6.0, 0.0), REL=0.12, kmin=4, prefix='X', hybrid=True),
           'KF': dict(disc=lambda: sector_discard(38, 100, 0.5, 1.8), rect=(-2.0, 2.0, -2.0, 0.0), REL=0.10, kmin=5, prefix='F')}
REGION = 'K1'
def REG(): return REGIONS[REGION]

def c_cells():
    """dyadic quadtree of [-6,6] x [-6,0]; returns (roots, outside) lists of (cx0, cy0, hx, hy)"""
    roots = []; outside = []
    disc = REG()['disc'](); REL = REG()['REL']
    x0, x1, y0, y1 = REG()['rect']
    # initial cells: squares of side (y1 - y0)/2 (exact binary for our rectangles)
    side = (y1 - y0)/2; stack = []
    nx = int(round((x1 - x0)/side)); ny = 2
    for i in range(nx):
        for j in range(ny):
            stack.append((x0 + side*(i + 0.5), y0 + side*(j + 0.5), side/2, side/2))
    while stack:
        cx, cy, hx, hy = stack.pop()
        p0 = [0.0, 0.0, cx, cy]; r = [0.0, 0.0, hx, hy]
        if disc(p0, r): outside.append((cx, cy, hx, hy)); continue
        # distance of the cell to B = (0,0) and D = (1,0)
        def dist(px, py):
            dx = max(abs(cx - px) - hx, 0.0); dy = max(abs(cy - py) - hy, 0.0); return math.hypot(dx, dy)
        d = min(dist(0.0, 0.0), dist(1.0, 0.0))
        if max(hx, hy) < 1e-4: raise RuntimeError(f'root quadtree did not terminate near ({cx}, {cy})')
        if max(hx, hy) <= REL*d:
            roots.append((cx, cy, hx, hy)); continue
        if hx >= hy:
            for s in (-1, 1): stack.append((cx + s*hx/2, cy, hx/2, hy))
        else:
            for s in (-1, 1): stack.append((cx, cy + s*hy/2, hx, hy/2))
    roots.sort(key=lambda c: (round(math.atan2(-c[1], c[0]), 6), round(math.hypot(c[0], c[1]), 6)))
    return roots, outside

def t_slabs():
    ts = [T1]
    while ts[-1] < T2: ts.append(min(ts[-1]*TRATIO, T2))
    return [(ts[i], ts[i + 1]) for i in range(len(ts) - 1)]

def all_roots(mode):
    cells, _ = c_cells()
    out = []
    if mode == 'COLLAR':
        for j, (cx, cy, hx, hy) in enumerate(cells):
            out.append((f"{REG()['prefix']}COL_{j}", [T1, 0.0, cx, cy], [0.0, 0.0, hx, hy]))
    else:
        for a, (t0, t1) in enumerate(t_slabs()):
            for j, (cx, cy, hx, hy) in enumerate(cells):
                out.append((f"{REG()['prefix']}SPK_{a}_{j}", [(t0 + t1)/2, 0.0, cx, cy], [(t1 - t0)/2, 0.0, hx, hy]))
    return out

def make_cover(mode, log):
    from cover3 import Cover3
    from cover4 import Cover4
    from charts import spike_problem, spike2_problem
    if REG().get('hybrid'):
        disc = REG()['disc'](); kmin = REG()['kmin']
        sc = lambda p: np.array([abs(p[0]) if mode == 'SPIKE' else 1.0, 1.0, max(abs(p[2]), abs(p[3]), 0.3), max(abs(p[2]), abs(p[3]), 0.3)])
        pf = spike_problem if mode == 'COLLAR' else (lambda k: spike2_problem(k, 1.0))
        return Cover4(log, pf, bo_kind='spike', bo_K=160, discard_fn=disc, defer_fn=None, kmin=kmin, kmax=6,
                      maxdepth=26, split_scale_fn=sc)
    disc = REG()['disc'](); kmin = REG()['kmin']
    def defer(p, rr, res, depth):
        xi = res.get('xi_h')
        if xi is not None and xi < XI_DEFER: return f'thin: float xi_h = {xi:.4f} < {XI_DEFER} (1-D cover)'
        return None
    if mode == 'COLLAR':
        return Cover3(log, spike_problem, discard_fn=disc, defer_fn=defer, kmin=kmin, kmax=6, maxdepth=24,
                      split_scale_fn=lambda p: np.array([1.0, 1.0, max(abs(p[2]), abs(p[3]), 0.3), max(abs(p[2]), abs(p[3]), 0.3)]))
    return Cover3(log, lambda k: spike2_problem(k, 1.0), discard_fn=disc, defer_fn=defer, kmin=kmin, kmax=6, maxdepth=28,
                  split_scale_fn=lambda p: np.array([abs(p[0]), 1.0, max(abs(p[2]), abs(p[3]), 0.3), max(abs(p[2]), abs(p[3]), 0.3)]))

if __name__ == '__main__':
    if sys.argv[1] in REGIONS: REGION = sys.argv.pop(1)
    mode, log = sys.argv[1], sys.argv[2]
    if log == 'manifest':
        cells, outside = c_cells()
        json.dump(dict(family=REGION + '_' + mode, mode=mode, region=REGION, t_slabs=t_slabs(), cells=cells, outside=outside,
                       roots=all_roots(mode)), open(sys.argv[3], 'w'))
        print(len(cells), 'cells', len(outside), 'outside', len(all_roots(mode)), 'roots'); sys.exit()
    roots = all_roots(mode)
    if sys.argv[3] == 'pilot':
        n, seed = int(sys.argv[4]), int(sys.argv[5])
        rng = np.random.default_rng(seed); pick = [roots[i] for i in rng.choice(len(roots), size=min(n, len(roots)), replace=False)]
    else:
        i, n = int(sys.argv[3]), int(sys.argv[4]); pick = [rt for j, rt in enumerate(roots) if j % n == i]
    cv = make_cover(mode, log)
    t0 = time.time()
    for rid, p0, r in pick:
        c0 = dict(cv.stats); tr = time.time()
        ok = cv.process_root(rid, p0, r)
        d = {k: cv.stats[k] - c0[k] for k in ('cert', 'discard', 'defer', 'FAIL', 'ldl', 'bo', 'fem') if k in cv.stats}
        print(json.dumps(dict(root=rid, ok=ok, p0=p0, r=r, dt=round(time.time() - tr, 2), **d, wall=round(time.time() - t0, 1))), flush=True)
    cv.close()
    print('DONE', json.dumps(cv.stats), round(time.time() - t0, 1), flush=True)
