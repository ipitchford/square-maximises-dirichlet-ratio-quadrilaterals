"""Compact bulk (README, family CC): A in the square (0, R] x (0, R] with tan b_A >= 0.3, C in the square [-R, R] x [-R, 0),
b_A <= b_C (mirror).  4-D roots = products of A-cells and C-cells (dyadic quadtrees refined until the cell half-width is
<= REL * distance to {B, D}, or the cell is tiny near B).  Every box is certified (BO first, FEM second: cover4.Cover4),
discarded by a rigorous predicate, or bisected; boxes refined below the floor are LOUD failures.
Discards (regions.check_reason re-checks them in replay):
  not a dart (f_D >= 0) | Prop 5.1 | mirror (b_C < b_A on the box) | tan b_A < 0.3 on the box (spike families)
  | Lemma 5.2 (b_C >= 174.29 deg, theta_B >= b_C)
Usage:  run_bulk.py LOG SHARD NSHARDS | run_bulk.py LOG pilot N SEED | run_bulk.py manifest OUT.json"""
import sys, json, time, math, numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts')
from regions import P_not_dart, P_prop51, P_mirror, P_tanbA_lt, P_bC_ge

R = 3.0; REL = 0.15; TMIN_CELL = 0.01

def disc(p0, r):
    if p0[1] - r[1] <= 0 or p0[3] + r[3] >= 0: return None        # never discard a box touching the axis (cannot certify either)
    if P_not_dart(p0, r): return 'not a dart (f_D >= 0 on box)'
    if P_prop51(p0, r): return 'Prop 5.1 (beta <= 0.0339)'
    if P_mirror(p0, r): return 'mirror (b_C < b_A on the box)'
    if P_tanbA_lt(p0, r, 0.3): return 'tan b_A < 0.3 on the box (spike families)'
    if P_bC_ge(p0, r, 174.29): return 'Lemma 5.2 (b_C >= 174.29 deg, theta_B >= b_C)'
    return None

def quadtree(x0, x1, y0, y1, drop):
    """dyadic quadtree of [x0,x1] x [y0,y1] (square, power-of-two side); drop(cx, cy, hx, hy) -> reason or None"""
    out = []; outside = []
    stack = [((x0 + x1)/2, (y0 + y1)/2, (x1 - x0)/2, (y1 - y0)/2)]
    while stack:
        cx, cy, hx, hy = stack.pop()
        why = drop(cx, cy, hx, hy)
        if why: outside.append(((cx, cy, hx, hy), why)); continue
        def dist(px, py):
            return math.hypot(max(abs(cx - px) - hx, 0.0), max(abs(cy - py) - hy, 0.0))
        d = min(dist(0.0, 0.0), dist(1.0, 0.0))
        if max(hx, hy) <= REL*d or max(hx, hy) <= TMIN_CELL:
            out.append((cx, cy, hx, hy)); continue
        for sx in (-1, 1):
            for sy in (-1, 1):
                stack.append((cx + sx*hx/2, cy + sy*hy/2, hx/2, hy/2))
    out.sort(); return out, outside

def a_cells():
    def drop(cx, cy, hx, hy):
        p0 = [cx, cy, 0.5, -0.5]; r = [hx, hy, 0.0, 0.0]
        if P_tanbA_lt(p0, r, 0.3): return 'tan b_A < 0.3 on the box (spike families)'
        return None
    return quadtree(0.0, 4.0, 0.0, 4.0, lambda cx, cy, hx, hy: ('A outside the square' if cx - hx >= R or cy - hy >= R else drop(cx, cy, hx, hy)))

def c_cells():
    def drop(cx, cy, hx, hy):
        if cx - hx >= R or cx + hx <= -R or cy + hy <= -R: return 'C outside the square'
        p0 = [1.0, 1.0, cx, cy]; r = [0.0, 0.0, hx, hy]
        if P_bC_ge(p0, r, 174.29): return 'Lemma 5.2 (b_C >= 174.29 deg, theta_B >= b_C)'
        return None
    return quadtree(-4.0, 4.0, -8.0, 0.0, drop)

def all_roots():
    A, _ = a_cells(); C, _ = c_cells()
    roots = []; skipped = []
    for i, (ax, ay, ahx, ahy) in enumerate(A):
        for j, (cx, cy, chx, chy) in enumerate(C):
            p0 = [ax, ay, cx, cy]; r = [ahx, ahy, chx, chy]
            why = disc(p0, r)
            if why: skipped.append((f'B_{i}_{j}', p0, r, why)); continue
            roots.append((f'B_{i}_{j}', p0, r))
    return roots, skipped

def make_cover(log):
    from cover4 import Cover4
    from charts import bd_problem
    sc = lambda p: np.array([max(math.hypot(p[0], p[1]), 0.3)]*2 + [max(math.hypot(p[2], p[3]), 0.3)]*2)
    return Cover4(log, bd_problem, bo_kind='bd', bo_K=120, discard_fn=disc, defer_fn=None, kmin=4, kmax=6, maxdepth=30,
                  split_scale_fn=sc)

if __name__ == '__main__':
    if sys.argv[1] == 'manifest':
        roots, skipped = all_roots()
        A, Aout = a_cells(); C, Cout = c_cells()
        json.dump(dict(family='BD_HYB', region='bulk', R=R, REL=REL, roots=roots, skipped=skipped, a_cells=A, a_out=Aout,
                       c_cells=C, c_out=Cout), open(sys.argv[2], 'w'))
        print(len(A), 'A cells', len(C), 'C cells', len(roots), 'roots', len(skipped), 'skipped products'); sys.exit()
    log = sys.argv[1]
    roots, _ = all_roots()
    if sys.argv[2] == 'pilot':
        n, seed = int(sys.argv[3]), int(sys.argv[4])
        rng = np.random.default_rng(seed); pick = [roots[i] for i in rng.choice(len(roots), size=min(n, len(roots)), replace=False)]
    else:
        i, n = int(sys.argv[2]), int(sys.argv[3]); pick = [rt for j, rt in enumerate(roots) if j % n == i]
    cv = make_cover(log)
    t0 = time.time()
    for rid, p0, r in pick:
        c0 = dict(cv.stats); tr = time.time()
        ok = cv.process_root(rid, p0, r)
        d = {k: cv.stats[k] - c0[k] for k in ('cert', 'discard', 'FAIL', 'bo', 'fem')}
        print(json.dumps(dict(root=rid, ok=ok, p0=p0, r=r, dt=round(time.time() - tr, 2), **d, wall=round(time.time() - t0, 1))), flush=True)
    cv.close()
    print('DONE', json.dumps(cv.stats), round(time.time() - t0, 1), flush=True)
