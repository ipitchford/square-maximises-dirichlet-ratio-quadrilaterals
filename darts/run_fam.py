"""Runner for cover5 families.  Roots: a dyadic product grid over the family's parameter ranges (GRIDS below); roots entirely
discarded by a family predicate are listed in the manifest as 'skipped' (with the reason; replay re-checks them), the others are
processed (certified / discarded / bisected).  Resumable; sharded by root index.
Usage:  run_fam.py FAMILY manifest OUT.json
        run_fam.py FAMILY LOG SHARD NSHARDS
        run_fam.py FAMILY LOG pilot N SEED"""
import sys, os, json, time, math, itertools, numpy as np
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import fam5

# grid: per coordinate (lo, hi, step) with lo, hi, step dyadic; a coordinate with step None is inactive (value lo, radius 0)
TCELLS = [(0.02*1.5**i, min(0.02*1.5**(i + 1), 0.3)) for i in range(7)]       # spike t cells (geometric, 0.02 .. 0.3)   # RIGOR: exact: definition of the nominal cell ends (floats); every cell is then encoded exactly by exact_cell and coverage is checked in accounting.py
GRIDS = {
    'PP':     [(0.0, 1.75, 0.25), (-2.5, 3.0, 0.5), (0.0, 3.25, 0.25), (-4.0, 3.0, 0.5)],
    'SP_col': [('cells', [(0.0, 0.02)]), (0.0, 0.0, None), (0.0, 3.25, 0.25), (-4.0, 3.0, 0.5)],
    'SP_spk': [('cells', TCELLS), (0.0, 0.0, None), (0.0, 3.25, 0.25), (-4.0, 3.0, 0.5)],
    'PF':     [(0.0, 1.75, 0.25), (-2.5, 3.0, 0.5), (0.0, 3.25, 0.25), (0.0, 0.0, None)],
    'FP':     [(0.0, 1.75, 0.25), (0.0, 0.0, None), (0.0, 3.25, 0.25), (-4.0, 3.0, 0.5)],
    'FF':     [(0.0, 1.75, 0.25), (0.0, 0.0, None), (0.0, 3.25, 0.25), (0.0, 0.0, None)],
    'SF_col': [('cells', [(0.0, 0.02)]), (0.0, 0.0, None), (0.0, 3.25, 0.25), (0.0, 0.0, None)],
    'SF_spk': [('cells', TCELLS), (0.0, 0.0, None), (0.0, 3.25, 0.25), (0.0, 0.0, None)],
}
FEM = {'PP': 'bd', 'SP_col': 'collar', 'SP_spk': 'spike'}

T_SPIKE_MAX = 0.30000000000000004       # the smallest float >= 3/10: the spike families cover t <= 3/10 (PP discards tan b_A < 3/10)
def exact_cell(a, b, keep):
    """exact-endpoint float encoding (c, h) of a closed cell containing [a, b] (floats): c - h and c + h are exact floats with
    c - h <= a and c + h >= b (checked in rational arithmetic).  The end named by keep ('lo', 'hi' or 'both') is kept exactly,
    the other one is moved outward by ulps until the encoding is exact."""
    from fractions import Fraction as Fr
    A, B = float(a), float(b)   # RIGOR: exact: float() of floats
    for _ in range(64):
        c = (Fr(A) + Fr(B))/2; h = (Fr(B) - Fr(A))/2
        cf, hf = float(c), float(h)   # RIGOR: exact: accepted only if the conversion is exact (checked with Fractions)
        if Fr(cf) == c and Fr(hf) == h and Fr(cf - hf) == Fr(A) and Fr(cf + hf) == Fr(B): return cf, hf
        if keep == 'lo': B = math.nextafter(B, math.inf)   # RIGOR: exact: moves the end outward by one ulp (the result is checked exactly)
        elif keep == 'hi': A = math.nextafter(A, -math.inf)   # RIGOR: exact: moves the end outward by one ulp (the result is checked exactly)
        else: raise ValueError(f'cell [{a!r}, {b!r}] has no exact encoding with both ends kept')
    raise ValueError(f'exact_cell failed for [{a!r}, {b!r}]')
def axis_cells(lo, hi=None, step=None):
    """per-axis cells (c, h) of the root grid, every one with EXACT float ends c -/+ h (represented closed boxes).  Dyadic step
    grids are exact as they stand; the spike t-cells (non-dyadic ends 0.02 * 1.5^i) are encoded outward: the first cell keeps its
    lower end, the last cell is extended to T_SPIKE_MAX >= 3/10 and keeps it, every other cell keeps its lower end."""
    if lo == 'cells':
        cells = list(hi); out = []
        for i, (a, b) in enumerate(cells):
            if i == len(cells) - 1 and b == 0.3: out.append(exact_cell(a, T_SPIKE_MAX, 'hi'))
            elif a == 0.0: out.append(exact_cell(a, b, 'both'))
            else: out.append(exact_cell(a, b, 'lo'))
        return out
    if step is None: return [(lo, 0.0)]
    n = int(round((hi - lo)/step))   # RIGOR: exact: the number of cells of a dyadic grid (exact quotient); coverage is checked exactly in accounting.py
    return [(lo + step*(i + 0.5), step/2) for i in range(n)]

def roots_for(name):
    fam = fam5.FAMILIES[name]()
    axes = [axis_cells(*g) for g in GRIDS[name]]
    roots = []; skipped = []
    for idx in itertools.product(*[range(len(a)) for a in axes]):
        p0 = [axes[i][idx[i]][0] for i in range(4)]; r = [axes[i][idx[i]][1] for i in range(4)]
        rid = name + '_' + '_'.join(map(str, idx))
        why = fam.discard(np.array(p0), np.array(r))   # RIGOR: exact: float64 arrays of Python floats (no rounding)
        if why: skipped.append((rid, p0, r, why))
        else: roots.append((rid, p0, r))
    return fam, roots, skipped

def fem_problem(name):
    from charts import bd_problem, spike_problem, spike2_problem
    f = FEM.get(name)
    if f == 'bd': return bd_problem
    if f == 'collar': return spike_problem
    if f == 'spike': return lambda k: spike2_problem(k, 1.0)
    return None

if __name__ == '__main__':
    name = sys.argv[1]
    if sys.argv[2] == 'manifest':
        fam, roots, skipped = roots_for(name)
        json.dump(dict(family=name, grid=GRIDS[name], roots=roots, skipped=skipped), open(sys.argv[3], 'w'))
        print(name, len(roots), 'roots', len(skipped), 'skipped'); sys.exit()
    log = sys.argv[2]
    # production roots come from the FROZEN manifest (so that a later code change cannot change the root list of a run)
    fam = fam5.FAMILIES[name]()
    man = json.load(open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'manifests', f'{name.lower()}.json')))
    roots = [(rid, p0, r) for rid, p0, r in man['roots']]
    if sys.argv[3] == 'pilot':
        n, seed = int(sys.argv[4]), int(sys.argv[5])
        rng = np.random.default_rng(seed); pick = [roots[i] for i in rng.choice(len(roots), size=min(n, len(roots)), replace=False)]   # RIGOR: heuristic: pilot sampling
    else:
        i, n = int(sys.argv[3]), int(sys.argv[4]); pick = [rt for j, rt in enumerate(roots) if j % n == i]
    from cover5 import Cover5
    cv = Cover5(log, fam, fem_problem=fem_problem(name))
    t0 = time.time()
    for rid, p0, r in pick:
        c0 = dict(cv.stats); tr = time.time()
        ok = cv.process_root(rid, p0, r)
        d = {k: cv.stats[k] - c0[k] for k in ('cert', 'discard', 'FAIL', 'bo', 'fem')}
        print(json.dumps(dict(root=rid, ok=ok, p0=p0, r=r, dt=round(time.time() - tr, 2), **d, wall=round(time.time() - t0, 1))), flush=True)   # RIGOR: info: progress output
    cv.close()
    print('DONE', json.dumps(cv.stats), round(time.time() - t0, 1), flush=True)   # RIGOR: info: progress output
