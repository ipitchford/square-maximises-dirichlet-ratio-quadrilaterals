"""Fail-closed check of the cached chart data (refcache/chart_*.pkl) that every certificate reads (adversarial review 2026-10-02):
for every chart constructor of charts.py and every mesh level k, if its cache file exists, the chart is rebuilt FROM THE EXACT
RATIONAL ASSEMBLY (chart.Chart(..., fresh=True): exact reference meshes, exact element integrals, one rounding to doubles) and
every cached attribute (CR / P2 matrices of every piece, mass diagonals, DOF counts and free-DOF lists, boundary classifications,
edge-piece incidences, reference edge vectors, exact reference triangles) must be IDENTICAL to the rebuilt one (exact equality,
no tolerance).  Every chart_*.pkl in refcache must belong to some constructor/level (no unexplained cached data).  The rebuilt
charts also pass the conformity check (no T-junctions; boundary mesh edges only on unshared piece edges).
Usage: check_charts.py [KMAX]   (exit 0 iff all checks pass)"""
import sys, os, glob, time
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import numpy as np, pickle, collections
import charts
from chart import cache_file, CACHE
import chart

CONSTRUCTORS = {'BD': charts.bd_chart, 'SPIKE-CR': charts.spike_cr_chart, 'C-only': charts.c_only_chart,
                'SPIKE2-CR': charts.spike2_cr_chart, 'FAR-CR': charts.far_cr_chart, 'DEEP-CR': charts.deep_cr_chart,
                'FAR2-P2': charts.far2_p2_chart}

def file_of(cons, k):
    class _Stub:
        def __init__(self, pieces, k, bc=None, kappa=2, tag=None, fresh=False): self.fn = cache_file(pieces, k, bc, kappa, tag)
    orig = charts.Chart; charts.Chart = _Stub
    try: return cons(k).fn
    finally: charts.Chart = orig

def same(a, b, path='$'):
    """exact structural equality; returns None or a description of the first difference"""
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        if not (isinstance(a, np.ndarray) and isinstance(b, np.ndarray)): return f'{path}: array vs non-array'
        if a.shape != b.shape or a.dtype != b.dtype: return f'{path}: shape/dtype {a.shape} {a.dtype} vs {b.shape} {b.dtype}'
        if a.dtype.kind == 'f':
            if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))): return f'{path}: non-finite entries'   # RIGOR: exact: boolean finiteness test
        if not np.array_equal(a, b): return f'{path}: array entries differ'
        return None
    if isinstance(a, dict) or isinstance(b, dict):
        if not (isinstance(a, dict) and isinstance(b, dict)): return f'{path}: dict vs non-dict'
        if set(a) != set(b): return f'{path}: keys differ'
        for k in a:
            d = same(a[k], b[k], f'{path}[{k!r}]')
            if d: return d
        return None
    if isinstance(a, (list, tuple)) or isinstance(b, (list, tuple)):
        if type(a) != type(b) or len(a) != len(b): return f'{path}: sequence type/length differ'
        for i, (x, y) in enumerate(zip(a, b)):
            d = same(x, y, f'{path}[{i}]')
            if d: return d
        return None
    if type(a) != type(b) and not (isinstance(a, (int, np.integer)) and isinstance(b, (int, np.integer))):
        return f'{path}: type {type(a).__name__} vs {type(b).__name__}'
    if isinstance(a, float) and not (np.isfinite(a) and np.isfinite(b)): return f'{path}: non-finite'
    return None if a == b else f'{path}: {a!r} != {b!r}'

def conformity(ch):
    pe = collections.Counter()
    for pc in ch.pieces_meta:
        L = pc['labels']
        for a, b in ((0, 1), (1, 2), (2, 0)): pe[frozenset((L[a], L[b]))] += 1
    bad = sum(1 for i in range(ch.nedge) if ch.cr_bc[i] in ('D', 'N') and len(ch.edge_pieces[i]) != 1)
    two = sum(1 for p in ch.edge_pieces if len(p) == 2)
    expected = sum(1 for e, c in pe.items() if c == 2)*(2**ch.k)
    return two == expected and bad == 0, f'mesh edges on shared piece edges {two} (expected {expected}), boundary-on-shared {bad}'

def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    ok = True; matched = set(); t00 = time.time()
    if chart.TRUSTED is not None:
        inv = chart.verify_cache_inventory()
        for b in inv: print('CACHE INVENTORY:', b)
        if inv: ok = False
    for name, cons in CONSTRUCTORS.items():
        for k in range(2, kmax + 1):
            fn = file_of(cons, k)
            if not os.path.exists(fn): continue
            matched.add(os.path.abspath(fn)); t0 = time.time()
            if chart.TRUSTED is not None:
                cached = chart.load_pinned(fn)                     # verified against the pin before unpickling
                if cached is None: print(f'UNPINNED cache file {os.path.basename(fn)}: not loaded'); ok = False; continue
            else:
                cached = pickle.load(open(fn, 'rb'))               # development use (no pin): trusts the local cache
            fresh = cons(k, fresh=True)
            d = same(cached, fresh.cached_state())
            conf, cinfo = conformity(fresh)
            good = d is None and conf
            ok &= good
            print(f'{name} k={k} {os.path.basename(fn)}: rebuilt in {time.time() - t0:.1f} s; '
                  f'{"IDENTICAL" if d is None else "DIFFERS: " + d}; conformity {"OK" if conf else "FAIL"} ({cinfo})', flush=True)
    for fn in sorted(glob.glob(os.path.join(CACHE, 'chart_*.pkl'))):
        if os.path.abspath(fn) not in matched:
            k = int(fn.rsplit('_k', 1)[1].split('.')[0])
            if k > kmax: print(f'NOT CHECKED (level above KMAX {kmax}): {os.path.basename(fn)}'); ok = False
            else: print(f'UNEXPLAINED cache file: {os.path.basename(fn)}'); ok = False
    print(f'check_charts: {"OK" if ok else "FAILED"} ({time.time() - t00:.0f} s)')
    sys.exit(0 if ok else 1)

if __name__ == '__main__':
    main()
