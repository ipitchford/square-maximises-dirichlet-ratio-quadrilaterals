"""Re-encode the cover5 manifests with the EXACT grid encoding of run_fam.axis_cells (adversarial review 2026-10-02, FATAL 3):
every cell (root or skipped) gets p0, r with exact float ends; the root / skipped classification and the reasons of the frozen
manifests are kept.  Checks (exit 1 on any failure):
  * the cell ids are exactly the grid cells;
  * every processed root keeps its dyadic snapping (cover2.dyadic_root of the new encoding == of the old one), so the existing
    logs (whose root_done boxes are the snapped boxes) stay consistent with the manifest;
  * every skipped cell's reason still holds on its new exact box (family predicate, current code).
Usage: remanifest.py FAMILY [FAMILY ...] [--write]   (without --write: check only)"""
import sys, json, itertools
sys.path.insert(0, '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts')
import numpy as np
import run_fam, fam5
from grid import dyadic_root
D = '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts/manifests/'

def main():
    fams = [a for a in sys.argv[1:] if not a.startswith('--')]; write = '--write' in sys.argv
    allok = True
    for name in fams:
        fn = D + name.lower() + '.json'
        man = json.load(open(fn))
        fam = fam5.FAMILIES[name]()
        axes = [run_fam.axis_cells(*g) for g in run_fam.GRIDS[name]]
        cells = {}
        for idx in itertools.product(*[range(len(x)) for x in axes]):
            rid = name + '_' + '_'.join(map(str, idx))
            cells[rid] = ([axes[i][idx[i]][0] for i in range(4)], [axes[i][idx[i]][1] for i in range(4)])
        ids = [c[0] for c in man['roots']] + [c[0] for c in man['skipped']]
        bad = []
        if sorted(ids) != sorted(cells) or len(set(ids)) != len(ids): bad.append('cell ids differ from the grid')
        roots = []; skipped = []; changed = 0
        for rid, p0, r in man['roots']:
            np0, nr = cells[rid]
            if dyadic_root(np0, nr) != dyadic_root(p0, r): bad.append(f'root {rid}: snapped box changes')
            changed += (np0 != p0 or nr != r)
            roots.append([rid, np0, nr])
        nfail = 0
        for rid, p0, r, why in man['skipped']:
            np0, nr = cells[rid]
            changed += (np0 != p0 or nr != r)
            if not fam.check(why, np.array(np0, float), np.array(nr, float)):   # RIGOR: exact: float64 arrays of Python floats (no rounding)
                nfail += 1; bad.append(f'skipped {rid}: {why!r} fails on the exact box {np0} +- {nr}')
            skipped.append([rid, np0, nr, why])
        print(f'{name}: {len(roots)} roots, {len(skipped)} skipped, {changed} cells re-encoded, {nfail} skipped predicates fail, '
              f'{len(bad)} problems', flush=True)
        for b in bad[:10]: print('   ', b)
        if bad: allok = False
        elif write:
            json.dump(dict(family=name, grid=man['grid'], roots=roots, skipped=skipped), open(fn, 'w'))
            print(f'   written {fn}')
    sys.exit(0 if allok else 1)

if __name__ == '__main__':
    main()
