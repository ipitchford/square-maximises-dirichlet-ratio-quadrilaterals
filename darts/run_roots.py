"""Run cover5 on given root ids of a family manifest into a separate log (timing tests of known-hard roots).
Usage: run_roots.py FAMILY LOG ROOT_ID [ROOT_ID ...]"""
import sys, json, time
sys.path.insert(0, '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts')
import fam5, run_fam
from cover5 import Cover5
name, log = sys.argv[1], sys.argv[2]
man = json.load(open(f'/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts/manifests/{name.lower()}.json'))
roots = {rid: (p0, r) for rid, p0, r in man['roots']}
cv = Cover5(log, fam5.FAMILIES[name](), fem_problem=run_fam.fem_problem(name))
for rid in sys.argv[3:]:
    t0 = time.time(); c0 = dict(cv.stats)
    ok = cv.process_root(rid, *roots[rid])
    d = {k: cv.stats[k] - c0[k] for k in ('cert', 'discard', 'FAIL', 'bo', 'fem', 'inc')}; d.update({k: round(cv.stats[k] - c0[k], 1) for k in ('t_bo', 't_inc', 't_fem')})
    print(json.dumps(dict(root=rid, ok=ok, dt=round(time.time() - t0, 1), **d)), flush=True)
cv.close()
