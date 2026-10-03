"""Cost benchmark (not part of the proof): run the first N roots of a family in the pool's pseudo-random order (an unbiased sample,
whatever their pool status) with the current cover code into a pilot log, and report seconds and certificates per root.
Usage: bench.py FAM N [OUTDIR]"""
import sys, os, json, time
HERE = '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts'
sys.path.insert(0, HERE)
import pool, fam5, run_fam
from cover5 import Cover5
fam, N = sys.argv[1], int(sys.argv[2])
out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(HERE, 'pilot_runs', 'bench')
os.makedirs(out, exist_ok=True)
items = [it for it in pool.order(['PP', 'SP_spk', 'SP_col', 'PF', 'FP']) if it[0] == fam][:N]
cv = Cover5(os.path.join(out, f'{fam.lower()}.jsonl'), fam5.FAMILIES[fam](), fem_problem=run_fam.fem_problem(fam))
tot = 0.0; n = 0
for f_, rid, p0, r in items:
    c0 = dict(cv.stats); t = time.time()
    ok = cv.process_root(rid, p0, r); dt = time.time() - t; tot += dt; n += 1
    d = {k: (round(cv.stats[k] - c0[k], 1) if isinstance(cv.stats[k], float) else cv.stats[k] - c0[k]) for k in cv.stats}
    print(json.dumps(dict(fam=fam, root=rid, ok=ok, dt=round(dt, 1), **d)), flush=True)
cv.close()
print('BENCH', json.dumps(dict(fam=fam, roots=n, total_s=round(tot, 1), mean_s=round(tot/max(n, 1), 1))), flush=True)
