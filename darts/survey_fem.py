"""FEM landscape (non-rigorous): xi over (thA, thC, alpha, beta), P2 graded k=4 (accuracy ~1e-3 on xi)."""
import sys, time, csv, numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts')
from dartgeom import dart
from fem_ref import xi
out = sys.argv[1]
degs = np.arange(5, 176, 10.0)
rows = []
t0 = time.time()
with open(out, 'w') as f:
    w = csv.writer(f); w.writerow(['thA', 'thC', 'alpha', 'beta', 'xi', 'l1', 'l2'])
    for a in degs:
        for c in degs:
            if c < a or a + c >= 178: continue   # symmetry thA<->thC (with alpha<->gamma) : keep thC >= thA, all D
            for be in (0.04, 0.08, 0.15, 0.25, 0.4, 0.55, 0.7, 0.85, 0.95):
                for fr in (0.02, 0.1, 0.3, 0.5, 0.7, 0.9, 0.98):
                    al = (1 - be)*fr
                    A, B, C, D = dart(np.radians(a), np.radians(c), al, be)
                    try:
                        r, ev = xi(A, B, C, D, k=4, kappa=2.5, deg=2)
                    except Exception as e:
                        r, ev = np.nan, [np.nan, np.nan]
                    w.writerow([a, c, al, be, r, ev[0], ev[1]])
            f.flush()
print('done', time.time() - t0)
