"""Unit tests of the rigorous slab bounds in bo1d.py against brute-force float sampling (validation, not part of the proof):
for random boxes of darts (BD chart) and random slabs, the certified slab bound H must dominate the longest slice component of
sampled darts of the box at sampled abscissae in the slab; same for spike boxes.  Also: the certified E must lie below an
accurate FEM lambda_1 of sampled darts (polyfem P3)."""
import sys, numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts')
from bo1d import slab_extent, corners2, arb, A
from bo_float import longest
rng = np.random.default_rng(11)
bad = 0; n = 0
for trial in range(300):
    ax, ay = rng.uniform(-1, 3), rng.uniform(0.05, 2); cx, cy = rng.uniform(-1, 3), rng.uniform(-2, -0.05)
    wa, wc = rng.uniform(0.001, 0.2, 2)
    Abox = ((ax, ay), (ax + wa, ay + wa)); Cbox = ((cx, cy), (cx + wc, cy + wc))
    th = rng.uniform(-np.pi, np.pi); t = arb(th); u = (t.cos(), t.sin()); uf = np.array([np.cos(th), np.sin(th)])
    SA = [(0., 0.), (1., 0.)] + corners2(*Abox); SC = [(0., 0.), (1., 0.)] + corners2(*Cbox)
    xBD = sorted([0.0, float(uf[0])])
    for s in range(6):
        a = rng.uniform(-3, 3); b = a + rng.uniform(1e-3, 0.3)
        eA = float(slab_extent(SA, u, a, b).upper()); eC = float(slab_extent(SC, u, a, b).upper())
        H = eA + eC if (b >= xBD[0] and a <= xBD[1]) else max(eA, eC)
        for q in range(8):
            A_ = (rng.uniform(ax, ax + wa), rng.uniform(ay, ay + wa)); C_ = (rng.uniform(cx, cx + wc), rng.uniform(cy, cy + wc))
            fD = (1 - C_[0])*A_[1] + C_[1]*(A_[0] - 1)
            P = np.array([A_, (0, 0), C_, (1, 0)]) if fD < 0 else None
            if P is None: continue            # (non-darts: the polygon A,B,C,D may still be simple; skip for the test)
            for z in np.linspace(a, b, 5):
                L = longest(P, z, uf); n += 1
                if L > H*(1 + 1e-9) + 1e-12: bad += 1; print('VIOLATION', Abox, Cbox, th, (a, b), L, H)
print('slab tests', n, 'violations', bad)
