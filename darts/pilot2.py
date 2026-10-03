"""Sampling pilot for the rigorous cover: random root boxes in each region, full adaptive cover of each root, statistics.
Usage: pilot2.py REGION NROOTS SEED      REGION in {MAIN, SPIKE, COLLAR}
Root size: half-widths 5% of |A| (ax, ay) and of |C| (cx, cy); SPIKE: t +- 10%; COLLAR: t fixed = T1."""
import sys, time, json, numpy as np, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts')
from cover2 import Cover2
from charts import bd_problem, spike2_problem, spike_problem
from regions import main_discard, fD_bounds, beta_max, R_MAX, T1, T2, X0

def sample_C(rng):
    while True:
        u = rng.uniform(np.log(0.3), np.log(R_MAX)); b = rng.uniform(0.05, np.radians(150))
        C = np.exp(u)*np.array([np.cos(b), -np.sin(b)])
        if C[1] < -0.02: return C

def run(region, nroots, seed, tlimit=600):
    rng = np.random.default_rng(seed)
    out = f'/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts/pilot_runs/pilot2_{region}_{seed}.jsonl'
    if region == 'MAIN':
        cv = Cover2(out, bd_problem, discard_fn=main_discard, kmin=4, kmax=6,
                    split_scale_fn=lambda p: np.array([abs(p[1]), abs(p[1]), abs(p[3]), abs(p[3])]))
    elif region == 'SPIKE':
        cv = Cover2(out, lambda k: spike2_problem(k, X0), discard_fn=None, kmin=4, kmax=6, maxdepth=24,
                    split_scale_fn=lambda p: np.array([abs(p[0]), 1.0, abs(p[3]), abs(p[3])]))
    else:
        cv = Cover2(out, spike_problem, discard_fn=None, kmin=4, kmax=6, maxdepth=20,
                    split_scale_fn=lambda p: np.array([1.0, 1.0, abs(p[3]), abs(p[3])]))
    rows = []
    i = 0; tries = 0
    while i < nroots and tries < 100000:
        tries += 1
        C = sample_C(rng)
        if region == 'MAIN':
            u = rng.uniform(0, np.log(R_MAX)); b = rng.uniform(np.arctan(T2), np.pi - 0.05)
            A = np.exp(u)*np.array([np.cos(b), np.sin(b)])
            p0 = np.array([A[0], A[1], C[0], C[1]])
            if np.linalg.norm(C) > R_MAX: continue
            r = 0.05*np.array([np.linalg.norm(A)]*2 + [np.linalg.norm(C)]*2)
            if main_discard(p0, r*1e-6) is not None: continue
            if p0[1] - r[1] <= 0 or p0[3] + r[3] >= 0: continue
        elif region == 'SPIKE':
            t = np.exp(rng.uniform(np.log(T1), np.log(T2)))
            if np.arctan(t) > np.arctan2(-C[1], C[0]): continue              # mirror domain b_A <= b_C
            if fD_bounds([X0, X0*t, C[0], C[1]], [0]*4)[1] >= 0: continue     # must be a dart at a_x = X0
            p0 = np.array([t, 0.0, C[0], C[1]]); r = np.array([0.1*t, 0.0, 0.05*np.linalg.norm(C), 0.05*np.linalg.norm(C)])
            if p0[3] + r[3] >= 0: continue
        else:
            if np.arctan(T1) > np.arctan2(-C[1], C[0]): continue
            p0 = np.array([T1, 0.0, C[0], C[1]]); r = np.array([0.0, 0.0, 0.05*np.linalg.norm(C), 0.05*np.linalg.norm(C)])
            if p0[3] + r[3] >= 0: continue
        t0 = time.time(); c0 = dict(cv.stats)
        ok = cv.process_root(f'{region}_{seed}_{i}', p0, r)
        dt = time.time() - t0
        d = {k: cv.stats[k] - c0[k] for k in cv.stats}
        rows.append(dict(i=i, p0=p0.tolist(), r=r.tolist(), ok=ok, time=dt, **d))
        print(json.dumps(rows[-1]), flush=True)
        i += 1
    return rows

if __name__ == '__main__':
    run(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
