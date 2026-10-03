import numpy as np, time
from proto import chartII_box_ratio
rng = np.random.default_rng(2); N = 150
# coords: y2 = s in (0.1,0.4]; t = y3/y2 in (0,1]; x2, x3 with 0<x3<x2<1
lo = np.array([0.1, 0.0, 0.0, 0.0]); hi = np.array([0.4, 1.0, 1.0, 1.0])
def pts(u): y2 = u[0]; y3 = u[1]*y2; return u[2], y2, u[3], y3
def valid(x2, y2, x3, y3):
    P = np.array([[0, 0], [1, 0], [x2, y2], [x3, y3]])
    if not (0 < x3 < x2 < 1): return False
    for i in range(4):
        for j in range(i+1, 4):
            if np.linalg.norm(P[i]-P[j]) > 1: return False
    c = [np.cross(P[(i+1) % 4]-P[i], P[(i+2) % 4]-P[i]) for i in range(4)]
    if min(c) <= 0: return False
    for i in range(4):
        if c[i] <= 0.035*c[(i+2) % 4]: return False   # Lemma A region excluded
    return True
def br(u, d):
    x2, y2, x3, y3 = pts(u)
    y2l, y2u = y2-d[0], y2+d[0]; tl, tu = u[1]-d[1], u[1]+d[1]
    return chartII_box_ratio(x2-d[2], x2+d[2], y2l, y2u, x3-d[3], x3+d[3], max(tl*y2l, 1e-9), tu*y2u, N)
acc, n_in, n_tot, worst, t0 = [], 0, 0, 0, time.time()
while n_in < 300:
    u = lo+(hi-lo)*rng.random(4); n_tot += 1
    if not valid(*pts(u)): continue
    n_in += 1; r0 = br(u, np.zeros(4)); worst = max(worst, r0)
    eps = 0.002; c = []
    for i in range(4):
        d = np.zeros(4); d[i] = eps; c.append(max((br(u, d)/r0-1)/eps, 0.05))
    L = 0.8*(2.5/r0-1); acc.append(16*np.prod(c)/L**4)
frac = n_in/n_tot
print("chart II: valid fraction %.3f, worst point ratio %.4f, est boxes %.3g (median-based %.3g), %.0f s"
      % (frac, worst, np.prod(hi-lo)*frac*np.mean(acc), np.prod(hi-lo)*frac*np.median(acc), time.time()-t0))
print("density percentiles", np.percentile(acc, [50, 90, 99, 100]))
