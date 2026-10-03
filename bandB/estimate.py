import numpy as np, time
from proto import chartI_box_ratio, width_chartI
rng = np.random.default_rng(1)
N = 150
lo = np.array([np.log(0.05), 0.035, 0.0, 0.0]); hi = np.array([np.log(0.4), 1.0, 1.0, 0.5])
vol_box = np.prod(hi-lo)
def chart(u):
    Y3 = np.exp(u[0]); y1 = u[1]*Y3; return u[2], y1, u[3], Y3
def valid(x1, y1, x3, y3):
    P = np.array([[0, 0], [x1, -y1], [1, 0], [x3, y3]])
    for i in range(4):
        for j in range(i+1, 4):
            if np.linalg.norm(P[i]-P[j]) > 1+1e-12: return False
    s = width_chartI(x1, y1, x3, y3)
    return 0.1 < s <= 0.4
def boxratio(u, d):
    x1, y1, x3, y3 = chart(u)
    Y3l, Y3u = np.exp(u[0]-d[0]), np.exp(u[0]+d[0])
    tl, tu = u[1]-d[1], u[1]+d[1]
    return chartI_box_ratio(x1-d[2], x1+d[2], tl*Y3l, tu*Y3u, x3-d[3], x3+d[3], Y3l, Y3u, N)
acc, n_in, n_tot, worst, t0 = [], 0, 0, 0, time.time()
while n_in < 300:
    u = lo + (hi-lo)*rng.random(4); n_tot += 1
    x1, y1, x3, y3 = chart(u)
    if not valid(x1, y1, x3, y3): continue
    n_in += 1
    r0 = boxratio(u, np.zeros(4)); worst = max(worst, r0)
    eps = 0.004; c = []
    for i in range(4):
        d = np.zeros(4); d[i] = eps
        c.append(max((boxratio(u, d)/r0 - 1)/eps, 0.05))
    L = 0.8*(2.5/r0 - 1)
    acc.append(16*np.prod(c)/L**4)
frac = n_in/n_tot
est = vol_box*frac*np.mean(acc)
print("valid fraction %.3f, worst point ratio %.4f, mean c-product density -> est boxes %.3g (median-based %.3g), %.0f s"
      % (frac, worst, est, vol_box*frac*np.median(acc), time.time()-t0))
q = np.percentile(acc, [50, 90, 99, 100]); print("density percentiles", q)
