# Monte Carlo composition of the E3T (deferred thin band) leaves of a transport cover log:
# volume fractions that are convex, outside E2 (eta = 0.035), and with 1/10 < w/D <= 2/5 (the part a cover must certify).
import json, sys, numpy as np
V0 = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float)
DV = np.zeros((4, 2, 4))
DV[0, 0] = [0.5, -0.5, 0, 0.5]; DV[0, 1] = [-0.5, 0, 0.5, 0.5]
DV[1, 0] = [-0.5, 0.5, 0, 0.5]; DV[1, 1] = [-0.5, 0, -0.5, -0.5]
DV[2, 0] = [-0.5, -0.5, 0, -0.5]; DV[2, 1] = [0.5, 0, 0.5, -0.5]
DV[3, 0] = [0.5, 0.5, 0, -0.5]; DV[3, 1] = [0.5, 0, -0.5, 0.5]
def verts(P):  # P: (N,4) -> (N,4,2)
    return V0[None] + np.einsum('jki,ni->njk', DV, P)
def cr(ax, ay, bx, by): return ax * by - ay * bx
def classify(P):
    V = verts(P); N = len(P)
    c = np.stack([cr(V[:, (i+1) % 4, 0] - V[:, i, 0], V[:, (i+1) % 4, 1] - V[:, i, 1],
                     V[:, (i+2) % 4, 0] - V[:, i, 0], V[:, (i+2) % 4, 1] - V[:, i, 1]) for i in range(4)], 1)
    convex = (c > 0).all(1)
    e2 = np.zeros(N, bool)
    for i in range(4): e2 |= c[:, i] <= 0.035 * c[:, (i+2) % 4]
    ws = []
    for i in range(4):
        ex = V[:, (i+1) % 4, 0] - V[:, i, 0]; ey = V[:, (i+1) % 4, 1] - V[:, i, 1]; L = np.hypot(ex, ey)
        ws.append(np.max(np.stack([np.abs(cr(ex, ey, V[:, j, 0] - V[:, i, 0], V[:, j, 1] - V[:, i, 1])) for j in range(4)], 1), 1) / L)
    w = np.min(np.stack(ws, 1), 1)
    D = np.max(np.stack([np.hypot(V[:, i, 0] - V[:, j, 0], V[:, i, 1] - V[:, j, 1]) for i in range(4) for j in range(i+1, 4)], 1), 1)
    return convex, e2, w / D
recs = [json.loads(l) for l in open(sys.argv[1]) if '"E3T"' in l or '"COMPOUND"' in l]
boxes = [(np.array(r['lo']), np.array(r['hi'])) for r in recs if r['status'] == 'E3T']
vols = np.array([np.prod(h - l) for l, h in boxes]); tot = vols.sum()
rng = np.random.default_rng(1); M = 400000
idx = rng.choice(len(boxes), M, p=vols / tot)
lo = np.array([boxes[i][0] for i in idx]); hi = np.array([boxes[i][1] for i in idx])
P = lo + (hi - lo) * rng.random((M, 4))
convex, e2, wd = classify(P)
need = convex & ~e2 & (wd > 0.1)
print(f"E3T leaves {len(boxes)}, volume {tot:.4e}")
print(f"convex {convex.mean():.3f}; convex&~E2 {(convex&~e2).mean():.3f}; convex&~E2&w/D>0.1 {need.mean():.3f} -> volume {need.mean()*tot:.4e}")
print(f"  of which w/D>0.35: {(need & (wd > 0.35)).mean():.3f}; max w/D among E3T samples {wd[convex].max():.4f}")
for r in (0.0025, 0.005, 0.01, 0.02):
    print(f"  boxes of radius {r}: {need.mean()*tot/(2*r)**4:.3e}")
