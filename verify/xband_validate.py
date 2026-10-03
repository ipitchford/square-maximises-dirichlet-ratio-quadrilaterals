import sys, random, math
from fractions import Fraction as Fr
sys.path.insert(0, "/Users/admin/claude-sandbox/quad-ppw-2026-09-30/verify"); sys.path.insert(0, "/Users/admin/claude-sandbox/quad-ppw-2026-09-30/explore")
import numpy as np
from xband import certify_box
from quadfem import eigs
rng = random.Random(5)
def width(V):
    w = 1e9
    for i in range(4):
        e = V[(i+1) % 4]-V[i]; n = np.array([-e[1], e[0]])/np.hypot(*e)
        w = min(w, max(abs((V-V[i]) @ n)))
    return w
def diam(V): return max(np.linalg.norm(V[i]-V[j]) for i in range(4) for j in range(4))
def convex(V):
    return all(np.cross(V[(i+1) % 4]-V[i], V[(i+2) % 4]-V[i]) > 0 for i in range(4))
done = {"I": 0, "II": 0}; worst = {"lo": 1e9, "up": 1e9}
while min(done.values()) < 20:
    chart = "I" if done["I"] <= done["II"] else "II"
    r = lambda a, b: Fr(rng.uniform(a, b)).limit_denominator(1000)
    if chart == "I":
        y3 = r(0.06, 0.4); t = r(0.05, 1); x1 = r(0.05, 0.95); x3 = r(0.05, 0.5); y1 = t*y3
        V = np.array([[0, 0], [float(x1), -float(y1)], [1, 0], [float(x3), float(y3)]])
        box = {"y3": (y3, y3), "t": (t, t), "x1": (x1, x1), "x3": (x3, x3)}
    else:
        y2 = r(0.1, 0.4); t = r(0.05, 1); x2 = r(0.1, 0.95); x3 = r(0.05, float(x2)-0.02); y3 = t*y2
        V = np.array([[0, 0], [1, 0], [float(x2), float(y2)], [float(x3), float(y3)]])
        box = {"y2": (y2, y2), "t": (t, t), "x2": (x2, x2), "x3": (x3, x3)}
    if not convex(V) or diam(V) > 1 + 1e-12: continue
    s = width(V)/diam(V)
    if not (0.1 < s <= 0.4): continue
    ok, E, R2 = certify_box(chart, box)
    lam = eigs([tuple(v) for v in V], k=6, nev=2)
    l1, l2 = float(lam[0]), float(lam[1])
    lo_slack = (l1 - float(E))/l1; up_slack = (float(R2.mid()) - l2)/l2
    worst["lo"] = min(worst["lo"], lo_slack); worst["up"] = min(worst["up"], up_slack)
    done[chart] += 1
    print(f"{chart} s={s:.3f} lam1_FEM={l1:.4f} E={float(E):.4f} (slack {lo_slack:+.4f}) lam2_FEM={l2:.4f} R2={float(R2.mid()):.4f} (slack {up_slack:+.4f}) cert_ratio={float(R2.mid())/float(E):.4f} true_xi={l2/l1:.4f}", flush=True)
print("WORST relative slacks: lower", worst["lo"], "upper", worst["up"], "(both must be > 0)")
