"""validate_fem.py — 2-D P2 FEM reference eigenvalues for random band shapes of both charts (BAND_SPEC §5 (a)).
Author: Anonymous.  Writes validation_fem.csv; validate_bounds.jl then compares the rigorous 1-D bounds against it.

Shapes: N random valid in-band points per chart (uniform in the chart root, exact rationals with denominator 10^6,
rejection to: convex, all vertex distances <= 1, 1/10 < s <= 2/5), plus near-worst shapes (thin-study adversarial
maximisers at caps 0.40/0.44 mapped into the charts, and near-triangles at the Lemma A boundary).
FEM: ../explore/quadfem.py eigs(V, k=6) (conforming P2, an upper approximation of lambda_1, lambda_2); k=7 recheck on
a subset to estimate the FEM error."""
import sys, time, csv
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, "../explore")
from quadfem import eigs

rng = np.random.default_rng(20261001)
NPER = int(sys.argv[1]) if len(sys.argv) > 1 else 60
DEN = 10**6

def verts(chart, p):
    if chart == 1:
        y3, t, x1, x3 = p
        return [(Fr(0), Fr(0)), (x1, -t * y3), (Fr(1), Fr(0)), (x3, y3)]
    y2, t, x2, x3 = p
    return [(Fr(0), Fr(0)), (Fr(1), Fr(0)), (x2, y2), (x3, t * y2)]

def cvals(V):
    out = []
    for i in range(4):
        a, b, c = V[i], V[(i + 1) % 4], V[(i + 2) % 4]
        out.append((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
    return out

def width2(V):
    c = cvals(V)
    return min(max(c[i], c[(i + 3) % 4])**2 / ((V[(i + 1) % 4][0] - V[i][0])**2 + (V[(i + 1) % 4][1] - V[i][1])**2) for i in range(4))

def valid(chart, p):
    V = verts(chart, p)
    if min(cvals(V)) <= 0: return False
    for i in range(4):
        for j in range(i + 1, 4):
            if (V[i][0] - V[j][0])**2 + (V[i][1] - V[j][1])**2 > 1: return False
    if chart == 2 and not V[3][0] < V[2][0]: return False
    return True

def lemmaA(V):
    c = cvals(V)
    return any(c[i] <= Fr(7, 200) * c[(i + 2) % 4] for i in range(4))

ROOT = {1: ([Fr(1, 20), 0, 0, 0], [Fr(2, 5), 1, 1, Fr(1, 2)]), 2: ([Fr(1, 10), 0, 0, 0], [Fr(2, 5), 1, 1, 1])}

def to_chart(Vin):
    """float vertices (any order/orientation, convex) -> (chart, p) via the reductions of BAND_SPEC §1"""
    V = [np.array(v, float) for v in Vin]
    a2 = sum(V[k][0] * V[(k + 1) % 4][1] - V[(k + 1) % 4][0] * V[k][1] for k in range(4))
    if a2 < 0: V = V[::-1]
    best = max(((np.linalg.norm(V[i] - V[j]), i, j) for i in range(4) for j in range(i + 1, 4)))
    _, i, j = best
    def tr(z, A, B):
        w = complex(*(z - A)) / complex(*(B - A)); return w.real, w.imag
    if (j - i) % 4 == 2:
        Pb = tr(V[(i + 1) % 4], V[i], V[j]); Pa = tr(V[(i + 3) % 4], V[i], V[j])
        x1, y1 = Pb[0], -Pb[1]; x3, y3 = Pa
        if y1 > y3: x1, y1, x3, y3 = x3, y3, x1, y1
        if x3 > 0.5: x1, x3 = 1 - x1, 1 - x3
        return 1, (y3, y1 / y3, x1, x3)
    a, b = (i, j) if (j - i) % 4 == 1 else (j, i)
    P2 = tr(V[(b + 1) % 4], V[a], V[b]); P3 = tr(V[(b + 2) % 4], V[a], V[b])
    x2, y2 = P2; x3, y3 = P3
    if y3 > y2: x2, y2, x3, y3 = 1 - x3, y3, 1 - x2, y2
    return 2, (y2, y3 / y2, x2, x3)

rows = []
def add(chart, p, tag):
    V = verts(chart, p)
    s = float(width2(V))**0.5
    Vf = [(float(v[0]), float(v[1])) for v in V]
    t = time.time(); w = eigs(Vf, k=6, nev=3); dt = time.time() - t
    rows.append(dict(chart=chart, tag=tag, p=[f"{x.numerator}/{x.denominator}" for x in p], s=s, lemmaA=int(lemmaA(V)),
                     l1=w[0], l2=w[1], l1k7="", l2k7="", sec=dt))

t0 = time.time()
for chart in (1, 2):
    lo, hi = ROOT[chart]; n = 0
    while n < NPER:
        p = tuple(Fr(int(round((float(lo[k]) + (float(hi[k]) - float(lo[k])) * rng.random()) * DEN)), DEN) for k in range(4))
        if not all(lo[k] <= p[k] <= hi[k] for k in range(4)) or p[1] == 0: continue
        if not valid(chart, p): continue
        s2 = width2(verts(chart, p))
        if not (Fr(1, 100) < s2 <= Fr(4, 25)): continue
        add(chart, p, "random"); n += 1
# near-worst shapes: thin-study adversarial maximisers (diameter-only 1-D bound 2.2250 at cap 0.40) and others
adv = [[(-0.471, -0.2197), (0.5282, -0.1782), (0.0419, 0.2019), (-0.099, 0.196)],
       [(-0.43707, -0.2022), (0.56293, -0.2022), (-0.03331, 0.20165), (-0.09255, 0.20275)],
       [(-0.53926, -0.35497), (0.29595, 0.14644), (0.31775, 0.16032), (-0.07444, 0.04822)],
       [(-0.50363, -0.02362), (0.00363, -0.02978), (0.49637, -0.02362), (0.00363, 0.07702)]]
for V in adv:
    chart, pf = to_chart(V)
    p = tuple(Fr(int(round(x * DEN)), DEN) for x in pf)
    if valid(chart, p): add(chart, p, "adversarial")
# shrink the adversarial cap-0.40 shape slightly into the band and add Lemma-A-boundary near-triangles
for scale in (0.97, 0.99, 0.995):
    chart, pf = to_chart(adv[0]); p = (Fr(int(round(pf[0] * scale * DEN)), DEN),) + tuple(Fr(int(round(x * DEN)), DEN) for x in pf[1:])
    if valid(chart, p) and width2(verts(chart, p)) <= Fr(4, 25): add(chart, p, f"adv0.40x{scale}")
for (y2, x2, x3c) in ((Fr(39, 100), Fr(1, 2), Fr(36, 1000)), (Fr(3, 10), Fr(7, 10), Fr(36, 1000)), (Fr(2, 5), Fr(3, 10), Fr(40, 1000))):
    t = Fr(1)  # P3 near the segment P0P2 so that c2 ~ 0.036 c0 (just outside Lemma A)
    x3 = t * x2 - x3c
    p = (y2, t * Fr(99, 100), x2, x3)
    if valid(2, p): add(2, p, "nearLemmaA")
# FEM error estimate: k=7 on every 10th row
for r in rows[::10]:
    V = verts(r["chart"], tuple(Fr(x) for x in r["p"]))
    w = eigs([(float(v[0]), float(v[1])) for v in V], k=7, nev=3)
    r["l1k7"], r["l2k7"] = w[0], w[1]
with open("validation_fem.csv", "w", newline="") as fh:
    wr = csv.writer(fh)
    wr.writerow(["chart", "tag", "p1", "p2", "p3", "p4", "s", "lemmaA", "lambda1_k6", "lambda2_k6", "lambda1_k7", "lambda2_k7"])
    for r in rows:
        wr.writerow([r["chart"], r["tag"], *r["p"], f"{r['s']:.6f}", r["lemmaA"], repr(float(r["l1"])), repr(float(r["l2"])), (repr(float(r["l1k7"])) if r["l1k7"] != "" else ""), (repr(float(r["l2k7"])) if r["l2k7"] != "" else "")])
print(f"{len(rows)} shapes written to validation_fem.csv in {time.time() - t0:.0f} s")
