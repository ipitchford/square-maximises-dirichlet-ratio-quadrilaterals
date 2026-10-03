"""xlocal.py — separately written check of the published records of the local step (Theorem 2).  Author: Anonymous.

Input: local/replay_v2/schur_cover2.json (written by local/cover2.py).  Exact rational arithmetic only (fractions.Fraction);
no code is shared with local/.  Checked:
  (1) the direction boxes lie on the eight faces of [-1,1]^4 and tile every face exactly (volumes sum to 8 per face and the
      interiors are pairwise disjoint), so every direction E with |E|_inf = 1 lies in some box;
  (2) every box certifies s in (0, 1/20] by ten s-intervals [q/200, (q+1)/200], q = 0..9, and |p|_2 <= 2s <= 1/10 = rho_0;
  (3) for every box and s-interval the logged shift level beta is one of the 17 levels and is admissible, i.e. beta <= lambda_2(Q_p)
      for every p = sE in the box, by one of two bounds re-derived here from the vertex formulas of paper §2.1:
        - circumscribed rectangle: Q_p lies in a rectangle of width 1 - a + max(|b|,|d|) and height 1 + a + max(|c|,|d|),
          so lambda_2(Q_p) >= pi^2 min(4/W^2 + 1/H^2, 1/W^2 + 4/H^2) (domain monotonicity), evaluated with an exact rational
          lower bound for pi^2 and exact rational upper bounds for W, H over the box;
        - Krahn-Szego: lambda_2(Q_p) >= 2 pi j01^2 / |Q_p| >= 2 pi j01^2 > 36 since |Q_p| = 1 - a^2 - d^2 <= 1;
  (4) the criterion is 'det' or 'tr' and the logged exact upper bound of det S/s^2 or tr S/s^2 is a negative rational
      (and agrees with the logged float to 1e-9 relative).
It does NOT recompute the Taylor-model enclosures themselves; local/replay_local_v2.sh does that.
The width/height formulas are first checked against the vertex coordinates on random exact points.
Exit status 0 iff every check passes.  Usage: xlocal.py [path/to/schur_cover2.json]"""
import sys, json, random, hashlib
from fractions import Fraction as F

PATH = sys.argv[1] if len(sys.argv) > 1 else "local/replay_v2/schur_cover2.json"
PI2_LO = F(986960440108, 10**11)          # pi^2 = 9.8696044010893586...; this rational is below it
LEVELS = {F(x) for x in ("36","37","38","39","40","41","42","43","44","45","46","93/2","47","95/2","48","97/2","49")}
problems = []

def verts(a, b, c, d):                    # paper §2.1 (Endo-Osting)
    h = F(1, 2)
    return [(h*(a+d-b), h*(-a+d+c)), (1+h*(-a+d+b), h*(-a-d-c)), (1+h*(-a-d-b), 1+h*(a-d+c)), (h*(a-d+b), 1+h*(a+d-c))]

# (0) the width/height formulas against the vertices, on random exact points
rng = random.Random(20261002)
for _ in range(20000):
    a, b, c, d = (F(rng.randint(-1000, 1000), 20000) for _ in range(4))
    v = verts(a, b, c, d); xs = [p[0] for p in v]; ys = [p[1] for p in v]
    if max(xs)-min(xs) != 1 - a + max(abs(b), abs(d)) or max(ys)-min(ys) != 1 + a + max(abs(c), abs(d)):
        problems.append(f"width/height formula fails at {(a, b, c, d)}"); break

def amax(lo, hi): return max(abs(lo), abs(hi))

raw = open(PATH, "rb").read()
data = json.loads(raw)
boxes = data["boxes"]
by_face = {}
for k, bx in enumerate(boxes):
    i, sgn = bx["face"]
    if i not in (0, 1, 2, 3) or sgn not in (1, -1): problems.append(f"box {k}: bad face {bx['face']}"); continue
    lo = [F(x) for x in bx["lo"]]; hi = [F(x) for x in bx["hi"]]
    if len(lo) != 3 or any(not (-1 <= l < h <= 1) for l, h in zip(lo, hi)): problems.append(f"box {k}: bad extent"); continue
    by_face.setdefault((i, sgn), []).append((lo, hi))
    smax = F(bx["smax"])
    if smax != F(1, 20): problems.append(f"box {k}: smax {smax} != 1/20")
    if 4 * smax**2 > F(1, 100): problems.append(f"box {k}: |p|_2 bound exceeds rho_0")    # |E|_2^2 <= 4 on the cube surface
    used = bx["used"]
    if len(used) != 10: problems.append(f"box {k}: {len(used)} s-intervals, expected 10"); continue
    E_lo = [None]*4; E_hi = [None]*4; E_lo[i] = E_hi[i] = F(sgn)
    for j, l, h in zip([j for j in range(4) if j != i], lo, hi): E_lo[j], E_hi[j] = l, h
    # coefficient of s in W and H, maximised over the box: W = 1 + s(-Ea + max(|Eb|,|Ed|)), H = 1 + s(Ea + max(|Ec|,|Ed|))
    kW = -E_lo[0] + max(amax(E_lo[1], E_hi[1]), amax(E_lo[3], E_hi[3]))
    kH = E_hi[0] + max(amax(E_lo[2], E_hi[2]), amax(E_lo[3], E_hi[3]))
    for q, (beta, crit, ub, ubf) in enumerate(used):
        slo, shi = F(q, 200), F(q + 1, 200)
        W = 1 + max(slo*kW, shi*kW); H = 1 + max(slo*kH, shi*kH)
        rect = PI2_LO * min(4/W**2 + 1/H**2, 1/W**2 + 4/H**2)
        beta = F(beta)
        if beta not in LEVELS: problems.append(f"box {k} s-int {q}: unknown level {beta}")
        if not (beta <= 36 or beta <= rect): problems.append(f"box {k} s-int {q}: level {beta} exceeds both lower bounds (rect {float(rect):.6f})")
        if crit not in ("det", "tr"): problems.append(f"box {k} s-int {q}: unknown criterion {crit}")
        ub = F(ub)
        if not ub < 0: problems.append(f"box {k} s-int {q}: bound {ub} is not negative")
        if abs(float(ub) - ubf) > 1e-9 * abs(ubf): problems.append(f"box {k} s-int {q}: float mismatch")

# (1) exact tiling of every face
for face in [(i, s) for i in range(4) for s in (1, -1)]:
    bl = by_face.get(face, [])
    vol = sum((h[0]-l[0])*(h[1]-l[1])*(h[2]-l[2]) for l, h in bl)
    if vol != 8: problems.append(f"face {face}: volume {vol} != 8")
    for x in range(len(bl)):
        for y in range(x + 1, len(bl)):
            (l1, h1), (l2, h2) = bl[x], bl[y]
            if all(max(l1[t], l2[t]) < min(h1[t], h2[t]) for t in range(3)): problems.append(f"face {face}: boxes {x},{y} overlap")
nb = len(boxes); nf = len(by_face)
crit = {"det": sum(u[1] == "det" for b in boxes for u in b["used"]), "tr": sum(u[1] == "tr" for b in boxes for u in b["used"])}
print(f"xlocal: {PATH} sha256 {hashlib.sha256(raw).hexdigest()}")
print(f"xlocal: {nb} boxes on {nf} faces; criteria {crit}; summary in file: {json.dumps(data['summary'].get('certified_boxes'))} boxes, "
      f"{data['summary'].get('failures')} failures")
if problems:
    print(f"XLOCAL FAIL: {len(problems)} problems"); [print("  ", p) for p in problems[:30]]; sys.exit(1)
print("XLOCAL PASS: exact tiling of all 8 faces, s in (0,1/20] covered, every shift level admissible, every logged bound negative")
