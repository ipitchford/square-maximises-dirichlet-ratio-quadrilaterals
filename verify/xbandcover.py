"""Independent exact checker for the Lemma B' cover (own code). Exclusion predicates in exact rationals, per BAND_SPEC.md §1–2.
Chart I box:  y3 in [y3l,y3u], t=y1/y3 in [tl,tu], x1 in [x1l,x1u], x3 in [x3l,x3u]; P0=(0,0), P1=(x1,-y1), P2=(1,0), P3=(x3,y3).
Chart II box: y2 in [y2l,y2u], t=y3/y2 in [tl,tu], x2 in [x2l,x2u], x3 in [x3l,x3u]; P0=(0,0), P1=(1,0), P2=(x2,y2), P3=(x3,y3).
Every c_i equals (positive scale y3 or y2) x (multilinear polynomial in (t, x, x')), so its extremes over a box are attained at
corners: 'holds on the whole box' is decided exactly by corner evaluation."""
from fractions import Fraction as Fr
from itertools import product
ETA = Fr(7, 200)

def corners(*ivs): return product(*ivs)

# ---------- chart I: c_i / y3 as multilinear functions of (t, x1, x3) ----------
def cI(t, x1, x3):
    c0 = t                       # y1/y3
    c1 = 1 + t - x1 - x3*t       # (y3 + y1 - x1 y3 - x3 y1)/y3
    c2 = Fr(1)
    c3 = x3*t + x1
    return (c0, c1, c2, c3)
# ---------- chart II: c_i / y2 as multilinear functions of (t, x2, x3) ----------
def cII(t, x2, x3):
    c0 = Fr(1)
    c1 = 1 - x3 - t*(1 - x2)
    c2 = x2*t - x3
    c3 = t
    return (c0, c1, c2, c3)

def lemma_A(chart, B):
    f = cI if chart == "I" else cII
    ivs = (B["t"], B["x1"], B["x3"]) if chart == "I" else (B["t"], B["x2"], B["x3"])
    for i in range(4):
        if all(f(*c)[i] - ETA*f(*c)[(i+2) % 4] <= 0 for c in corners(*ivs)): return True
    return False

def nonconvex(chart, B):
    f = cI if chart == "I" else cII
    ivs = (B["t"], B["x1"], B["x3"]) if chart == "I" else (B["t"], B["x2"], B["x3"])
    return any(all(f(*c)[i] <= 0 for c in corners(*ivs)) for i in range(4))

def sqmin(l, h): return Fr(0) if l <= 0 <= h else min(l*l, h*h)
def dist2_min_exceeds_one(chart, B):
    """True if some pair of vertices is at distance > 1 at every point of the box (exact lower bounds of squared distances)."""
    if chart == "I":
        (yl, yu), (tl, tu), (x1l, x1u), (x3l, x3u) = B["y3"], B["t"], B["x1"], B["x3"]
        y1l, y1u = tl*yl, tu*yu
        lows = [sqmin(x1l, x1u) + y1l**2,                    # |P1|^2
                sqmin(1-x1u, 1-x1l) + y1l**2,                # |P1-P2|^2
                sqmin(x3l, x3u) + yl**2,                     # |P3|^2
                sqmin(1-x3u, 1-x3l) + yl**2,                 # |P3-P2|^2
                sqmin(x3l-x1u, x3u-x1l) + (yl + y1l)**2]     # |P1-P3|^2
    else:
        (yl, yu), (tl, tu), (x2l, x2u), (x3l, x3u) = B["y2"], B["t"], B["x2"], B["x3"]
        y3l = tl*yl
        lows = [sqmin(x2l, x2u) + yl**2, sqmin(1-x2u, 1-x2l) + yl**2,
                sqmin(x3l, x3u) + y3l**2, sqmin(1-x3u, 1-x3l) + y3l**2,
                sqmin(x2l-x3u, x2u-x3l) + Fr(0)]             # |P2-P3|^2 >= (x2-x3)^2 (y-difference may vanish)
    return any(v > 1 for v in lows)

def s_upper(chart, B):   # upper bound of width/diam on the box (diam = 1 in the chart)
    if chart == "I": return B["t"][1]*B["y3"][1] + B["y3"][1]          # s <= y1 + y3
    return B["y2"][1]                                                  # s = y2 exactly
def s_lower(chart, B):
    return B["y3"][0] if chart == "I" else B["y2"][0]                  # s >= y3 (chart I), s = y2 (chart II)

def check_exclusion(chart, B, status):
    if status == "LEMMA_A": return lemma_A(chart, B)
    if status == "LEMMA_B": return s_upper(chart, B) <= Fr(1, 10)
    if status == "OUT_OF_BAND": return s_lower(chart, B) > Fr(2, 5)
    if status == "OUTSIDE_CHART":
        if nonconvex(chart, B) or dist2_min_exceeds_one(chart, B): return True
        if chart == "II" and B["x3"][0] >= B["x2"][1]: return True
        if chart == "I" and B["x3"][0] > Fr(1, 2): return True
        if B["t"][0] > 1: return True
        return False
    raise ValueError(status)

if __name__ == "__main__":
    # self-test on obvious boxes
    B = {"y3": (Fr(1, 10), Fr(2, 10)), "t": (Fr(0), Fr(1, 100)), "x1": (Fr(1, 4), Fr(1, 2)), "x3": (Fr(1, 4), Fr(1, 2))}
    print("chart I tiny-t box is LEMMA_A:", lemma_A("I", B))
    B2 = {"y3": (Fr(1, 100), Fr(2, 100)), "t": (Fr(1, 2), Fr(1)), "x1": (Fr(1, 4), Fr(1, 2)), "x3": (Fr(1, 4), Fr(1, 2))}
    print("chart I very thin box is LEMMA_B:", check_exclusion("I", B2, "LEMMA_B"))
    B3 = {"y2": (Fr(2, 10), Fr(3, 10)), "t": (Fr(1, 2), Fr(1)), "x2": (Fr(1, 10), Fr(2, 10)), "x3": (Fr(6, 10), Fr(7, 10))}
    print("chart II box with x3 > x2 everywhere is OUTSIDE_CHART:", check_exclusion("II", B3, "OUTSIDE_CHART"))
    # cross-check the multilinear c_i formulas against direct cross products at random points
    import random
    rng = random.Random(1)
    def cross(o, a, b): return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])
    for _ in range(200):
        y3 = Fr(rng.randint(1, 400), 1000); t = Fr(rng.randint(1, 1000), 1000); x1 = Fr(rng.randint(1, 999), 1000); x3 = Fr(rng.randint(1, 999), 1000)
        P = [(0, 0), (x1, -t*y3), (1, 0), (x3, y3)]
        direct = [cross(P[i], P[(i+1) % 4], P[(i+2) % 4]) for i in range(4)]
        assert all(direct[i] == y3*cI(t, x1, x3)[i] for i in range(4)), "chart I c_i mismatch"
        y2 = y3; x2 = x1
        P = [(0, 0), (1, 0), (x2, y2), (x3, t*y2)]
        direct = [cross(P[i], P[(i+1) % 4], P[(i+2) % 4]) for i in range(4)]
        assert all(direct[i] == y2*cII(t, x2, x3)[i] for i in range(4)), "chart II c_i mismatch"
    print("c_i formulas verified exactly at 200 random points per chart")


# =====================================================================================================================
# Full independent check of the Lemma B' logs (fail-closed): exact tiling from path letters, exact exclusions, sampled
# certificates with verify/xband.py.  Box fields: chart I (y3, t, x1, x3); chart II (y2, t, x2, x3).
import json, sys, math, random
ROOTS = {1: [(Fr(1, 20), Fr(2, 5)), (Fr(0), Fr(1)), (Fr(0), Fr(1)), (Fr(0), Fr(1, 2))],
         2: [(Fr(1, 10), Fr(2, 5)), (Fr(0), Fr(1)), (Fr(0), Fr(1)), (Fr(0), Fr(1))]}
NAMES = {1: ("y3", "t", "x1", "x3"), 2: ("y2", "t", "x2", "x3")}
STATUSES = {"CERT", "LEMMA_A", "LEMMA_B", "OUTSIDE_CHART", "OUT_OF_BAND"}

def box_from_path(chart, path):
    B = [list(iv) for iv in ROOTS[chart]]
    for ch in path:
        k = "abcd".index(ch.lower()); m = (B[k][0] + B[k][1])/2
        if ch.islower(): B[k][1] = m
        else: B[k][0] = m
    return [tuple(iv) for iv in B]

def verts_corner(chart, z):
    if chart == 1:
        y3, t, x1, x3 = z; return [(Fr(0), Fr(0)), (x1, -t*y3), (Fr(1), Fr(0)), (x3, y3)]
    y2, t, x2, x3 = z; return [(Fr(0), Fr(0)), (Fr(1), Fr(0)), (x2, y2), (x3, t*y2)]

def cr(o, a, b): return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])

def width_bounds(chart, B):
    """rigorous (lower, upper) bounds of width(Q) over the box: width = min_i max_j dist(v_j, line e_i).
    Cross products are multilinear in the box coordinates (corners give exact ranges); |e_i| bounded via component ranges."""
    Cs = [verts_corner(chart, z) for z in product(*B)]
    lows, ups = [], []
    for i in range(4):
        ex = [V[(i+1) % 4][0]-V[i][0] for V in Cs]; ey = [V[(i+1) % 4][1]-V[i][1] for V in Cs]
        e2max = max(min(ex)**2, max(ex)**2) + max(min(ey)**2, max(ey)**2)
        e2min = sqmin(min(ex), max(ex)) + sqmin(min(ey), max(ey))
        dlo, dhi = [], []
        for j in range(4):
            if j in (i, (i+1) % 4): continue
            vals = [cr(V[i], V[(i+1) % 4], V[j]) for V in Cs]; lo_, hi_ = min(vals), max(vals)
            amin = Fr(0) if lo_ <= 0 <= hi_ else min(abs(lo_), abs(hi_)); amax = max(abs(lo_), abs(hi_))
            dlo.append((amin, e2max)); dhi.append((amax, e2min))
        # extent_i >= max_j amin/sqrt(e2max) ; extent_i <= max_j amax/sqrt(e2min)
        lows.append(max(a*a/e for a, e in dlo) if all(e > 0 for _, e in dlo) else Fr(0))          # squared
        ups.append(max(a*a/e for a, e in dhi) if all(e > 0 for _, e in dhi) else None)            # squared
    wlo2 = min(lows)                                            # width^2 >= min_i extent_i_lower^2
    wup2 = min(u for u in ups if u is not None) if any(u is not None for u in ups) else None
    return wlo2, wup2

def exclusion_ok(chart, B, status):
    Bd = dict(zip(NAMES[chart], B)); ch = "I" if chart == 1 else "II"
    if status == "LEMMA_A": return lemma_A(ch, Bd)
    if status == "LEMMA_B":
        if s_upper(ch, Bd) <= Fr(1, 10): return True
        _, wup2 = width_bounds(chart, B)            # s <= width (diam >= |P0P2| = 1 or |P0P1| = 1)
        return wup2 is not None and wup2 <= Fr(1, 100)
    if status == "OUT_OF_BAND":
        if s_lower(ch, Bd) > Fr(2, 5): return True
        wlo2, _ = width_bounds(chart, B)            # in-chart points have diam = 1, so s = width
        return wlo2 > Fr(4, 25)
    if status == "OUTSIDE_CHART": return check_exclusion(ch, Bd, "OUTSIDE_CHART")
    return False

def sha256(fn):
    import hashlib; return hashlib.sha256(open(fn, "rb").read()).hexdigest()

def check_log(fn, chart, n_cert_sample=0, seed=3, out=None):
    problems = []; leaves = {}; nrec = 0; subtrees = {}; summaries = []
    for ln, line in enumerate(open(fn), 1):
        try: r = json.loads(line)
        except Exception: problems.append((ln, "malformed")); continue
        if "subtree" in r:
            if set(r) != {"subtree", "done", "leaves", "sec"} or r["done"] is not True: problems.append((ln, "bad subtree marker"))
            if r["subtree"] in subtrees: problems.append((ln, "duplicate subtree marker"))
            subtrees[r["subtree"]] = r; continue
        if "summary" in r:
            summaries.append(r["summary"]); continue
        nrec += 1
        if r.get("status") not in STATUSES: problems.append((ln, f"unknown status {r.get('status')!r}")); continue
        if r.get("chart") != chart: problems.append((ln, "wrong chart")); continue
        p = r.get("path", "")
        if not p or any(c not in "abcdABCD" for c in p): problems.append((ln, "bad path")); continue
        if p in leaves: problems.append((ln, "duplicate path")); continue
        B = box_from_path(chart, p)
        lo = [Fr(x) for x in r["lo"]]; hi = [Fr(x) for x in r["hi"]]
        if [b[0] for b in B] != lo or [b[1] for b in B] != hi: problems.append((ln, "box mismatch")); continue
        leaves[p] = (r["status"], B, r)
    # tiling: prefix-free + complete binary tree (every proper prefix of a leaf has both children represented)
    paths = set(leaves)
    for p in paths:
        for k in range(1, len(p)):
            if p[:k] in paths: problems.append((p, "leaf is prefix of another leaf")); break
    prefixes = {p[:k] for p in paths for k in range(len(p))}
    from collections import defaultdict as _dd
    kids = _dd(set)
    for q in prefixes | paths:
        if q: kids[q[:-1]].add(q[-1])
    for pre in prefixes:
        L = kids.get(pre, set())
        if len(L) != 2 or len({c.lower() for c in L}) != 1 or not any(c.islower() for c in L):
            problems.append((pre, f"incomplete or inconsistent node children {sorted(L)}"))
    # metadata: 256 subtree markers whose leaf counts match the leaves under each 8-letter prefix; one summary matching the counts
    if len(subtrees) != 256: problems.append(("meta", f"{len(subtrees)} subtree markers (expected 256)"))
    from collections import Counter as _C
    per = _C(p[:8] for p in leaves)
    for k, r in subtrees.items():
        if per.get(k, 0) != r["leaves"]: problems.append(("meta", f"subtree {k}: marker says {r['leaves']} leaves, log has {per.get(k, 0)}"))
    if len(summaries) != 1: problems.append(("meta", f"{len(summaries)} summary records (expected 1)"))
    else:
        cnt = _C(st for st, _, _ in leaves.values())
        if summaries[0].get("leaves") != len(leaves) or {k: v for k, v in summaries[0].get("counts", {}).items() if v} != dict(cnt):
            problems.append(("meta", "summary counts do not match the leaves"))
    root_vol = math.prod(h-l for l, h in ROOTS[chart]); vol = sum(math.prod(h-l for l, h in B) for _, B, _ in leaves.values())
    if vol != root_vol: problems.append(("volume", f"{vol} != {root_vol}"))
    # exclusions
    nexc = 0
    for p, (st, B, r) in leaves.items():
        if st != "CERT":
            nexc += 1
            if not exclusion_ok(chart, B, st): problems.append((p, f"exclusion {st} ({r.get('rule')}) not confirmed"))
    # sampled certificates (independent 1-D implementation)
    certs = [p for p, (st, _, _) in leaves.items() if st == "CERT"]
    rng = random.Random(seed); sample = rng.sample(certs, min(n_cert_sample, len(certs))) if n_cert_sample else []
    ncert_ok = 0
    if sample:
        from xband import certify_box
        f = open(out, "w") if out else None
        for i, p in enumerate(sample):
            st, B, r = leaves[p]; Bd = dict(zip(NAMES[chart], B))
            for nmesh in (240, 600, 1200):            # finer meshes only change the test discretisation, not the rigour
                ok, E, R2 = certify_box("I" if chart == 1 else "II", Bd, n=nmesh)
                if ok: break
            ncert_ok += ok
            if f: f.write(json.dumps({"path": p, "ok": bool(ok), "mesh": nmesh, "E": float(E), "R2": float(R2.mid()) if R2 is not None else None}) + "\n"); f.flush()
            if not ok: problems.append((p, "certificate NOT reproduced by independent 1-D implementation"))
            if (i+1) % 50 == 0: print(f"  [{i+1}/{len(sample)}] independent certificates ok {ncert_ok}", flush=True)
    print(f"chart {chart}: records {nrec}, leaves {len(leaves)}, exclusions {nexc} checked, volume {'exact' if vol == root_vol else 'MISMATCH'}, "
          f"cert sample {ncert_ok}/{len(sample)}, problems {len(problems)} {problems[:8]}")
    return not problems

if __name__ == "__main__" and len(sys.argv) > 1:
    here = __file__.rsplit("/", 1)[0]
    print("PROVENANCE", json.dumps({"log_chart1": sha256(sys.argv[1]), "log_chart2": sha256(sys.argv[2]),
          "xbandcover.py": sha256(here + "/xbandcover.py"), "xband.py": sha256(here + "/xband.py"), "args": sys.argv[1:]}), flush=True)
    ok = True
    for fn, chart in ((sys.argv[1], 1), (sys.argv[2], 2)):
        ok &= check_log(fn, chart, n_cert_sample=int(sys.argv[3]) if len(sys.argv) > 3 else 0, out=f"xband_sample_chart{chart}.jsonl")
    print("XBANDCOVER", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
