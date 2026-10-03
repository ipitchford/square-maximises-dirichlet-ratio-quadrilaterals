"""Independent pieces of a Monte Carlo coverage check for darts (own code).  BD chart: B=(0,0), D=(1,0), A above, C below.
A = rhoA (cos bA, sin bA), C = rhoC (cos bC, -sin bC).  A dart with reflex vertex D: D strictly inside triangle ABC."""
import math, random
import numpy as np

def verts(bA, rA, bC, rC):
    return (0.0, 0.0), (1.0, 0.0), (rA*math.cos(bA), rA*math.sin(bA)), (rC*math.cos(bC), -rC*math.sin(bC))

def cross(o, a, b): return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])

from fractions import Fraction as _F
def _fin(*xs):
    for x in xs:
        if not math.isfinite(x): raise FloatingPointError("non-finite value in xdart")
def is_dart_direct(bA, rA, bC, rC):
    """D=(1,0) strictly inside triangle ABC: orientation tests in EXACT rational arithmetic on the float vertices."""
    _fin(bA, rA, bC, rC)
    B, D, A, C = verts(bA, rA, bC, rC); _fin(*A, *C)
    B, D, A, C = [(_F(x), _F(y)) for x, y in (B, D, A, C)]
    s1, s2, s3 = cross(A, B, D), cross(B, C, D), cross(C, A, D)
    return (s1 > 0 and s2 > 0 and s3 > 0) or (s1 < 0 and s2 < 0 and s3 < 0)

def f_lemma131(bA, rA, bC, rC):
    return math.sin(bA + bC)/(math.sin(bA)/rC + math.sin(bC)/rA)

def beta(bA, rA, bC, rC):
    """notch depth via the stable polar form 1 - beta = (sin bC/rA + sin bA/rC)/sin(thetaB) (Lemma 13.1; checked against the
    direct formula on 20,000 moderate configurations, max difference 2.6e-15)."""
    v = 1 - (math.sin(bC)/rA + math.sin(bA)/rC)/math.sin(bA + bC); _fin(v); return v

def angle(P, Q, R):   # angle at Q
    v1 = (P[0]-Q[0], P[1]-Q[1]); v2 = (R[0]-Q[0], R[1]-Q[1])
    return math.atan2(abs(v1[0]*v2[1]-v1[1]*v2[0]), v1[0]*v2[0]+v1[1]*v2[1])

def inradius(P, Q, R):
    a = math.dist(Q, R); b = math.dist(P, R); c = math.dist(P, Q); area = abs(cross(P, Q, R))/2
    return 2*area/(a+b+c)

def dist_point_segment(P, S0, S1):
    vx, vy = S1[0]-S0[0], S1[1]-S0[1]; t = max(0.0, min(1.0, ((P[0]-S0[0])*vx + (P[1]-S0[1])*vy)/(vx*vx+vy*vy)))
    return math.hypot(P[0]-S0[0]-t*vx, P[1]-S0[1]-t*vy)

def _h_tri(b, r, x):
    """vertical slice height at x in [0,1] of the triangle (B=(0,0), D=(1,0), X = r(cos b, sin b)) — stable polar form.
    Upper boundary over [0,1]: min of line BX (y = x tan b, only if X_x > 0 ... handled by cases) and line DX."""
    Xx = r*math.cos(b)
    if Xx >= 1:                      # apex beyond D: upper boundary is BX over [0,1]
        return x*math.tan(b)
    if Xx <= 0:                      # apex behind B: upper boundary is DX over [0,1]; slope = sin b/(cos b - 1/r)
        return (1 - x)*math.sin(b)/(1/r - math.cos(b))
    Xy = r*math.sin(b)               # apex above (0,1): tent
    return Xy*x/Xx if x <= Xx else Xy*(1 - x)/(1 - Xx)

def cstar(bA, rA, bC, rC):
    """max over x in [0,1] of h_A(x) + h_C(x) (piecewise linear; kinks at the apex abscissae)."""
    xs = [0.0, 1.0] + [v for v in (rA*math.cos(bA), rC*math.cos(bC)) if 0 < v < 1]
    v = max(_h_tri(bA, rA, x) + _h_tri(bC, rC, x) for x in xs); _fin(v); return v

L95 = [(40, 1318), (41, 418), (42, 230), (43, 155), (45, 92), (47.5, 61), (50, 46), (55, 31), (60, 25), (65, 21), (70, 19), (80, 17), (90, 16)]

def lemmas(bA, rA, bC, rC):
    """names of the analytic results that apply (own implementation of each hypothesis; small safety margins on floats)."""
    out = []
    B, D, A, C = verts(bA, rA, bC, rC)
    th = bA + bC
    if th <= 0.1772*(1-1e-9): out.append("L5.3 thin wedge at B")
    if math.tan(bA) + math.tan(bC) <= 0.1*(1-1e-9) and bA < math.pi/2 and bC < math.pi/2: out.append("L5.1")
    th_ = bA + bC
    angA = math.atan2(rC*math.sin(th_), rA - rC*math.cos(th_)); angC = math.atan2(rA*math.sin(th_), rC - rA*math.cos(th_))
    _fin(angA, angC)
    if max(angA, angC) <= math.atan(0.1)*(1-1e-9): out.append("L5.2 flat hull")
    if beta(bA, rA, bC, rC) <= 0.0339*(1-1e-9): out.append("P5.1 shallow")
    if dist_point_segment(A, B, D) <= 0.035*inradius(C, B, D)*(1-1e-9) or dist_point_segment(C, B, D) <= 0.035*inradius(A, B, D)*(1-1e-9):
        out.append("L3.1 collapse")
    if cstar(bA, rA, bC, rC) <= 0.05*(1-1e-9): out.append("L12.1 thin stem")
    for (b, r) in ((bA, rA), (bC, rC)):
        bd = math.degrees(b)
        for bk, Lk in L95:
            if bk <= bd <= 180 - bk and r >= Lk: out.append("L9.5 long steep arm"); break
    return out

def sample(rng, logr=(-6.0, 5.0)):
    """random dart, broad: angles uniform, radii log-uniform; returns (bA, rA, bC, rC) with bA <= bC after reflection."""
    while True:
        bA = rng.uniform(0, math.pi); bC = rng.uniform(0, math.pi)
        if bA + bC >= math.pi: continue
        rA = math.exp(rng.uniform(*logr)); rC = math.exp(rng.uniform(*logr))
        if not is_dart_direct(bA, rA, bC, rC): continue
        if bA > bC: bA, rA, bC, rC = bC, rC, bA, rA          # reflection in BD swaps A and C
        return bA, rA, bC, rC

if __name__ == "__main__":
    rng = random.Random(1); nd = 0; agree = 0; maxerr = 0.0
    for _ in range(20000):
        bA = rng.uniform(0.01, 3.1); bC = rng.uniform(0.01, 3.1); rA = math.exp(rng.uniform(-3, 3)); rC = math.exp(rng.uniform(-3, 3))
        direct = is_dart_direct(bA, rA, bC, rC)
        lem = (bA + bC < math.pi) and f_lemma131(bA, rA, bC, rC) > 1
        agree += (direct == lem); nd += direct
        if direct: maxerr = max(maxerr, abs((1 - 1/f_lemma131(bA, rA, bC, rC)) - beta(bA, rA, bC, rC)))
    print(f"Lemma 13.1 dart criterion agrees with direct test on {agree}/20000 configurations ({nd} darts); max |1-1/f - beta| = {maxerr:.2e}")
    # lemma predicates on known cases
    print("equilateral-hull shallow dart:", lemmas(math.radians(30), 1.0, math.radians(30), 1.0))
    print("thin wedge (5 deg):", lemmas(math.radians(2.5), 2.0, math.radians(2.5), 2.0))
    print("thin stem (A,C close to BD):", lemmas(math.radians(1), 1.5, math.radians(1.5), 0.02))
    S = [sample(rng) for _ in range(5)]; print("samples:", [tuple(round(v, 3) for v in s) for s in S])


# ---------------------------------------------------------------------------------------------- leaf index for dart logs
import json
class LeafIndex:
    """CERT leaves of one dart log as closed boxes p0 +- r (coordinates as logged); coordinates with r = 0 are fixed values."""
    def __init__(self, fn):
        P0, R, n_other = [], [], 0
        self.statuses = {}
        for l in open(fn):
            rec = json.loads(l)
            if rec.get("type") != "leaf": continue
            self.statuses[rec["status"]] = self.statuses.get(rec["status"], 0) + 1
            if rec["status"] != "cert": n_other += 1; continue
            P0.append(rec["p0"]); R.append(rec["r"])
        self.p0 = np.array(P0, dtype=float); self.r = np.array(R, dtype=float); self.fn = fn
    def contains(self, x, tol=0.0):
        """True if some CERT leaf contains the coordinate vector x (closed box; tol for float round-off of logged centres)"""
        if len(self.p0) == 0: return False
        x = np.asarray(x, dtype=float)
        if not np.all(np.isfinite(x)): raise FloatingPointError("non-finite coordinate in containment test")
        d = np.abs(self.p0 - x) <= self.r + tol
        return bool(np.any(np.all(d, axis=1)))

import os
DARTS = os.environ.get("XDART_DARTS") or os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "darts"))
DARTS = DARTS.rstrip("/") + "/"

if __name__ == "__main__" and "--test-index" in __import__("sys").argv:
    ix = LeafIndex(DARTS + "certs2/kf_spike.jsonl")
    print("kf_spike leaves by status:", ix.statuses, "cert boxes:", len(ix.p0))
    rng = np.random.default_rng(0); ok_in = 0
    for k in rng.choice(len(ix.p0), 200, replace=False):
        x = ix.p0[k] + ix.r[k]*rng.uniform(-1, 1, 4)
        ok_in += ix.contains(x)
    far = [ix.contains([5.0, 0.0, 3.0, 4.0]), ix.contains([0.0, 0.0, -1.0, 0.0])]
    print(f"points inside random cert leaves found: {ok_in}/200;  far-away points found: {far} (must be False)")


# ---------------------------------------------------------------------------------------------- coverage of one dart
import os
FAMDIR = os.environ.get("XDART_CERTS", "certs3")      # family logs: certs3 (pre-repair) or certs4 (re-certified)
assert FAMDIR in ("certs3", "certs4"), FAMDIR
LOGS = {"crit_collar": "certs/crit_collar.jsonl", "crit_spike_X1": "certs/crit_spike_X1.jsonl", "main_cbox": "certs/main_cbox.jsonl",
        "shallow_crit": "certs/shallow_crit.jsonl", "kf_collar": "certs2/kf_collar.jsonl", "kf_spike": "certs2/kf_spike.jsonl",
        "PP": FAMDIR + "/pp.merged.jsonl", "SP_col": FAMDIR + "/sp_col.merged.jsonl", "SP_spk": FAMDIR + "/sp_spk.merged.jsonl",
        "PF": FAMDIR + "/pf.merged.jsonl", "FP": FAMDIR + "/fp.merged.jsonl", "FF": FAMDIR + "/ff.merged.jsonl",
        "SF_col": FAMDIR + "/sf_col.merged.jsonl", "SF_spk": FAMDIR + "/sf_spk.merged.jsonl"}
def load_indices(names=None):
    ix = {}
    for k, rel in LOGS.items():
        if names and k not in names: continue
        if os.path.exists(DARTS + rel): ix[k] = LeafIndex(DARTS + rel)
    return ix

def covered(d, ix):
    """first reason that covers the dart d = (bA, rA, bC, rC) (bA <= bC), or None."""
    bA, rA, bC, rC = d
    lem = lemmas(bA, rA, bC, rC)
    if lem: return "lemma " + lem[0]
    ax, ay = rA*math.cos(bA), rA*math.sin(bA); cx, cy = rC*math.cos(bC), -rC*math.sin(bC)
    lA, lC = math.log(rA), math.log(rC)
    if ax > 0:
        t = math.tan(bA)
        if t <= 0.02:
            for k in ("crit_collar", "kf_collar"):
                if k in ix and ix[k].contains([0.02, 0.0, cx, cy]): return k
            if "SP_col" in ix and ix["SP_col"].contains([t, 0.0, bC, lC]): return "SP_col"
            if rC >= 20 and "SF_col" in ix and ix["SF_col"].contains([t, 0.0, bC, 0.0]): return "SF_col"
        if ax >= 1:
            for k in ("crit_spike_X1", "kf_spike"):
                if k in ix and ix[k].contains([t, 0.0, cx, cy]): return k
            if "SP_spk" in ix and ix["SP_spk"].contains([t, 0.0, bC, lC]): return "SP_spk"
            if rC >= 20 and "SF_spk" in ix and ix["SF_spk"].contains([t, 0.0, bC, 0.0]): return "SF_spk"
    for k in ("main_cbox", "shallow_crit"):
        if k in ix and ix[k].contains([ax, ay, cx, cy]): return k
    if "PP" in ix and ix["PP"].contains([bA, lA, bC, lC]): return "PP"
    if rC >= 20 and "PF" in ix and ix["PF"].contains([bA, lA, bC, 0.0]): return "PF"
    if rA >= 20 and "FP" in ix and ix["FP"].contains([bA, 0.0, bC, lC]): return "FP"
    if rA >= 20 and rC >= 20 and "FF" in ix and ix["FF"].contains([bA, 0.0, bC, 0.0]): return "FF"
    return None

if __name__ == "__main__" and "--test-cover" in __import__("sys").argv:
    ix = load_indices(); print("loaded:", {k: len(v.p0) for k, v in ix.items()})
    # known answers: (1) centre of a kf_spike cert leaf, with a = 1.5 (a >= 1): must be covered (by kf_spike or earlier reason)
    k = ix["kf_spike"]; p = k.p0[17]; t, cx, cy = p[0], p[2], p[3]
    bA = math.atan(t); rA = 1.5*math.hypot(1, t); bC = math.atan2(-cy, cx); rC = math.hypot(cx, cy)
    print("kf_spike centre dart:", is_dart_direct(bA, rA, bC, rC), covered((bA, rA, bC, rC), ix))
    # (2) same C, but a = 0.5 (< 1): the spike runs must NOT claim it
    rA2 = 0.5*math.hypot(1, t); print("same with a=0.5 (dart?)", is_dart_direct(bA, rA2, bC, rC), "->", covered((bA, rA2, bC, rC), ix))
    # (3) thin wedge
    print("thin wedge:", covered((math.radians(2), 3.0, math.radians(3), 3.0), ix))
    # (4) an FF dart: both arms long, 60/60 deg deep notch
    print("FF 60/60 deg, rho 50/50:", is_dart_direct(math.radians(60), 50, math.radians(60), 50), covered((math.radians(60), 50.0, math.radians(60), 50.0), ix))


# ---------------------------------------------------------------------------------------------- targeted boundary sampling
def jitter(rng, x, rel=0.03): return x*math.exp(rng.uniform(-rel, rel))
def boundary_samples(rng, n_each=2000):
    """darts concentrated near every boundary of the decomposition (Theorem 13.2) and of the lemma regions."""
    out = []
    def add(bA, rA, bC, rC):
        if 0 < bA < math.pi and 0 < bC < math.pi and bA + bC < math.pi and rA > 0 and rC > 0 and is_dart_direct(bA, rA, bC, rC):
            if bA > bC: bA, rA, bC, rC = bC, rC, bA, rA
            out.append(("", (bA, rA, bC, rC)))
    rnd = lambda a, b: rng.uniform(a, b); lrnd = lambda a, b: math.exp(rng.uniform(math.log(a), math.log(b)))
    kinds = {
      "t~0.02":   lambda: (math.atan(jitter(rng, 0.02)), lrnd(0.05, 60), rnd(0.02, 3.1), lrnd(0.01, 60)),
      "t~0.3":    lambda: (math.atan(jitter(rng, 0.3)), lrnd(0.5, 60), rnd(0.3, 2.8), lrnd(0.01, 60)),
      "ax~1":     lambda: (lambda t, ax: (math.atan(t), ax*math.hypot(1, t), rnd(0.02, 3.1), lrnd(0.01, 60)))(rnd(0.02, 0.3), jitter(rng, 1.0, 0.05)),
      "rhoA~e-2.5": lambda: (rnd(0.01, 1.6), jitter(rng, math.exp(-2.5), 0.1), rnd(0.01, 3.1), lrnd(0.01, 60)),
      "rhoC~e-4": lambda: (rnd(0.001, 0.3), lrnd(0.05, 60), rnd(0.01, 3.1), jitter(rng, math.exp(-4), 0.1)),
      "rho~20":   lambda: (rnd(0.01, 1.6), jitter(rng, 20.0, 0.05), rnd(0.01, 3.1), lrnd(0.01, 60)),
      "rhoC~20":  lambda: (rnd(0.01, 1.6), lrnd(0.05, 60), rnd(0.01, 3.1), jitter(rng, 20.0, 0.05)),
      "thetaB~0.1772": lambda: (lambda th, f: (th*f, lrnd(0.3, 60), th*(1-f), lrnd(0.3, 60)))(jitter(rng, 0.1772, 0.05), rnd(0.05, 0.5)),
      "thetaB~174.29deg": lambda: (lambda th, f: (th*f, lrnd(0.05, 60), th*(1-f), lrnd(0.05, 60)))(math.radians(jitter(rng, 174.29, 0.004)), rnd(0.01, 0.5)),
      "KF sector edge": lambda: (math.atan(rnd(0.0, 0.3)), lrnd(0.5, 60), math.radians(rng.choice([38, 100])*math.exp(rng.uniform(-0.03, 0.03))), rng.choice([0.5, 1.8])*math.exp(rng.uniform(-0.05, 0.05))),
    }
    for name, gen in kinds.items():
        k0 = len(out); tries = 0
        while len(out) - k0 < n_each and tries < 200*n_each:
            tries += 1; add(*gen())
            if len(out) > k0 and out[-1][0] == "": out[-1] = (name, out[-1][1])
    return out

def mc_report(ix, rng, n_uniform=50000, n_each=2000, logr=(-6.0, 5.0)):
    import collections
    res = collections.Counter(); unc = []
    for _ in range(n_uniform):
        d = sample(rng, logr); r = covered(d, ix); res["uniform:" + ("covered" if r else "UNCOVERED")] += 1
        if r is None: unc.append(("uniform", d))
    for name, d in boundary_samples(rng, n_each):
        r = covered(d, ix); res[name + ":" + ("covered" if r else "UNCOVERED")] += 1
        if r is None: unc.append((name, d))
    return res, unc


def sample_hull(rng):
    """shape-natural dart: random hull triangle (angles uniform on the simplex), reflex vertex D uniform inside, opposite vertex B
    chosen uniformly among the three hull vertices; mapped to the BD chart by a similarity; reflected so that bA <= bC."""
    while True:
        e = [rng.expovariate(1.0) for _ in range(3)]; s = sum(e); ang = [math.pi*x/s for x in e]
        if min(ang) < 1e-6: continue
        # triangle with angles ang at P0, P1, P2 (law of sines)
        P0 = (0.0, 0.0); P1 = (1.0, 0.0)
        a2 = ang[0]; L02 = math.sin(ang[1])/math.sin(ang[2])
        P2 = (L02*math.cos(a2), L02*math.sin(a2))
        T = [P0, P1, P2]
        w = [rng.expovariate(1.0) for _ in range(3)]; sw = sum(w)
        Dp = (sum(w[i]*T[i][0] for i in range(3))/sw, sum(w[i]*T[i][1] for i in range(3))/sw)
        k = rng.randrange(3); Bp = T[k]; others = [T[(k+1) % 3], T[(k+2) % 3]]
        # similarity: B -> 0, D -> (1,0)
        vx, vy = Dp[0]-Bp[0], Dp[1]-Bp[1]; L = math.hypot(vx, vy)
        if L < 1e-9: continue
        c, s_ = vx/L, vy/L
        def tf(P): x, y = P[0]-Bp[0], P[1]-Bp[1]; return ((c*x + s_*y)/L, (-s_*x + c*y)/L)
        X1, X2 = tf(others[0]), tf(others[1])
        if X1[1] * X2[1] >= 0: continue          # numerical degeneracy
        A, C = (X1, X2) if X1[1] > 0 else (X2, X1)
        bA = math.atan2(A[1], A[0]); rA = math.hypot(*A); bC = math.atan2(-C[1], C[0]); rC = math.hypot(*C)
        if not (0 < bA < math.pi and 0 < bC < math.pi) or not is_dart_direct(bA, rA, bC, rC): continue
        if bA > bC: bA, rA, bC, rC = bC, rC, bA, rA
        return bA, rA, bC, rC
