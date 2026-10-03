"""Independent coverage + discard checker for the global transport log (exact rational arithmetic; own code).
Checks: (1) the 24 root boxes are rebuilt exactly as documented (a-priori box cut at a,b,c = S and d = +-S);
(2) the bisection tree is complete and every leaf's logged box equals the reconstructed one (exact float equality),
so the leaves tile the a-priori box exactly; (3) every discard reason holds on the whole box, re-derived with own exact
polynomial bounds: E1 (cube), SHELL (inside the shell), E2 (c_i <= 0.035 c_{i+2}), E4 (some c_i <= 0), E3 (thin <= 1/10),
E3T (thin <= 2/5; deferred band), COMPOUND (recursive tree of these). Certified leaves are only counted here."""
import json, sys, math
from fractions import Fraction as Fr
from collections import defaultdict

def f_up(x):
    """smallest double >= exact rational x"""
    f = float(x)
    return f if Fr(f) >= x else math.nextafter(f, math.inf)
APLO = [0.0, 0.0, 0.0, -f_up(Fr(905, 1000))]; APHI = [f_up(Fr(905, 1000)), f_up(Fr(436, 100)), f_up(Fr(436, 100)), f_up(Fr(905, 1000))]
S = float(Fr(1, 5))           # Float64(1//5), round to nearest
DIM = "abcd"

MODE = "transport"
def f_down(x):
    f = float(x)
    return f if Fr(f) <= x else math.nextafter(f, -math.inf)
def roots():
    if MODE == "shell":
        # cover.jl --region 0,0.225,0,0.225,0,0.225,-0.225,0.225 ; cuts at +-Float64(49/1000, RoundDown); names R1..R24
        r = float("0.225"); c = f_down(Fr(49, 1000))
        lo0 = [0.0, 0.0, 0.0, -r]; hi0 = [r, r, r, r]
        def cut2(l, h, pts):
            xs = [l] + [x for x in pts if l < x < h] + [h]; return [(xs[k], xs[k+1]) for k in range(len(xs)-1)]
        ivs = [cut2(lo0[i], hi0[i], [c]) for i in range(3)] + [cut2(lo0[3], hi0[3], [-c, c])]
        out = []; k = 0
        for ia in ivs[0]:
            for ib in ivs[1]:
                for ic in ivs[2]:
                    for id_ in ivs[3]:
                        k += 1; out.append((f"R{k}", [ia[0], ib[0], ic[0], id_[0]], [ia[1], ib[1], ic[1], id_[1]]))
        return out
    def cut(l, h, pts):
        xs = [l] + [x for x in pts if l < x < h] + [h]; return [(xs[k], xs[k+1]) for k in range(len(xs)-1)]
    ivs = [cut(APLO[i], APHI[i], [S]) for i in range(3)] + [cut(APLO[3], APHI[3], [-S, S])]
    out = []; k = 0
    for ia in ivs[0]:
        for ib in ivs[1]:
            for ic in ivs[2]:
                for id_ in ivs[3]:
                    k += 1; out.append((f"T{k}", [ia[0], ib[0], ic[0], id_[0]], [ia[1], ib[1], ic[1], id_[1]]))
    return out

# ---- exact geometry: vertices affine in p; polynomials as dicts {exponent tuple: Fraction} -------------------------
H = Fr(1, 2)
VX = [([H, -H, 0, H], 0), ([-H, H, 0, H], 1), ([-H, -H, 0, -H], 1), ([H, H, 0, -H], 0)]
VY = [([-H, 0, H, H], 0), ([-H, 0, -H, -H], 0), ([H, 0, H, -H], 1), ([H, 0, -H, H], 1)]
def lin(coefs, const): 
    p = {(0, 0, 0, 0): Fr(const)}
    for i, c in enumerate(coefs):
        if c: e = [0]*4; e[i] = 1; p[tuple(e)] = p.get(tuple(e), 0) + Fr(c)
    return p
def add(a, b, s=1):
    out = dict(a)
    for k, v in b.items(): out[k] = out.get(k, 0) + s*v
    return {k: v for k, v in out.items() if v != 0}
def mul(a, b):
    out = defaultdict(Fr)
    for ka, va in a.items():
        for kb, vb in b.items(): out[tuple(x+y for x, y in zip(ka, kb))] += va*vb
    return {k: v for k, v in out.items() if v != 0}
def scal(a, s): return {k: v*s for k, v in a.items()}
VXP = [lin(*VX[i]) for i in range(4)]; VYP = [lin(*VY[i]) for i in range(4)]
def cross(ax, ay, bx, by): return add(mul(ax, by), mul(ay, bx), -1)
def cpoly(i):
    j, k = (i+1) % 4, (i+2) % 4
    return cross(add(VXP[j], VXP[i], -1), add(VYP[j], VYP[i], -1), add(VXP[k], VXP[i], -1), add(VYP[k], VYP[i], -1))
CP = [cpoly(i) for i in range(4)]

def shift_bound(poly, lo, hi, upper=True):
    """rigorous upper (or lower) bound over the box: expand around the centre (exact), |delta^g| <= r^g."""
    c = [(Fr(l)+Fr(h))/2 for l, h in zip(lo, hi)]; r = [(Fr(h)-Fr(l))/2 for l, h in zip(lo, hi)]
    # substitute p = c + delta and expand
    out = defaultdict(Fr)
    for e, v in poly.items():
        terms = {(0, 0, 0, 0): v}
        for i in range(4):
            for _ in range(e[i]):
                nt = defaultdict(Fr)
                for k, w in terms.items():
                    nt[k] += w*c[i]; kk = list(k); kk[i] += 1; nt[tuple(kk)] += w
                terms = nt
        for k, w in terms.items(): out[k] += w
    tot = out.get((0, 0, 0, 0), Fr(0)); sgn = 1 if upper else -1
    for k, w in out.items():
        if sum(k) == 0 or w == 0: continue
        tot += sgn*abs(w)*math.prod(r[i]**k[i] for i in range(4))
    return tot

def nonconvex_necessary(lo, hi):
    """every convex Q_p (all c_i > 0) has |a|,|d| < 1 and |b|,|c| < 1 + |a| (paper section 2.4); True if the box violates
    this at every point, i.e. the box contains no convex parameter."""
    lo = [Fr(x) for x in lo]; hi = [Fr(x) for x in hi]
    amax = max(abs(lo[0]), abs(hi[0])); dist = lambda l, h: Fr(0) if l <= 0 <= h else min(abs(l), abs(h))
    if dist(lo[0], hi[0]) >= 1 or dist(lo[3], hi[3]) >= 1: return True
    return dist(lo[1], hi[1]) >= 1 + amax or dist(lo[2], hi[2]) >= 1 + amax

def pred(name, lo, hi):
    lo_f = [Fr(x) for x in lo]; hi_f = [Fr(x) for x in hi]
    if name == "N": return nonconvex_necessary(lo, hi)
    if name == "E1":
        return all(-Fr(49, 1000) <= l and h <= Fr(49, 1000) for l, h in zip(lo_f, hi_f))
    if name == "SHELL":
        return all(Fr(0) <= lo_f[i] and hi_f[i] <= Fr(225, 1000) for i in range(3)) and -Fr(225, 1000) <= lo_f[3] and hi_f[3] <= Fr(225, 1000)
    if name == "E4":
        return any(shift_bound(CP[i], lo, hi) <= 0 for i in range(4))
    if name == "E2":
        return any(shift_bound(add(CP[i], CP[(i+2) % 4], -Fr(35, 1000)), lo, hi) <= 0 for i in range(4))
    if name in ("E3", "E3T"):
        s2 = Fr(1, 100) if name == "E3" else Fr(4, 25)
        E = [(add(VXP[(i+1) % 4], VXP[i], -1), add(VYP[(i+1) % 4], VYP[i], -1)) for i in range(4)]
        for k in range(4):
            for l in range(k+1, 4):
                dx = add(VXP[k], VXP[l], -1); dy = add(VYP[k], VYP[l], -1); d2 = add(mul(dx, dx), mul(dy, dy))
                for i in range(4):
                    e2 = add(mul(E[i][0], E[i][0]), mul(E[i][1], E[i][1])); rhs = scal(mul(e2, d2), s2)
                    ok = True
                    for j in ((i+2) % 4, (i+3) % 4):
                        a = cross(E[i][0], E[i][1], add(VXP[j], VXP[i], -1), add(VYP[j], VYP[i], -1))
                        if shift_bound(add(mul(a, a), rhs, -1), lo, hi) > 0: ok = False; break
                    if ok: return True
        return False
    raise ValueError(name)

def any_exclusion(lo, hi, depth=0, allow=("E1", "SHELL", "N", "E4", "E2", "E3", "E3T"), maxdepth=8):
    """own refinement: True if every point of the box satisfies one allowed exclusion (bisect widest side)."""
    for name in allow:
        if pred(name, lo, hi): return True
    if depth >= maxdepth: return False
    w = [h - l for l, h in zip(lo, hi)]; d = w.index(max(w)); mid = (lo[d] + hi[d]) / 2
    hi1 = list(hi); hi1[d] = mid; lo2 = list(lo); lo2[d] = mid
    return any_exclusion(lo, hi1, depth+1, allow, maxdepth) and any_exclusion(lo2, hi, depth+1, allow, maxdepth)

def compound_ok(tree, lo, hi):
    """tree grammar: leaf = condition name; node = '<dim>(<left>,<right>)' with dim a letter or index"""
    tree = tree.strip()
    if "(" not in tree: return pred(tree, lo, hi) or pred("N", lo, hi) or any_exclusion(lo, hi, allow=(tree, "N", "E4"))
    dimc = tree[:tree.index("(")]; inner = tree[tree.index("(")+1:-1]
    depth = 0
    for k, ch in enumerate(inner):
        if ch == "(": depth += 1
        elif ch == ")": depth -= 1
        elif ch == "," and depth == 0: left, right = inner[:k], inner[k+1:]; break
    d = DIM.index(dimc) if dimc in DIM else int(dimc) - 1
    mid = (lo[d] + hi[d]) / 2
    hi1 = list(hi); hi1[d] = mid; lo2 = list(lo); lo2[d] = mid
    return compound_ok(left, lo, hi1) and compound_ok(right, lo2, hi)

def main(fn, check_discards=True, max_discards=None):
    recs = {}; bad = []
    KNOWN = {"split", "certified", "E1", "E2", "E3", "E3T", "E4", "SHELL", "COMPOUND"}
    for ln, line in enumerate(open(fn), 1):
        try:
            r = json.loads(line)
        except Exception:
            bad.append((f"line {ln}", "malformed JSON")); continue
        if not isinstance(r, dict) or "path" not in r or "status" not in r:
            bad.append((f"line {ln}", "missing path/status")); continue
        if r["status"] not in KNOWN: bad.append((r["path"], f"unknown status {r['status']!r}"))
        if r["path"] in recs: bad.append((r["path"], "duplicate path"))
        if r["status"] == "split" and r.get("dim") not in (1, 2, 3, 4): bad.append((r["path"], "split without valid dim"))
        if r["status"] == "COMPOUND" and not isinstance(r.get("tree"), str): bad.append((r["path"], "COMPOUND without tree"))
        if r["status"] != "split" and ("lo" not in r or "hi" not in r): bad.append((r["path"], "leaf without box"))
        recs[r["path"]] = r
    reached = set()
    stats = defaultdict(int); vol = defaultdict(Fr); pending = []; nd = 0
    stack = [(name, lo, hi) for name, lo, hi in roots()]
    while stack:
        path, lo, hi = stack.pop()
        r = recs.get(path)
        if r is None: pending.append(path); continue
        reached.add(path)
        st = r["status"]
        if st == "split":
            d = r["dim"] - 1; mid = (lo[d] + hi[d]) / 2
            assert lo[d] < mid < hi[d]
            hi1 = list(hi); hi1[d] = mid; lo2 = list(lo); lo2[d] = mid
            stack.append((f"{path}.{DIM[d]}0", list(lo), hi1)); stack.append((f"{path}.{DIM[d]}1", lo2, list(hi)))
            continue
        if r["lo"] != lo or r["hi"] != hi: bad.append((path, "box mismatch")); continue
        if st not in ("certified", "E1", "E2", "E3", "E3T", "E4", "SHELL", "COMPOUND"):
            bad.append((path, f"leaf with non-leaf/unknown status {st!r}")); continue
        stats[st] += 1; vol[st] += math.prod(Fr(h)-Fr(l) for l, h in zip(lo, hi))
        if check_discards and st in ("E1", "SHELL", "E2", "E3", "E3T", "E4", "COMPOUND") and (max_discards is None or nd < max_discards):
            nd += 1
            ok = compound_ok(r["tree"], lo, hi) if st == "COMPOUND" else (pred(st, lo, hi) or pred("N", lo, hi) or any_exclusion(lo, hi, allow=(st, "N", "E4")))
            if not ok: bad.append((path, f"discard {st} not verified"))
    rts = roots(); total = sum(math.prod(Fr(h) - Fr(l) for l, h in zip(lo, hi)) for _, lo, hi in rts)
    print("leaf counts:", dict(stats)); print("pending (no record):", len(pending), pending[:5])
    print("volume fractions:", {k: float(v/total) for k, v in vol.items()}, " sum =", float(sum(vol.values())/total))
    orphans = [p for p in recs if p not in reached]
    for p in orphans[:1000]: bad.append((p, "orphan record (not in any tree)"))
    if check_discards and max_discards is not None: bad.append(("(run)", f"discards only partially checked (max {max_discards})"))
    print("orphan records:", len(orphans))
    print("problems:", len(bad), bad[:10])
    ok = not bad and not pending
    print("XCOVER", "PASS" if ok else "FAIL")
    return stats, bad, pending, ok

def _sha256(fn):
    import hashlib; return hashlib.sha256(open(fn, "rb").read()).hexdigest()

if __name__ == "__main__":
    if "--shell" in sys.argv: MODE = "shell"; sys.argv.remove("--shell")
    print("PROVENANCE", json.dumps({"log": _sha256(sys.argv[1]), "xcover.py": _sha256(__file__), "mode": MODE, "args": sys.argv[1:]}), flush=True)
    _, _, _, ok = main(sys.argv[1], check_discards=(len(sys.argv) < 3 or sys.argv[2] != "nodiscard"),
                       max_discards=int(sys.argv[3]) if len(sys.argv) > 3 else None)
    sys.exit(0 if ok else 1)
