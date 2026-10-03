"""Stage 4: certified cover of the punctured ball {0 < |p|_2 <= rho0} by face-direction boxes, using the multivariate
Schur-complement Taylor models (schur_mtm.py).  Parametrisation p = s*E with E on a face of the cube [-1,1]^4
(E_i = +-1 for the face index i), s = |p|_inf > 0.  For a box of E's we certify s in (0, s_max], s_max <= rho0/max|E|,
so |p|_2 <= rho0 (validity of the Taylor models).  Per (box, s-interval) we pick beta = the largest precomputed level
<= the bounding-rectangle lower bound for lambda_2 (or 36 <= 2 pi j01^2/|Q_p|, Krahn-Szego), and verify
detS/s^2 < 0 or trS/s^2 < 0 by a second-order centred form in E (exact point value and gradient at the box centre,
interval Hessian over the box) plus the remainder r s^(m-1) |E|^(m+1)."""
import sys, time, json, pickle, math
from fractions import Fraction as Fr
from flint import arb
from mtm import deg
from schur_mtm import deser
from localcert import fr2arb

PI = arb.pi()
def iv(lo, hi): lo, hi = Fr(lo), Fr(hi); return fr2arb((lo+hi)/2) + fr2arb((hi-lo)/2)*arb("[+/- 1]")
def sym(h): return fr2arb(Fr(h)) * arb("[+/- 1]")

class Poly:
    """f(E,s) = sum_{|g|>=2} c_g E^g s^{|g|-2}  (+ remainder r |p|^{m+1}/s^2)."""
    def __init__(self, t):
        self.m = t.m; self.r = t.r
        self.terms = [(g, v) for g, v in t.c.items()]
        assert all(deg(g) >= 2 for g, _ in self.terms)
    def parts(self, c, X):
        """per degree k: value and gradient at point c, Hessian over interval box X (all arb)."""
        m = self.m
        cp = [[arb(1)] + [None]*m for _ in range(4)]; Xp = [[arb(1)] + [None]*m for _ in range(4)]
        for i in range(4):
            for e in range(1, m+1): cp[i][e] = cp[i][e-1]*c[i]; Xp[i][e] = Xp[i][e-1]*X[i]
        V = [arb(0) for _ in range(m+1)]; G = [[arb(0)]*4 for _ in range(m+1)]; H = [[[arb(0)]*4 for _ in range(4)] for _ in range(m+1)]
        for g, v in self.terms:
            k = deg(g)
            V[k] += v * cp[0][g[0]]*cp[1][g[1]]*cp[2][g[2]]*cp[3][g[3]]
            for i in range(4):
                if g[i] == 0: continue
                gi = list(g); gi[i] -= 1
                G[k][i] = G[k][i] + v * g[i] * cp[0][gi[0]]*cp[1][gi[1]]*cp[2][gi[2]]*cp[3][gi[3]]
                for j in range(i, 4):
                    gj = list(gi)
                    if gj[j] == 0: continue
                    coef = g[i] * gj[j]; gj[j] -= 1
                    term = v * coef * Xp[0][gj[0]]*Xp[1][gj[1]]*Xp[2][gj[2]]*Xp[3][gj[3]]
                    H[k][i][j] = H[k][i][j] + term
        return V, G, H
    def enclose(self, parts, h, slo, shi, Emax):
        V, G, H = parts; m = self.m
        Sp = [iv(slo**(k-2), shi**(k-2)) if k >= 2 else None for k in range(m+1)]
        val = arb(0)
        for k in range(2, m+1):
            if V[k] == 0 and all(x == 0 for x in G[k]): continue
            lin = V[k] + sum((G[k][i]*sym(h[i]) for i in range(4) if h[i] != 0), arb(0))
            quad = arb(0)
            for i in range(4):
                for j in range(i, 4):
                    if h[i] == 0 or h[j] == 0: continue
                    if i == j: quad += H[k][i][i] * iv(0, Fr(h[i])**2) / 2      # 1/2 H_ii d_i^2, d_i^2 in [0,h^2]
                    else: quad += H[k][i][j] * sym(Fr(h[i])*Fr(h[j]))            # H_ij d_i d_j (i<j counted once: 2*1/2)
            val += (lin + quad) * Sp[k]
        rem = self.r * fr2arb(shi)**(m-1) * Emax**(m+1)
        return val + rem*arb("[+/- 1]")

def beta_rect_lower(E_lo, E_hi, slo, shi):
    """lower bound of pi^2 min(4/A^2+1/B^2, 1/A^2+4/B^2), A = 1 + s(-Ea + max(|Eb|,|Ed|)), B = 1 + s(Ea + max(|Ec|,|Ed|))."""
    amax = lambda lo, hi: max(abs(lo), abs(hi))
    kA = -E_lo[0] + max(amax(E_lo[1], E_hi[1]), amax(E_lo[3], E_hi[3]))
    kB = E_hi[0] + max(amax(E_lo[2], E_hi[2]), amax(E_lo[3], E_hi[3]))
    Aup = 1 + max(slo*kA, shi*kA); Bup = 1 + max(slo*kB, shi*kB)
    assert Aup > 0 and Bup > 0
    A2 = fr2arb(Aup)**2; B2 = fr2arb(Bup)**2
    v1 = PI**2*(4/A2 + 1/B2); v2 = PI**2*(1/A2 + 4/B2)
    return min(arb_lower_fraction(v1), arb_lower_fraction(v2))       # exact rational lower endpoints (no rounding; referee audit 2026-10-02)

def arb_lower_fraction(x):
    """the exact rational value of a rigorous lower endpoint of the ball x"""
    m, e = x.lower().man_exp()
    return Fr(int(m)) * Fr(2)**int(e)

def arb_upper_fraction(x):
    """the exact rational value of a rigorous upper endpoint of the ball x"""
    m, e = x.upper().man_exp()
    return Fr(int(m)) * Fr(2)**int(e)

def run(schur_fn, k_init=4, max_depth=18, nsplit=6, log_every=2000):
    S = pickle.load(open(schur_fn, 'rb')); meta = S['meta']; rho0 = Fr(meta['rho0']); rho = fr2arb(rho0)
    levels = sorted(Fr(b) for b in S['res'])
    polys = {b: (Poly(deser(S['res'][str(b)]['detS'], rho)), Poly(deser(S['res'][str(b)]['trS'], rho))) for b in levels}
    t0 = time.time(); queue = []
    for i in range(4):
        for sgn in (1, -1):
            for idx in range(k_init**3):
                a, b, c = idx // k_init**2, (idx // k_init) % k_init, idx % k_init
                lo = [Fr(-1) + Fr(2*x, k_init) for x in (a, b, c)]; hi = [l + Fr(2, k_init) for l in lo]
                queue.append(((i, sgn), lo, hi, 0))
    nok = nfail = neval = 0; maxd = 0; min_smax = None; min_ratio = None; cert = []; crit_count = {'det': 0, 'tr': 0}
    while queue:
        face, lo3, hi3, depth = queue.pop(); neval += 1
        i, sgn = face; free = [j for j in range(4) if j != i]
        E_lo = [None]*4; E_hi = [None]*4; E_lo[i] = E_hi[i] = Fr(sgn)
        for j, l, h in zip(free, lo3, hi3): E_lo[j], E_hi[j] = l, h
        c = [fr2arb((E_lo[j]+E_hi[j])/2) for j in range(4)]; hw = [(E_hi[j]-E_lo[j])/2 for j in range(4)]
        X = [fr2arb((E_lo[j]+E_hi[j])/2) + fr2arb(hw[j])*arb("[+/- 1]") for j in range(4)]
        Emax2 = sum(max(E_lo[j]**2, E_hi[j]**2) for j in range(4)); Emin2 = sum((0 if E_lo[j] <= 0 <= E_hi[j] else min(E_lo[j]**2, E_hi[j]**2)) for j in range(4))
        Emax = fr2arb(Emax2).sqrt()
        smax = rho0 / 2                                   # exact: the certified set is the cube 0 < |p|_inf <= rho0/2
        assert smax**2 * Emax2 <= rho0**2, "s_max violates the Taylor-model ball"   # exact rational check (|E|^2 <= 4)
        cache = {}; ok_box = True; used = []
        for q in range(nsplit):
            slo, shi = smax*q/nsplit, smax*(q+1)/nsplit
            br = beta_rect_lower(E_lo, E_hi, slo, shi); lv = [b for b in levels if b <= max(br, 36)]
            assert lv, "no admissible beta level"; beta = lv[-1]
            if beta not in cache: cache[beta] = tuple(p.parts(c, X) for p in polys[beta])
            pdet, ptr = polys[beta]
            ed = pdet.enclose(cache[beta][0], hw, slo, shi, Emax)
            if ed < 0: used.append((str(beta), 'det', str(arb_upper_fraction(ed)), float(ed.upper()))); continue
            et = ptr.enclose(cache[beta][1], hw, slo, shi, Emax)
            if et < 0: used.append((str(beta), 'tr', str(arb_upper_fraction(et)), float(et.upper()))); continue
            ok_box = False; break
        if ok_box:
            nok += 1; maxd = max(maxd, depth)
            for u in used: crit_count[u[1]] += 1
            min_smax = smax if min_smax is None else min(min_smax, smax)
            ratio = math.sqrt(float(Emin2)/float(Emax2)); min_ratio = ratio if min_ratio is None else min(min_ratio, ratio)
            cert.append(dict(face=list(face), lo=[str(x) for x in lo3], hi=[str(x) for x in hi3], smax=str(smax), used=used))
        elif depth >= max_depth:
            nfail += 1; print(f"FAIL face {face} box {[float(x) for x in lo3]}..{[float(x) for x in hi3]} at s-interval {q}", flush=True)
        else:
            w = [h - l for l, h in zip(lo3, hi3)]; j = w.index(max(w)); mid = (lo3[j] + hi3[j])/2
            l1, h1 = list(lo3), list(hi3); h1[j] = mid; l2, h2 = list(lo3), list(hi3); l2[j] = mid
            queue += [(face, l1, h1, depth+1), (face, l2, h2, depth+1)]
        if neval % log_every == 0: print(f"  evals {neval} ok {nok} queue {len(queue)} maxdepth {maxd} {time.time()-t0:.0f}s", flush=True)
    summary = dict(schur=schur_fn, meta=meta, levels=[str(b) for b in levels], certified_boxes=nok, failures=nfail, evaluations=neval,
                   max_depth=maxd, seconds=round(time.time()-t0, 1), min_smax=str(min_smax), cube_halfwidth=str(min_smax),
                   crit_count=crit_count)
    print("SUMMARY", json.dumps(summary), flush=True)
    json.dump(dict(summary=summary, boxes=cert), open(schur_fn.replace('.pkl', '_cover2.json'), 'w'))
    return nfail == 0

if __name__ == "__main__":
    ok = run(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 4, int(sys.argv[3]) if len(sys.argv) > 3 else 18, int(sys.argv[4]) if len(sys.argv) > 4 else 6)
    sys.exit(0 if ok else 1)
