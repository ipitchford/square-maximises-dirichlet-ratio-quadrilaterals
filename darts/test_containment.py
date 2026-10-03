"""Containment tests of the bound functions used by the certificates (adversarial review 2026-10-02, gate 2).
Every test compares a production bound against an EXACT rational value (fractions.Fraction) or a 300-bit ball-arithmetic
reference (python-flint at prec 300, an independent evaluation) and requires strict containment: lower bound <= every point of
the reference, upper bound >= every point of the reference; no tolerance anywhere.  Inputs: random values, edge cases
(balls with zero midpoint, huge and tiny radii, exact zeros), and the reviewer's witnesses
  pp.merged.jsonl:66327 (Schur), :189234 (strain factors), :80526 (FARL tail / cut), sp_spk.merged.jsonl:7 (atan endpoint),
  fp.merged.jsonl:321 (TM range with a quadratic coefficient containing 0), :325 (FARL tail).
Sections: helpers (fdown/fup/amax/amin/enclose/exact_box/cuts/dec, dyadic snapping), TM range, Schur upper bound, strain
factors, atan / cot endpoints, FARL tail and cut validity, BO2 projection extrema, spike tail.
Usage: test_containment.py [--quick]   (exit 0 iff every check passes)"""
import sys, os, json, math, random, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from flint import arb, ctx
import nanguard as ng
from nanguard import fdown, fup, amax, amin, enclose, exact_box, cuts, dec, InexactBox
PROD = 106
QUICK = '--quick' in sys.argv
STATS = {}
def check(section, cond, msg):
    s = STATS.setdefault(section, [0, 0])
    if cond: s[0] += 1
    else:
        s[1] += 1
        if s[1] <= 5: print(f'FAIL [{section}] {msg}', flush=True)

def F(x):
    """exact Fraction of a float or of an EXACT arb point"""
    if isinstance(x, arb):
        if x.rad() != 0: raise ValueError('not an exact arb point')
        m, e = x.mid().man_exp(); m, e = int(m), int(e)   # RIGOR: exact: the exact midpoint of a radius-0 ball as mantissa and exponent
        return Fr(m)*Fr(2)**e   # RIGOR: exact: Fraction ** int
    return Fr(float(x))   # RIGOR: exact: float() of floats (test inputs)
def lo_of(x): return F(x.lower())
def hi_of(x): return F(x.upper())
class prec:
    def __init__(self, p): self.p = p
    def __enter__(self): self.old = ctx.prec; ctx.prec = self.p
    def __exit__(self, *a): ctx.prec = self.old
def ref(fn):
    """evaluate fn() at 300 bits (independent, rigorous reference ball)"""
    with prec(300): return fn()

rng = random.Random(20261002)
def rand_float(scale=None):
    k = rng.choice([0, 1, 2, 3, 4])
    if k == 0: return 0.0
    if k == 1: return rng.uniform(-1, 1)
    if k == 2: return rng.uniform(-1, 1)*10**rng.uniform(-30, 30)
    if k == 3: return rng.choice([1.0, -1.0])*2.0**rng.randint(-60, 60)
    return rng.uniform(-1e3, 1e3)
def rand_ball():
    m = rand_float(); k = rng.choice([0, 1, 2, 3])
    r = [0.0, 10**rng.uniform(-40, -10), 10**rng.uniform(-5, 5), abs(m)*10**rng.uniform(-17, -1)][k]
    return arb(m, r)

# ------------------------------------------------------------------------------------------------ helpers
def test_helpers():
    N = 3000 if not QUICK else 600
    for _ in range(N):
        x = rand_ball()
        try:
            d, u = fdown(x), fup(x)
        except ng.NaNError:
            check('helpers', False, f'fdown/fup raised on finite {x}'); continue
        check('helpers', F(d) <= lo_of(x) and F(u) >= hi_of(x), f'fdown/fup containment {x}: {d} {u}')
        # tightness: the next float outward would no longer be <= / >= (directed rounding of the exact endpoint)
        check('helpers', F(math.nextafter(d, math.inf)) > lo_of(x) or d == lo_of(x), f'fdown not the largest float <= {x}')   # RIGOR: exact: test
        check('helpers', F(math.nextafter(u, -math.inf)) < hi_of(x) or u == hi_of(x), f'fup not the smallest float >= {x}')   # RIGOR: exact: test
        y = rand_ball()
        A, B = amax(x, y), amin(x, y)
        check('helpers', F(A) >= hi_of(x) and F(A) >= hi_of(y) and F(B) <= lo_of(x) and F(B) <= lo_of(y), 'amax/amin')
        a, b = sorted([rand_float(), rand_float()])
        if math.isfinite(a) and math.isfinite(b) and abs(a) < 1e300 and abs(b) < 1e300:
            c, h = enclose(a, b)
            check('helpers', F(c) - F(h) <= F(a) and F(c) + F(h) >= F(b) and F(c - h) == F(c) - F(h) and F(c + h) == F(c) + F(h),
                  f'enclose {a} {b} -> {c} {h}')
            if a < b:
                n = rng.randint(1, 9); xs = cuts(a, b, n)
                check('helpers', xs[0] == a and xs[-1] == b and len(xs) == n + 1, 'cuts ends')
    for bad in ([0.025], [0.004999999999999999]), ([0.1], [0.3]):
        try: exact_box(*bad); check('helpers', False, f'exact_box accepted {bad}')
        except InexactBox: check('helpers', True, '')
    for s in ('0.0339', '1.8', '0.3', '174.29', '0.1772', '0.05', '38', '100', '0.5', '0.02'):
        d = dec(s); q = Fr(s)
        check('helpers', lo_of(d) <= q <= hi_of(d), f'dec({s}) does not contain {s}')
    from grid import dyadic_root
    for _ in range(N):
        p0 = [rng.uniform(-5, 5) for _ in range(4)]; r = [rng.choice([0.0, rng.uniform(1e-9, 1)]) for _ in range(4)]
        sp, sr = dyadic_root(p0, r)
        for i in range(4):
            check('helpers', F(sp[i]) - F(sr[i]) <= F(p0[i]) - F(r[i]) and F(sp[i]) + F(sr[i]) >= F(p0[i]) + F(r[i]), 'dyadic_root containment')

# ------------------------------------------------------------------------------------------------ TM range
from tm import TM, MON, NC, NV, IDX
def tm_eval_exact(coefs, t):
    s = Fr(0)
    for k, m in enumerate(MON):
        v = coefs[k]
        for i in m: v = v*t[i]
        s += v
    return s
def tm_samples(T, nsamp):
    """exact values of polynomials of the model at sample points: coefficients at ball ends / mid, remainder at +-R, 0"""
    out = []
    for _ in range(nsamp):
        coefs = []
        for c in T.c:
            if c.rad() == 0: coefs.append(F(c.mid()))   # RIGOR: exact: the exact midpoint of a radius-0 ball
            else: coefs.append(rng.choice([lo_of(c), hi_of(c), F(c.mid())]))   # RIGOR: exact: exact ball endpoints / midpoint (all points of the ball)
        t = [Fr(rng.choice([-1, 1, 0, Fr(rng.randint(-1000, 1000), 1000)])) for _ in range(NV)]
        # also the 1-D vertices of the quadratic parts (where the old code formed a NaN candidate)
        if rng.random() < 0.3:
            i = rng.randrange(NV); a = coefs[IDX[(i, i)]]; b = coefs[IDX[(i,)]]
            if a != 0 and abs(b) <= 2*abs(a): t[i] = -b/(2*a)
        e = rng.choice([Fr(0), hi_of(T.R), -hi_of(T.R)])
        out.append(tm_eval_exact(coefs, t) + e)
    return out
def rand_tm(zero_quad=False):
    c = []
    for k in range(NC):
        r = rng.random()
        if r < 0.25: c.append(arb(0))
        elif r < 0.45: c.append(arb(0, 10**rng.uniform(-20, 1)))              # zero-midpoint ball
        elif r < 0.6: c.append(arb(rng.uniform(-3, 3), 10**rng.uniform(-30, -5)))
        elif r < 0.7: c.append(arb(rng.uniform(-1e6, 1e6), 10**rng.uniform(0, 6)))   # huge radius
        else: c.append(arb(rng.uniform(-3, 3)))
    if zero_quad:
        for i in range(NV): c[IDX[(i, i)]] = arb(rng.uniform(-1e-3, 1e-3), 2e-3)
    return TM(c, arb(abs(rng.uniform(0, 1))*10**rng.uniform(-30, 0)))
def test_tm():
    N = 400 if not QUICK else 80
    for j in range(N):
        T = rand_tm(zero_quad=(j % 2 == 0))
        try: lo, hi = T.range()
        except ng.NaNError: check('tm', False, 'range raised NaNError on a finite model'); continue
        check('tm', lo.is_finite() and hi.is_finite(), 'non-finite range')
        L, H = lo_of(lo), hi_of(hi)
        for v in tm_samples(T, 40): check('tm', L <= v <= H, f'TM range does not contain a sample: {float(v)} not in [{float(L)}, {float(H)}]')   # RIGOR: info: message

# ------------------------------------------------------------------------------------------------ Schur upper bound
from schur import schur_upper
def small_tm(scale):
    """a Taylor model with coefficients |mid| <= scale, including zero-midpoint balls and huge / tiny relative radii"""
    c = []
    for k in range(NC):
        r = rng.random()
        if r < 0.3: c.append(arb(0))
        elif r < 0.5: c.append(arb(0, scale*10**rng.uniform(-25, -1)))                 # zero-midpoint ball
        elif r < 0.6: c.append(arb(scale*rng.uniform(-1, 1), scale*10**rng.uniform(-3, -1)))   # large relative radius
        else: c.append(arb(scale*rng.uniform(-1, 1), scale*10**rng.uniform(-30, -12)))       # tiny radius
    return TM(c, arb(scale*10**rng.uniform(-20, -2)))
def test_schur():
    N = 400 if not QUICK else 80
    nok = 0
    for j in range(N):
        A = [[None]*2 for _ in range(2)]
        A[0][0] = small_tm(rng.choice([0.5, 5.0])) - rng.uniform(10, 60)
        A[1][1] = small_tm(rng.choice([0.5, 5.0])) - rng.uniform(10, 60)
        A[0][1] = small_tm(rng.choice([0.1, 3.0, 20.0])); A[1][0] = A[0][1]
        ok, schur, why = schur_upper(A)
        if schur is None: continue
        n11u = hi_of(A[0][0].upper())
        lo, hi = A[0][1].range(); m = max(abs(lo_of(lo)), abs(hi_of(hi)))   # RIGOR: exact: max of Fractions
        exact = hi_of(A[1][1].upper()) + m*m/abs(n11u)
        check('schur', hi_of(schur) >= exact and lo_of(schur) <= exact, 'Schur bound does not contain the exact rational value')
        if not ok: continue
        nok += 1
        S = hi_of(schur)
        for _ in range(30):
            t = [Fr(rng.choice([-1, 1, rng.randint(-100, 100)]), 100 if rng.random() < 0.5 else 1) for _ in range(NV)]
            t = [min(max(x, Fr(-1)), Fr(1)) for x in t]   # RIGOR: exact: Fractions
            vals = []
            for (i, k) in ((0, 0), (0, 1), (1, 1)):
                T = A[i][k]
                cs = [F(c.mid()) if c.rad() == 0 else rng.choice([lo_of(c), hi_of(c)]) for c in T.c]   # RIGOR: exact: exact midpoint of a radius-0 ball, else an exact ball endpoint
                vals.append(tm_eval_exact(cs, t) + rng.choice([Fr(0), hi_of(T.R), -hi_of(T.R)]))
            a11, a12, a22 = vals
            check('schur', a11 <= n11u < 0, 'N11 sample above the bound')
            if a11 < 0: check('schur', a22 + a12*a12/abs(a11) <= S, f'Schur sample {float(a22 + a12*a12/abs(a11))} > bound {float(S)}')   # RIGOR: info: message
    check('schur', nok >= N//10, f'too few certified random Schur cases ({nok})')

def test_witness_tm():
    """fp.merged:321 (FARL FEM leaf): every Taylor-model range formed while recomputing the certificate in which some quadratic
    coefficient is a ball containing 0 (the case that produced NaN candidates before) must contain sampled exact values"""
    import fam5, gcert, tm as tmmod
    from gsigma import sigma_box
    L = rec('fp', 321); fam = fam5.FAMILIES['FP']()
    st = fam.fem_setup(np.array(L['p0']), np.array(L['r']), L['k'])   # RIGOR: exact: float64 arrays of JSON floats (test inputs)
    pb, pc, rc = st
    caught = []
    orig = tmmod.TM.range
    def spy(self):
        res = orig(self)
        if len(caught) < 300:
            for i in range(NV):
                a = self.c[IDX[(i, i)]]
                if a.lower() <= 0 <= a.upper() and not (a == 0): caught.append((self.copy(), res)); break
        return res
    tmmod.TM.range = spy
    try:
        gc = gcert.GCert(pb); sol = gc.solve(np.array(pc), np.array(rc)); sc, _ = sigma_box(gc, sol, np.array(pc), np.array(rc))   # RIGOR: exact: float64 arrays of JSON floats (test inputs)
        res = gc.certify(np.array(pc), np.array(rc), sc, sol=sol)   # RIGOR: exact: float64 arrays of JSON floats (test inputs)
    finally: tmmod.TM.range = orig
    check('witness-tm', bool(res.get('ok')), f'fp:321 no longer certifies: {res.get("why")}')
    for T, (lo, hi) in caught:
        L_, H_ = lo_of(lo), hi_of(hi)
        for v in tm_samples(T, 30): check('witness-tm', L_ <= v <= H_, 'range misses a sample')
    check('witness-tm', len(caught) > 0, 'fp:321: no Taylor-model range with a zero-containing quadratic coefficient was exercised (vacuous test)')
    print(f'   witness fp:321: {len(caught)} ranges with a quadratic coefficient ball containing 0 checked; certificate ok={res.get("ok")}')

# ------------------------------------------------------------------------------------------------ witnesses: load records
def rec(fam, line):
    with open(os.path.join(HERE, 'certs3', f'{fam}.merged.jsonl')) as f:
        for i, l in enumerate(f, 1):
            if i == line: return json.loads(l)
def test_witness_schur():
    """pp.merged:66327 (FEM leaf): recompute the certificate with the current code; check the Schur step of every N against
    the exact rational recomputation from the exact range endpoints"""
    import fam5, gcert
    from gsigma import sigma_box
    L = rec('pp', 66327); fam = fam5.FAMILIES['PP']()
    st = fam.fem_setup(np.array(L['p0']), np.array(L['r']), L['k'])   # RIGOR: exact: float64 arrays of JSON floats (test inputs)
    check('witness-schur', st is not None, 'FEM setup failed'); pb, pc, rc = st
    seen = []
    orig = gcert.schur_upper
    def spy(N):
        res = orig(N); seen.append((N, res)); return res
    gcert.schur_upper = spy
    try:
        gc = gcert.GCert(pb); sol = gc.solve(np.array(pc), np.array(rc)); sc, _ = sigma_box(gc, sol, np.array(pc), np.array(rc))   # RIGOR: exact: float64 arrays of JSON floats (test inputs)
        res = gc.certify(np.array(pc), np.array(rc), sc, sol=sol)   # RIGOR: exact: float64 arrays of JSON floats (test inputs)
    finally: gcert.schur_upper = orig
    check('witness-schur', bool(res.get('ok')), f'witness leaf no longer certifies: {res.get("why")}')
    for N, (ok, schur, why) in seen:
        if schur is None: continue
        n11 = hi_of(N[0][0].upper()); lo, hi = N[0][1].range(); m = max(abs(lo_of(lo)), abs(hi_of(hi)))   # RIGOR: exact: max of Fractions
        exact = hi_of(N[1][1].upper()) + m*m/abs(n11)
        check('witness-schur', hi_of(schur) >= exact and lo_of(schur) <= exact, 'Schur bound below the exact rational value')
    print(f'   witness pp:66327: {len(seen)} Schur steps checked, certificate ok={res.get("ok")}')

# ------------------------------------------------------------------------------------------------ strain factors
def test_strain(witness=True):
    """piece_factors (mu_q, j_q): for sample points p of the box, v^T G_q(p) v >= mu_q v^T G_q(pc) v and J_q(p) <= j_q J_q(pc),
    i.e. mu_q <= (J_q(p)/J_q(pc)) lambda_min(T T^T), T = E(pc) E(p)^{-1}, and j_q >= J_q(p)/J_q(pc) (300-bit references)"""
    import fam5, gcert
    from gsigma import piece_factors
    from poly import Poly
    cases = []
    if witness:
        L = rec('pp', 189234); cases.append(('pp:189234', L))
    with open(os.path.join(HERE, 'certs3', 'pp.merged.jsonl')) as f:
        fem = [json.loads(l) for i, l in zip(range(400000), f) if '"FEM"' in l]
    rng.shuffle(fem); cases += [('pp:FEM', x) for x in fem[:(6 if QUICK else 20)]]
    fam = fam5.FAMILIES['PP']()
    for name, L in cases:
        st = fam.fem_setup(np.array(L['p0']), np.array(L['r']), L['k'])   # RIGOR: exact: float64 arrays of JSON floats (test inputs)
        if st is None: check('strain', False, f'{name}: setup failed'); continue
        pb, pc, rc = st
        gc = gcert.GCert(pb)
        fac = piece_factors(gc, np.array(pc), np.array(pc), np.array(rc))   # RIGOR: exact: float64 arrays of JSON floats (test inputs)
        for _ in range(10 if QUICK else 25):
            t = [rng.choice([-1.0, 1.0, rng.uniform(-1, 1)]) for _ in range(4)]
            p = [pc[i] + rc[i]*t[i] for i in range(4)]
            def geo_at(pt):
                Z = [Poly.const(arb(x)) for x in pt]
                return pb.geometry(pb.crc, Z)
            with prec(300):
                g1 = geo_at(p); g0 = geo_at(pc)
                for q, ((mu, jq), a, b) in enumerate(zip(fac, g1, g0)):
                    E1 = [list(e.d.values())[0] if e.d else arb(0) for e in a[5]]; E0 = [list(e.d.values())[0] if e.d else arb(0) for e in b[5]]
                    d1 = E1[0]*E1[3] - E1[1]*E1[2]; d0 = E0[0]*E0[3] - E0[1]*E0[2]
                    # T = E0 E1^{-1}
                    i11, i12, i21, i22 = E1[3]/d1, -E1[1]/d1, -E1[2]/d1, E1[0]/d1
                    T11 = E0[0]*i11 + E0[1]*i21; T12 = E0[0]*i12 + E0[1]*i22
                    T21 = E0[2]*i11 + E0[3]*i21; T22 = E0[2]*i12 + E0[3]*i22
                    S11 = T11*T11 + T12*T12; S12 = T11*T21 + T12*T22; S22 = T21*T21 + T22*T22
                    tr = S11 + S22; det = S11*S22 - S12*S12
                    lmin = (tr - (tr*tr - 4*det).sqrt())/2
                    jr = abs(d1)/abs(d0)
                    check('strain', F(mu) <= lo_of(jr*lmin), f'{name} piece {q}: mu {mu} > reference {jr*lmin}')
                    check('strain', F(jq) >= hi_of(jr), f'{name} piece {q}: j {jq} < reference {jr}')
    print(f'   strain factors: {len(cases)} FEM setups (incl. witness pp:189234)')

# ------------------------------------------------------------------------------------------------ atan / cot endpoints
def test_angles():
    import fam5
    ts = [0.281951904296875, 0.02, 0.3, 0.30000000000000004, 1e-300, 5e-324] + [rng.uniform(0, 0.4) for _ in range(300 if not QUICK else 50)]
    for t in ts:
        r = ref(lambda: arb(t).atan())
        check('angles', F(fdown(arb(t).atan())) <= lo_of(r) and F(fup(arb(t).atan())) >= hi_of(r), f'atan({t})')
    # witness sp_spk:7: the BO2 arm of the leaf must contain atan of the whole t-range
    L = rec('sp_spk', 7); fam = fam5.FAMILIES['SP_spk']()
    A_, C_ = fam.arms(np.array(L['p0']), np.array(L['r']))   # RIGOR: exact: float64 arrays of JSON floats (test inputs)
    t_lo = F(L['p0'][0]) - F(L['r'][0]); t_hi = F(L['p0'][0]) + F(L['r'][0])
    with prec(300):
        a_lo = (arb(t_lo.numerator)/t_lo.denominator).atan(); a_hi = (arb(t_hi.numerator)/t_hi.denominator).atan()
    check('angles', F(A_.b1) <= lo_of(a_lo) and F(A_.b2) >= hi_of(a_hi), f'witness sp_spk:7 arm [{A_.b1}, {A_.b2}] misses atan range')
    print(f'   witness sp_spk:7: arm angles [{A_.b1!r}, {A_.b2!r}] contain [atan {float(t_lo)}, atan {float(t_hi)}]')   # RIGOR: info: message
    # SP predicates use the same directed conversion (tb in SP.preds): check via the mirror predicate's angle ball
    for bf in [rng.uniform(0.05, 1.5707) for _ in range(300 if not QUICK else 50)] + [1.5707963267948966, 0.1]:
        u_lo = fdown(1/arb(bf).tan()); r = ref(lambda: 1/arb(bf).tan())
        check('angles', F(u_lo) <= lo_of(r) and F(fup(1/arb(bf).tan())) >= hi_of(r), f'cot({bf})')

# ------------------------------------------------------------------------------------------------ FARL tail and cut
def exact_cut_max(lo, hi, k0):
    """max over u in [lo, hi] (Fractions) of g(u^2) = (1 - k0 u^2)^2 + k0^2 u^2, exactly (convex quadratic in s = u^2)"""
    def g(s): return (1 - k0*s)**2 + k0*k0*s   # RIGOR: exact: Fraction ** int (the exact reference)
    s_hi = max(lo*lo, hi*hi); s_lo = Fr(0) if lo <= 0 <= hi else min(lo*lo, hi*hi)   # RIGOR: exact: Fractions
    return max(g(s_lo), g(s_hi))   # RIGOR: exact: Fractions
def test_farl():
    import fam5, charts
    from fractions import Fraction
    # random closed-form checks of charts.cut_tail
    for _ in range(400 if not QUICK else 80):
        c = rng.uniform(-3, 6); h = rng.choice([0.0, 10**rng.uniform(-12, 0)]); k0 = Fraction(rng.randint(1, 4096), 4096)
        tail = charts.cut_tail(c, h, k0)
        m = exact_cut_max(F(c) - F(h), F(c) + F(h), k0)
        r = ref(lambda: arb.pi()*arb.pi()/(arb(m.numerator)/m.denominator))
        check('farl', F(tail) <= lo_of(r), f'cut_tail({c}, {h}, {k0}) = {tail} > reference')
        # brute-force samples: pi^2/g(u^2) >= tail at sample points
        for _ in range(5):
            u = F(c) + F(h)*Fr(rng.randint(-1000, 1000), 1000)
            gu = (1 - k0*u*u)**2 + k0*k0*u*u   # RIGOR: exact: Fraction ** int (the exact reference)
            rr = ref(lambda: arb.pi()*arb.pi()/(arb(gu.numerator)/gu.denominator))
            check('farl', F(tail) <= lo_of(rr), 'cut_tail above a sample value')
    # witnesses pp:80526 (PP, FARL) and fp:325 (FP, FARL): fresh setups; the tail is <= the exact value; the cut u-range contains
    # cot of the whole angle range; validity re-checked at sample u in exact rational arithmetic
    for fam_name, line in (('pp', 80526), ('fp', 325), ('fp', 321)):
        L = rec(fam_name, line); fam = fam5.FAMILIES[fam_name.upper()]()
        p0 = np.array(L['p0']); r = np.array(L['r'])   # RIGOR: exact: float64 arrays of JSON floats (test inputs)
        st = fam.fem_setup(p0, r, L['k'])
        if st is None: check('farl', False, f'{fam_name}:{line}: no FARL setup'); continue
        pb, pc, rc = st
        tail = pb.tails_fn(pc, rc)[0]
        k0 = None
        # recover k0 from the problem name (FARL_k.._RA.._L.._k0..) and recompute the exact maximum
        import re
        k0 = pb.farl[2] if hasattr(pb, 'farl') else None                 # the exact rational k0 of the FARL problem
        if k0 is not None:
            mx = exact_cut_max(F(pc[0]) - F(rc[0]), F(pc[0]) + F(rc[0]), k0)
            rr = ref(lambda: arb.pi()*arb.pi()/(arb(mx.numerator)/mx.denominator))
            check('farl', F(tail) <= lo_of(rr), f'{fam_name}:{line}: tail {tail} above the exact value')
            print(f'   witness {fam_name}:{line}: tail {tail!r} <= exact {rr.lower()} (old logged {L.get("tails")})')
        # u-range: the chart box contains [cot b2, cot b1] of the far arm's angle range
        if fam_name == 'pp': ib = 0 if fam._long(p0, r) == 'A' else 2
        else: ib = 0 if fam.Afar else 2
        b1 = F(p0[ib]) - F(r[ib]); b2 = F(p0[ib]) + F(r[ib])
        with prec(300):
            c1 = 1/(arb(b1.numerator)/b1.denominator).tan(); c2 = 1/(arb(b2.numerator)/b2.denominator).tan()
        check('farl', F(pc[0]) - F(rc[0]) <= lo_of(c2) and F(pc[0]) + F(rc[0]) >= hi_of(c1), f'{fam_name}:{line}: u-range misses cot range')

# ------------------------------------------------------------------------------------------------ FARL validity H1-H3 (Lemma 9.6)
def farl_conditions_exact(RA, L, k0, u, cpts):
    """Lemma 9.6 at one exact u: (H1) R > L > k0 u; with d = Q2 - Q1 = (1 - k0 u^2, -k0 u) and S(P) = det(d, P - Q1):
    (H2) S(B), S(D) and S(corner) for every corner of the C-box are nonzero with the sign of S(B); (H3) S(R (u, 1)) has the
    opposite sign.  Also det(d, (u, 1)) = 1 (the identity behind S(B) = -L, S(D) = -(L - k0 u), S(s (u, 1)) = s - L).
    Exact rational arithmetic."""
    q1 = (L*u, L); d = (1 - k0*u*u, -k0*u)
    S = lambda P: d[0]*(P[1] - q1[1]) - d[1]*(P[0] - q1[0])
    h1 = RA > L > k0*u
    sB = S((Fr(0), Fr(0)))
    h2 = sB != 0 and all(S(P)*sB > 0 for P in [(Fr(1), Fr(0))] + cpts)
    h3 = S((RA*u, RA))*sB < 0
    ident = (d[0]*1 - d[1]*u == 1) and sB == -L and S((Fr(1), Fr(0))) == -(L - k0*u)
    return h1, h2, h3, ident
def farl_conditions_interval(RA, L, k0, u_lo, u_hi, cpts):
    """Lemma 9.6 over the WHOLE interval u in [u_lo, u_hi] and the whole C-box (corners cpts), exactly: with
    S_P(u) = det(Q2 - Q1, P - Q1) = -k0 P_y u^2 + k0 P_x u + (P_y - L) (affine in P, so its maximum over the box is at a corner;
    a quadratic in u, so its maximum over the interval is at an end or at the vertex u* = P_x / (2 P_y) when that lies inside):
    (H1) R > L > k0 u_hi;  (H2) S_B = -L < 0, max_u S_D(u) = k0 u_hi - L < 0 and max over corners and u of S_P(u) < 0;
    (H3) S_{R (u,1)}(u) = R - L > 0 for every u.  Returns (H1, H2, H3, worst S over the C-box)."""
    def S(P, u): return -k0*P[1]*u*u + k0*P[0]*u + (P[1] - L)
    def smax(P):
        cand = [u_lo, u_hi]
        if P[1] != 0:
            us = P[0]/(2*P[1])
            if u_lo <= us <= u_hi: cand.append(us)
        return max(S(P, u) for u in cand)   # RIGOR: exact: max of Fractions
    h1 = RA > L > k0*u_hi and k0 > 0
    worst = max(smax(P) for P in cpts)   # RIGOR: exact: max of Fractions
    h2 = -L < 0 and k0*u_hi - L < 0 and worst < 0
    h3 = RA - L > 0
    return h1, h2, h3, worst

def test_farl_validity():
    import fam5
    cases = [('pp', 80526), ('fp', 325), ('fp', 321)]
    for fam_name in ('pp', 'fp'):
        with open(os.path.join(HERE, 'certs3', f'{fam_name}.merged.jsonl')) as f:
            farl = [json.loads(l) for i, l in zip(range(700000), f) if '"FEM"' in l and '"pc": [' in l]
        farl = [d for d in farl if d['pc'][1] == 0.0]
        rng.shuffle(farl); cases += [(fam_name, d) for d in farl[:(5 if QUICK else 25)]]
    n_setups = 0
    for fam_name, line in cases:
        L_ = rec(fam_name, line) if isinstance(line, int) else line
        fam = fam5.FAMILIES[fam_name.upper()]()
        st = fam.fem_setup(np.array(L_['p0']), np.array(L_['r']), L_['k'])   # RIGOR: exact: float64 arrays of JSON floats (test inputs)
        if st is None: check('farl-H', False, f'{fam_name}: no FARL setup'); continue
        pb, pc, rc = st
        if not hasattr(pb, 'farl'): check('farl-H', False, f'{fam_name}: not a FARL problem ({pb.name})'); continue
        RA, Lc, k0 = pb.farl; n_setups += 1
        u_lo, u_hi = F(pc[0]) - F(rc[0]), F(pc[0]) + F(rc[0])
        cx = (F(pc[2]) - F(rc[2]), F(pc[2]) + F(rc[2])); cy = (F(pc[3]) - F(rc[3]), F(pc[3]) + F(rc[3]))
        cpts = [(x, y) for x in cx for y in cy]
        i1, i2, i3, worst = farl_conditions_interval(RA, Lc, k0, u_lo, u_hi, cpts)
        check('farl-H', i1 and i2 and i3, f'{fam_name}: Lemma 9.6 fails on the interval: H1 {i1} H2 {i2} H3 {i3} (worst S {float(worst)})')   # RIGOR: info: message
        # cross-check of the closed form: the direct determinant at the ends and at sampled interior points
        us = [u_lo, u_hi] + [u_lo + (u_hi - u_lo)*Fr(rng.randint(0, 1000), 1000) for _ in range(8)]
        for u in us:
            h1, h2, h3, ident = farl_conditions_exact(RA, Lc, k0, u, cpts)
            check('farl-H', h1 and h2 and h3 and ident, f'{fam_name}: Lemma 9.6 condition fails at u = {float(u)}: H1 {h1} H2 {h2} H3 {h3} identities {ident}')   # RIGOR: info: message
            q1 = (Lc*u, Lc); d = (1 - k0*u*u, -k0*u)
            for P in cpts: check('farl-H', d[0]*(P[1] - q1[1]) - d[1]*(P[0] - q1[0]) == -k0*P[1]*u*u + k0*P[0]*u + (P[1] - Lc), 'closed form of S_P')
    check('farl-H', n_setups >= 3, f'only {n_setups} FARL setups were tested (vacuous)')
    print(f'   FARL validity (H1-H3 exact over the whole u-interval and C-box via the analytic extremes; closed form cross-checked at 10 '
          f'values of u): {n_setups} setups incl. pp:80526, fp:325, fp:321')

# ------------------------------------------------------------------------------------------------ BO2 projection extrema
def test_bo2_proj():
    import fam5, bo2
    from bo2 import proj_range
    with open(os.path.join(HERE, 'certs3', 'pp.merged.jsonl')) as f:
        recs = [json.loads(l) for i, l in zip(range(300000), f) if '"BO2"' in l]
    rng.shuffle(recs)
    fam = fam5.FAMILIES['PP']()
    for L in recs[:(40 if QUICK else 200)]:
        armA, armC = bo2.arm_variant(*fam.arms(np.array(L['p0']), np.array(L['r'])), L.get('arms', 'polar'))   # RIGOR: exact: float64 arrays of JSON floats (test inputs)
        th = arb(L['theta']); u = (th.cos(), th.sin())
        pts = [(arb(0), arb(0)), (arb(1), arb(0))]
        for a in (armA, armC):
            if a is not None and a.kind == 'polar': pts += a.outer_points()
        lo, hi, prc = proj_range(pts, u)
        for x in prc: check('bo2', F(lo) <= lo_of(x) and F(hi) >= hi_of(x), 'projection outside the directed range')

# ------------------------------------------------------------------------------------------------ spike tail
def test_spike_tail():
    import charts
    for _ in range(300 if not QUICK else 60):
        t = rng.uniform(1e-4, 0.4); h = rng.choice([0.0, rng.uniform(0, t/2)])
        tail = charts.spike_tail([t], [h])
        th = F(t) + F(h)
        r = ref(lambda: arb.pi()*arb.pi()/((arb(th.numerator)/th.denominator)**2))   # RIGOR: info: reference only (positive ball)
        check('spike', F(tail) <= lo_of(r), f'spike_tail({t}, {h})')

def main():
    import time
    t0 = time.time()
    for name, fn in (('helpers', test_helpers), ('tm', test_tm), ('schur', test_schur), ('angles', test_angles),
                     ('farl', test_farl), ('bo2', test_bo2_proj), ('spike', test_spike_tail), ('witness-schur', test_witness_schur),
                     ('witness-tm', test_witness_tm), ('farl-H', test_farl_validity),
                     ('strain', test_strain)):
        t1 = time.time(); fn()
        p, f = STATS.get(name, [0, 0])
        print(f'{name}: {p} checks passed, {f} failed ({time.time() - t1:.0f} s)', flush=True)
    total_fail = sum(f for p, f in STATS.values())
    print(f'test_containment: {"PASS" if total_fail == 0 else "FAIL"} ({sum(p for p, f in STATS.values())} checks, {total_fail} failures, '
          f'{time.time() - t0:.0f} s)')
    sys.exit(0 if total_fail == 0 else 1)

if __name__ == '__main__':
    main()
