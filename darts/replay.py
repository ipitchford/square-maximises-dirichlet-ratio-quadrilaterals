"""Independent replay of a cover log (JSONL written by cover2.Cover2).

For every root: (1) the leaves tile the root exactly: every leaf is a dyadic sub-box of the root (aligned, power-of-two
half-widths), no two leaves overlap, and the exact (Fraction) volumes add up to the root volume; (2) every DISCARDED leaf is
re-checked with the region predicate; (3) every CERTIFIED leaf is re-certified from scratch: fresh float solve at the box centre,
fresh verified sigma2 (gsigma.sigma_box: LDL^T of the piecewise Loewner minorant), fresh rigorous certificate (gcert.GCert.certify).
Usage: replay.py LOG REGION [--sample N]   REGION in {MAIN, SPIKE, COLLAR}  (selects the problem family and discard predicate)."""
import sys, json, argparse, numpy as np, warnings, time
from fractions import Fraction as Fr
warnings.filterwarnings('ignore')
sys.path.insert(0, '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts')
from gcert import GCert
from gsigma import sigma_box
from charts import bd_problem, spike2_problem, spike_problem
from regions import main_discard, X0, check_reason

def problem_for(region):
    if region == 'MAIN': return bd_problem, main_discard
    if region == 'SPIKE': return (lambda k: spike2_problem(k, X0)), None
    if region == 'SPIKE_X1': return (lambda k: spike2_problem(k, 1.0)), None
    if region == 'COLLAR': return spike_problem, None
    if region in ('MAIN_CBOX', 'BOX'): return bd_problem, None
    raise ValueError(region)

def tiles(root, leaves):
    p0 = [Fr(x) for x in root['p0']]; r = [Fr(x) for x in root['r']]
    act = [i for i in range(4) if r[i] != 0]
    vol_root = np.prod([2*r[i] for i in act]) if act else Fr(1)
    boxes = []
    for L in leaves:
        q = [Fr(x) for x in L['p0']]; s = [Fr(x) for x in L['r']]
        for i in range(4):
            if r[i] == 0:
                if s[i] != 0 or q[i] != p0[i]: return False, 'inactive coordinate moved'
                continue
            ratio = r[i]/s[i]
            if ratio.denominator != 1 or (ratio.numerator & (ratio.numerator - 1)) != 0: return False, 'non-dyadic width'
            off = (q[i] - s[i] - (p0[i] - r[i]))/(2*s[i])
            if off.denominator != 1 or off < 0 or off >= ratio: return False, 'misaligned or outside root'
        boxes.append((q, s))
    # pairwise disjoint (dyadic boxes: overlap iff the intervals overlap in every active coordinate)
    for a in range(len(boxes)):
        qa, sa = boxes[a]
        for b in range(a + 1, len(boxes)):
            qb, sb = boxes[b]
            if all(abs(qa[i] - qb[i]) < sa[i] + sb[i] for i in act): return False, f'overlap {a} {b}'
    vol = sum(np.prod([2*s[i] for i in act]) for (q, s) in boxes) if act else Fr(len(boxes))
    if vol != vol_root: return False, f'volume {float(vol)} != {float(vol_root)}'
    return True, 'ok'

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('log'); ap.add_argument('region'); ap.add_argument('--sample', type=int, default=0)
    a = ap.parse_args()
    pfn, disc = problem_for(a.region)
    roots = {}; leaves = {}
    for line in open(a.log):
        d = json.loads(line)
        if d.get('type') == 'leaf': leaves.setdefault(d['root'], []).append(d)
        elif d.get('type') == 'root_done': roots[d['root']] = d
    gcs = {}
    n_ok = 0; n_bad = 0; t0 = time.time()
    rng = np.random.default_rng(0)
    for nr, (rid, rd) in enumerate(roots.items()):
        if nr % 32 == 0: print(json.dumps(dict(progress=nr, of=len(roots), ok=n_ok, bad=n_bad, t=round(time.time() - t0, 1))), flush=True)
        Ls = leaves.get(rid, [])
        # root box = recorded in root_done (cover2 records it) or reconstructed as the bounding box of the leaves
        root = rd.get('box')
        if root is None:
            lo = np.min([np.array(L['p0']) - np.array(L['r']) for L in Ls], axis=0); hi = np.max([np.array(L['p0']) + np.array(L['r']) for L in Ls], axis=0)
            root = dict(p0=((lo + hi)/2).tolist(), r=((hi - lo)/2).tolist())
        ok, why = tiles(root, Ls)
        if not ok: print('TILING FAIL', rid, why); n_bad += 1
        if not rd.get('ok', False): print('ROOT NOT OK in log', rid); n_bad += 1
        todo = Ls if not a.sample else [Ls[i] for i in rng.choice(len(Ls), size=min(a.sample, len(Ls)), replace=False)]
        for L in todo:
            p0 = np.array(L['p0']); r = np.array(L['r'])
            if L['status'] == 'discard':
                if not check_reason(L['reason'], p0, r): print('DISCARD FAIL', rid, L['p0'], L['r'], L['reason']); n_bad += 1
                else: n_ok += 1
                continue
            if L['status'] != 'cert': print('NON-CERT LEAF', rid, L['status']); n_bad += 1; continue
            k = L['k']
            if k not in gcs: gcs[k] = GCert(pfn(k))
            gc = gcs[k]
            sol = gc.solve(p0, r)
            sc, info = sigma_box(gc, sol, p0, r)
            res = gc.certify(p0, r, sc, sol=sol) if sc is not None else dict(ok=False, why='sigma')
            if res.get('ok'): n_ok += 1
            else: print('REPLAY FAIL', rid, L['p0'], L['r'], res.get('why')); n_bad += 1
    print(json.dumps(dict(log=a.log, region=a.region, roots=len(roots), leaves_checked=n_ok + n_bad, ok=n_ok, bad=n_bad,
                          time_s=round(time.time() - t0, 1))))

if __name__ == '__main__':
    main()
