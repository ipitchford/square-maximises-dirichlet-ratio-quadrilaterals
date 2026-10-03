"""Fail-closed independent replay of a cover log (cover2 legacy logs, cover3/cover4 path logs).

Usage: replay2.py LOG FAMILY MANIFEST [--shard I/N] [--out RESULTS.jsonl] [--no-recert]
  FAMILY   one of FAMILIES below (fixes the problem charts, the certificate methods allowed and the discard predicates)
  MANIFEST JSON file with key "roots": [[root_id, p0, r], ...] -- the roots the run claims (from the run script's manifest mode)

Checks (any violation -> exit status 1, nothing is "skipped"):
 1. every line parses as a JSON object; record types are whitelisted (leaf, split, esc, root_done); per-type/per-status schemas
    (required keys with types; no unknown keys); leaf statuses whitelisted (cert, discard) -- defer, FAIL or anything else fails;
 2. no duplicate leaf/split/esc paths (legacy logs: no duplicate (root, p0, r)); no duplicate root_done; every manifest root has
    exactly one root_done with ok = true and the recorded box equal to the dyadic snapping of the manifest box; no record refers to
    a root outside the manifest;
 3. the leaves of every root tile its box exactly (dyadic sub-boxes, pairwise disjoint, exact rational volumes sum to the root);
 4. every discard is re-decided from scratch by regions.check_reason (rigorous predicate named in the record);
 5. every certificate is recomputed from scratch with the method named in the record (FEM: fresh float solve, verified sigma2,
    rigorous certificate; BO: fresh rigorous 1-D bound with the recorded slicing axis and cell count, fresh rigorous P2 test).
Every leaf box must have exact float ends (nanguard.exact_box); the manifest has a strict schema (its roots are checked against
their generators, and their snapped boxes against the nominal cells, in accounting.py).
Output (--out): one JSON line per checked leaf: root, path, status, ok, record sha256 (canonical JSON of the log record),
binding sha256 (binding.py: transitive sources, chart cache, manifest, numerics environment); then a summary line (with the
sha256 of the sorted multiset of checked record hashes).  Exit 0 iff all ok."""
import sys, os, json, time, hashlib, argparse, math, numpy as np, warnings
from fractions import Fraction as Fr
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)

MANIFEST_KEYS = {'COLLAR': {'family', 'roots'}, 'SPIKE_X1': {'family', 'roots'}, 'MAIN_CBOX': {'family', 'roots'}, 'BOX': {'family', 'roots'},
                 'KF_COLLAR': {'family', 'mode', 'region', 't_slabs', 'cells', 'outside', 'roots'},
                 'KF_SPIKE': {'family', 'mode', 'region', 't_slabs', 'cells', 'outside', 'roots'}}

def rec_hash(d): return hashlib.sha256(json.dumps(d, sort_keys=True, separators=(',', ':')).encode()).hexdigest()

# ------------------------------------------------------------------------------------------------ families
def _fam():
    from charts import bd_problem, spike_problem, spike2_problem
    return {
        'COLLAR':     dict(pf=spike_problem, bo=None, legacy=True),
        'SPIKE_X1':   dict(pf=lambda k: spike2_problem(k, 1.0), bo=None, legacy=True),
        'MAIN_CBOX':  dict(pf=bd_problem, bo=None, legacy=True),
        'BOX':        dict(pf=bd_problem, bo=None, legacy=True),
        'KF_COLLAR':  dict(pf=spike_problem, bo=None, legacy=False),
        'KF_SPIKE':   dict(pf=lambda k: spike2_problem(k, 1.0), bo=None, legacy=False),
        'K1X_COLLAR': dict(pf=spike_problem, bo='spike', legacy=False),
        'K1X_SPIKE':  dict(pf=lambda k: spike2_problem(k, 1.0), bo='spike', legacy=False),
        'BD_HYB':     dict(pf=bd_problem, bo='bd', legacy=False),
    }

NUM = (int, float)
def is_vec4(v): return isinstance(v, list) and len(v) == 4 and all(isinstance(x, NUM) and np.isfinite(x) for x in v)
SCHEMA = {
    ('leaf', 'discard'): dict(req={'type': str, 'root': str, 'status': str, 'reason': str, 'p0': is_vec4, 'r': is_vec4},
                              opt={'path': str}),
    ('leaf', 'cert', 'FEM'): dict(req={'type': str, 'root': str, 'status': str, 'p0': is_vec4, 'r': is_vec4, 'k': int,
                                       'L1': list, 'schur': NUM, 'h': NUM, 'temple': NUM},
                                  opt={'path': str, 'method': str, 'sigma': (list, type(None)), 'tails': (list, type(None))}),
    ('leaf', 'cert', 'BO'): dict(req={'type': str, 'root': str, 'status': str, 'method': str, 'p0': is_vec4, 'r': is_vec4,
                                      'k': int, 'E': NUM, 'K': int, 'schur': NUM},
                                 opt={'path': str, 'theta': (NUM, type(None)), 'lam2p2': (NUM, type(None))}),
    ('split',): dict(req={'type': str, 'path': str, 'root': str, 'dirn': int, 'knew': int}, opt={}),
    ('esc',): dict(req={'type': str, 'path': str, 'root': str, 'knew': int}, opt={}),
    ('root_done',): dict(req={'type': str, 'root': str, 'ok': bool, 'box': dict}, opt={'time': NUM, 'counts': dict}),
}
def check_schema(d):
    t = d.get('type')
    if t == 'leaf':
        st = d.get('status')
        if st == 'discard': key = ('leaf', 'discard')
        elif st == 'cert': key = ('leaf', 'cert', d.get('method', 'FEM'))
        else: return f'leaf status {st!r} not whitelisted'
    elif t in ('split', 'esc', 'root_done'): key = (t,)
    else: return f'record type {t!r} not whitelisted'
    if key not in SCHEMA: return f'no schema for {key}'
    sc = SCHEMA[key]
    for k, ty in sc['req'].items():
        if k not in d: return f'missing key {k}'
        v = d[k]
        ok = ty(v) if callable(ty) and not isinstance(ty, type) and not isinstance(ty, tuple) else isinstance(v, ty)
        if isinstance(v, bool) and ty in (int, NUM) : ok = False
        if not ok: return f'bad type for {k}: {v!r}'
    for k in d:
        if k not in sc['req'] and k not in sc['opt']: return f'unknown key {k}'
        if k in sc['opt'] and not isinstance(d[k], sc['opt'][k]): return f'bad type for optional {k}'
    if t == 'leaf' and d['status'] == 'cert' and d.get('method', 'FEM') not in ('FEM', 'BO'): return 'unknown method'
    return None

def dyadic_box(p0, r):
    from grid import dyadic_root
    return dyadic_root(p0, r)

def tiles(root, leaves):
    p0 = [Fr(x) for x in root['p0']]; r = [Fr(x) for x in root['r']]
    act = [i for i in range(4) if r[i] != 0]
    vol_root = math.prod([2*r[i] for i in act]) if act else Fr(1)              # Fractions: exact
    boxes = []
    for L in leaves:
        q = [Fr(x) for x in L['p0']]; s = [Fr(x) for x in L['r']]
        for i in range(4):
            if r[i] == 0:
                if s[i] != 0 or q[i] != p0[i]: return False, 'inactive coordinate moved'
                continue
            if s[i] <= 0: return False, 'nonpositive half-width'
            ratio = r[i]/s[i]
            if ratio.denominator != 1 or (ratio.numerator & (ratio.numerator - 1)) != 0: return False, 'non-dyadic width'
            off = (q[i] - s[i] - (p0[i] - r[i]))/(2*s[i])
            if off.denominator != 1 or off < 0 or off >= ratio: return False, 'misaligned or outside root'
        boxes.append((q, s))
    # disjointness: sort by first active lower corner, sweep (dyadic boxes overlap iff they overlap in every active coordinate)
    if act:
        order = sorted(range(len(boxes)), key=lambda j: boxes[j][0][act[0]] - boxes[j][1][act[0]])
        for ii in range(len(order)):
            qa, sa = boxes[order[ii]]
            for jj in range(ii + 1, len(order)):
                qb, sb = boxes[order[jj]]
                if qb[act[0]] - sb[act[0]] >= qa[act[0]] + sa[act[0]]: break
                if all(abs(qa[i] - qb[i]) < sa[i] + sb[i] for i in act): return False, 'overlapping leaves'
    vol = sum((math.prod([2*s[i] for i in act]) for (q, s) in boxes), Fr(0)) if act else Fr(len(boxes))
    if vol != vol_root: return False, f'volume {float(vol)} != {float(vol_root)}'   # RIGOR: info: message
    return True, 'ok'

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('log'); ap.add_argument('family'); ap.add_argument('manifest')
    ap.add_argument('--shard', default='0/1'); ap.add_argument('--out', default=None); ap.add_argument('--no-recert', action='store_true')
    a = ap.parse_args()
    FAM = _fam()
    if a.family not in FAM: print('unknown family', a.family); sys.exit(1)
    fam = FAM[a.family]
    si, sn = [int(x) for x in a.shard.split('/')]
    from binding import binding
    from nanguard import exact_box, InexactBox
    import fpenv, replay_common as rc
    sc = rc.start_checks()                                   # flags, FP environment, pinned chart caches (before any load)
    if sc: print('START CHECKS FAILED:', sc); sys.exit(1)
    if not (0 <= si < sn): print('bad shard', a.shard); sys.exit(1)
    msnap = rc.snapshot(a.manifest); lsnap = rc.snapshot(a.log)      # the replay parses exactly these bytes
    vh, vdet = binding('replay2.py', msnap['sha256'], a.family)
    bad = []
    def fail(msg):
        bad.append(msg)
        if len(bad) <= 50: print('FAIL:', msg, flush=True)
    try: man = json.loads(msnap['bytes'])
    except Exception: man = None
    if not isinstance(man, dict) or set(man) != MANIFEST_KEYS.get(a.family, {'family', 'roots'}): fail('manifest: unexpected keys'); man = dict(family=None, roots=[])
    if man.get('family') != a.family: fail(f'manifest family {man.get("family")!r} differs from the replayed family {a.family!r}')
    mroots = {}
    for c in man['roots']:
        if not (isinstance(c, list) and len(c) == 3 and isinstance(c[0], str) and is_vec4(c[1]) and is_vec4(c[2])):
            fail('manifest: malformed root entry'); continue
        rid, p0, r = c
        if rid in mroots: fail(f'duplicate manifest root {rid}')
        if not all(x >= 0 for x in r): fail(f'manifest root {rid}: negative radius')
        # a processed root is represented by its exact outward dyadic snapping (grid.dyadic_root), which its leaves tile; its own
        # float encoding need not have exact ends (accounting.py checks that the snapped boxes cover the nominal cells)
        mroots[rid] = (p0, r)
    leaves = {}; seen_paths = set(); seen_keys = set(); done = {}; nlines = 0
    if True:
        for ln, line in enumerate(lsnap['bytes'].decode().split('\n')[:-1] if lsnap['bytes'].endswith(b'\n') else lsnap['bytes'].decode().split('\n'), 1):
            nlines += 1
            s = line.strip()
            if not s: fail(f'line {ln}: empty line'); continue
            try: d = json.loads(s)
            except Exception: fail(f'line {ln}: malformed JSON'); continue
            if not isinstance(d, dict): fail(f'line {ln}: not an object'); continue
            why = check_schema(d)
            if why: fail(f'line {ln}: {why}'); continue
            root = d['root']
            if root not in mroots: fail(f'line {ln}: root {root} not in manifest'); continue
            t = d['type']
            if t == 'root_done':
                if root in done: fail(f'line {ln}: duplicate root_done {root}')
                done[root] = d; continue
            if 'path' in d:
                pk = (root, d['path'])
                if pk in seen_paths: fail(f'line {ln}: duplicate path {d["path"]} in root {root}')
                seen_paths.add(pk)
                if not d['path'].startswith(root) or any(c not in '01k' for c in d['path'][len(root):]):
                    fail(f'line {ln}: path {d["path"]} is not root + bisection tags')
            elif not fam['legacy']:
                fail(f'line {ln}: record without path in a path log')
            if t == 'leaf':
                try: exact_box(d['p0'], d['r'])
                except InexactBox as e: fail(f'line {ln}: leaf box {e}'); continue
                key = (root, tuple(d['p0']), tuple(d['r']))
                if key in seen_keys: fail(f'line {ln}: duplicate leaf box')
                seen_keys.add(key)
                d['_line'] = ln; d['_sha'] = rec_hash({k: v for k, v in d.items() if not k.startswith('_')})
                leaves.setdefault(root, []).append(d)
    # roots: completion, recorded boxes, tiling
    for rid, (p0, r) in mroots.items():
        if rid not in done: fail(f'root {rid}: no root_done (pending)'); continue
        rd = done[rid]
        if rd['ok'] is not True: fail(f'root {rid}: root_done ok = {rd["ok"]}')
        sp0, sr = dyadic_box(p0, r)
        box = rd['box']
        if [float(x) for x in box.get('p0', [])] != [float(x) for x in sp0] or [float(x) for x in box.get('r', [])] != [float(x) for x in sr]:   # RIGOR: exact: float() of JSON floats; exact equality
            fail(f'root {rid}: recorded box differs from the snapped manifest box')
        ok, why = tiles(dict(p0=sp0, r=sr), leaves.get(rid, []))
        if not ok: fail(f'root {rid}: tiling: {why}')
    # re-verification of leaves (sharded by root order)
    from regions import check_reason
    out = open(a.out, 'x') if a.out else None                # exclusive creation: never overwrites an output
    gcs = {}
    def gc_for(k):
        from gcert import GCert
        if k not in gcs: gcs[k] = GCert(fam['pf'](k))
        return gcs[k]
    nchk = 0; nok = 0; t0 = time.time(); shas = []
    for j, rid in enumerate(sorted(mroots)):
        if j % sn != si: continue
        for L in leaves.get(rid, []):
            p0 = np.array(L['p0'], float); r = np.array(L['r'], float)   # RIGOR: exact: float64 arrays of Python floats (no rounding)
            res_ok = False; info = ''
            try:
                if L['status'] == 'discard':
                    res_ok = bool(check_reason(L['reason'], p0, r)); info = L['reason']
                elif a.no_recert:
                    res_ok = True; info = 'not recertified (--no-recert)'
                elif L.get('method', 'FEM') == 'FEM':
                    from gsigma import sigma_box
                    gc = gc_for(L['k']); sol = gc.solve(p0, r)
                    sc, sinfo = sigma_box(gc, sol, p0, r)
                    res = gc.certify(p0, r, sc, sol=sol) if sc is not None else dict(ok=False, why='sigma')
                    res_ok = bool(res.get('ok')); info = res.get('why', '') if not res_ok else f"schur {res['schur']:.6g}"
                else:
                    if fam['bo'] is None: raise ValueError('BO certificate in a family without BO')
                    from bocert import bo_cert_bd, bo_cert_spike
                    gc = gc_for(L['k'])
                    if fam['bo'] == 'bd':
                        if L.get('theta') is None: raise ValueError('BO bd leaf without theta')
                        from bocert import bo_lower_bd
                        from p2up import solve_p2, p2_test
                        E, th = bo_lower_bd(p0, r, L['K'], dirs=[L['theta']])
                        sol = solve_p2(gc, p0, r)
                        res = p2_test(gc, p0, r, sol, E) if E > 0 else dict(ok=False, why='E <= 0')
                    else:
                        res = bo_cert_spike(gc, p0, r, L['K'])
                    res_ok = bool(res.get('ok')); info = res.get('why', '') if not res_ok else f"E {res.get('E', 0):.6g}"
            except Exception as e:
                res_ok = False; info = 'exception: ' + repr(e)[:200]
            try: fpenv.check()
            except fpenv.FPEnvError as e: res_ok = False; info = str(e)[:200]
            nchk += 1; nok += res_ok; shas.append(L['_sha'])
            if not res_ok: fail(f'root {rid} path {L.get("path")} line {L["_line"]}: {L["status"]} not re-verified: {info}')
            if out:
                out.write(json.dumps(dict(root=rid, path=L.get('path'), line=L['_line'], status=L['status'],
                                          method=L.get('method', 'FEM' if L['status'] == 'cert' else None), ok=res_ok, info=info,
                                          record_sha256=L['_sha'], binding_sha256=vh, family=a.family)) + '\n')
    summ = dict(log=a.log, family=a.family, shard=a.shard, partial=False, recert=not a.no_recert, lines=nlines, roots=len(mroots), leaves_checked=nchk, ok=nok,
                failures=len(bad), binding_sha256=vh, manifest_sha256=vdet.get('manifest_sha256'), environment=vdet['environment'],
                leaves_sha256=hashlib.sha256(''.join(sorted(shas)).encode()).hexdigest(),
                log_sha256=lsnap['sha256'],
                time_s=round(time.time() - t0, 1))   # RIGOR: info: timing
    for w in rc.end_checks('replay2.py', msnap['sha256'], a.family, vh, [msnap, lsnap], vdet): fail('end of run: ' + w)
    summ['failures'] = len(bad); summ['late_loaded'] = dict(rc.LATE); summ['niceness'] = os.nice(0)   # (recorded, not bound)
    summ['binding_sources'] = vdet['sources']; summ['binding_refcache'] = vdet['refcache']
    print(json.dumps(summ), flush=True)
    if out: out.write(json.dumps(dict(summary=summ)) + '\n'); out.close()
    sys.exit(0 if not bad else 1)

if __name__ == '__main__':
    main()
