"""Fail-closed independent replay of a cover5 log (families of fam5.py).
Usage: replay3.py LOG FAMILY MANIFEST [--shard I/N] [--out RESULTS.jsonl] [--no-recert]

Checks (any violation -> exit status 1):
 0. strict manifest schema (keys family / grid / roots / skipped only; grid = run_fam.GRIDS[family]; roots [id, p0, r], skipped
    [id, p0, r, reason]); the roots and skipped cells are exactly the cells of the family grid: same ids, and p0 AND r equal
    (float equality) to the grid encoding run_fam.axis_cells, whose boxes have exact float ends (nanguard.exact_box); each
    skipped cell's reason is re-checked with the family predicate on that exact box;
 1. every log line is a JSON object of a whitelisted type (leaf, split, esc, root_done) with the per-type / per-status / per-method
    schema (required keys and types, no unknown keys); leaf statuses: cert (methods BO2, INC, FEM) or discard -- anything else fails;
 2. no duplicate (root, path); no duplicate root_done / handoff / item_done; every manifest root has exactly one root_done with
    ok = true (or, for roots split into work items by pool.py, an item_done of the root item, every handed-off item done, all
    item_done ok = true; completeness is then the exact tiling of point 3) and the
    recorded box equal to the dyadic snapping of the manifest box; no record outside the manifest roots;
 3. the leaves of every root tile its (snapped) box exactly (dyadic, disjoint, exact volume);
 4. discards re-decided with the family predicate named in the record; BO2 certificates recomputed (fresh rigorous 1-D bound with
    the recorded slicing axis and cell count, fresh inner dart, fresh rigorous P2 test at the recorded mesh level); FEM
    certificates recomputed (fresh Cartesian enclosure, which must equal the recorded one, fresh solve / sigma2 / certificate);
    INC certificates (family PP only): the recorded outer / inner apexes re-checked against fresh rigorous angle ranges of the box
    (inc.verify_inclusion), a fresh L1 lower bound at the outer dart (CR level k) that must be >= the recorded E, and a fresh P2
    test of the inner dart (P2 level kp) with the recorded E.
Every leaf box must have exact float ends (nanguard.exact_box), so that every predicate / certificate works with the true box.
Output: per checked record a JSON line with the record sha256 and the binding sha256 (binding.py: transitive sources, chart
cache, manifest, numerics environment); summary line (with the sha256 of the sorted multiset of checked record hashes);
exit 0 iff all ok."""
import sys, os, json, time, hashlib, argparse, numpy as np, warnings
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
from replay2 import tiles, rec_hash
from binding import binding
from nanguard import exact_box, InexactBox

NUM = (int, float)
def vec4(v): return isinstance(v, list) and len(v) == 4 and all(isinstance(x, NUM) and not isinstance(x, bool) and np.isfinite(x) for x in v)
def vec2(v): return isinstance(v, list) and len(v) == 2 and all(isinstance(x, NUM) and np.isfinite(x) for x in v)
def opt_vec2(v): return v is None or vec2(v)
SCHEMA = {
    ('leaf', 'discard'): ({'type': str, 'path': str, 'root': str, 'status': str, 'reason': str, 'p0': vec4, 'r': vec4}, {}),
    ('leaf', 'cert', 'BO2'): ({'type': str, 'path': str, 'root': str, 'status': str, 'method': str, 'p0': vec4, 'r': vec4,
                               'k': int, 'E': NUM, 'theta': NUM, 'K': int, 'inner': int, 'Ain': opt_vec2, 'Cin': vec2,
                               'schur': NUM, 'lam2p2': NUM}, {'arms': str}),
    ('leaf', 'cert', 'FEM'): ({'type': str, 'path': str, 'root': str, 'status': str, 'method': str, 'p0': vec4, 'r': vec4,
                               'k': int, 'pc': vec4, 'rc': vec4, 'L1': list, 'schur': NUM, 'h': NUM, 'temple': NUM},
                              {'sigma': (list, type(None)), 'tails': (list, type(None))}),
    ('leaf', 'cert', 'INC'): ({'type': str, 'path': str, 'root': str, 'status': str, 'method': str, 'p0': vec4, 'r': vec4,
                               'k': int, 'kp': int, 'Aout': vec2, 'Cout': vec2, 'Ain': vec2, 'Cin': vec2, 'E': NUM, 'sigma': NUM,
                               'schur': NUM, 'lam2p2': NUM, 'lam1h': NUM, 'h': NUM, 'temple': NUM}, {'ratio': NUM}),
    ('split',): ({'type': str, 'path': str, 'root': str, 'dirn': int, 'knew': int}, {}),
    ('esc',): ({'type': str, 'path': str, 'root': str, 'knew': int}, {}),
    ('root_done',): ({'type': str, 'root': str, 'ok': bool, 'box': dict, 'time': NUM, 'counts': dict}, {}),
    # work splitting (pool.py): a pending sub-box handed to the queue as an independent item, and the end of an item
    ('handoff',): ({'type': str, 'root': str, 'path': str, 'p0': vec4, 'r': vec4, 'depth': int, 'k': int}, {}),
    ('item_done',): ({'type': str, 'root': str, 'path': str, 'ok': bool, 'counts': dict, 'handed_off': int}, {}),
}
def ok_type(v, ty):
    if callable(ty) and not isinstance(ty, type) and not isinstance(ty, tuple): return bool(ty(v))
    if isinstance(v, bool) and ty in (int, NUM): return False
    return isinstance(v, ty)
def check_schema(d):
    t = d.get('type')
    if t == 'leaf':
        st = d.get('status')
        if st == 'discard': key = ('leaf', 'discard')
        elif st == 'cert': key = ('leaf', 'cert', d.get('method'))
        else: return f'leaf status {st!r} not whitelisted'
    elif t in ('split', 'esc', 'root_done', 'handoff', 'item_done'): key = (t,)
    else: return f'record type {t!r} not whitelisted'
    if key not in SCHEMA: return f'no schema for {key}'
    req, opt = SCHEMA[key]
    for k, ty in req.items():
        if k not in d: return f'missing key {k}'
        if not ok_type(d[k], ty): return f'bad type for {k}'
    for k in d:
        if k not in req and k not in opt: return f'unknown key {k}'
        if k in opt and not ok_type(d[k], opt[k]): return f'bad type for optional {k}'
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('log'); ap.add_argument('family'); ap.add_argument('manifest')
    ap.add_argument('--shard', default='0/1'); ap.add_argument('--out', default=None); ap.add_argument('--no-recert', action='store_true')
    ap.add_argument('--partial', action='store_true', help='TEST ONLY: check only the roots present in the log (pilots); a partial replay never counts as a proof')
    a = ap.parse_args()
    import fam5, run_fam
    from grid import dyadic_root
    if a.family not in fam5.FAMILIES: print('unknown family'); sys.exit(1)
    fam = fam5.FAMILIES[a.family]()
    si, sn = [int(x) for x in a.shard.split('/')]
    import fpenv, replay_common as rc
    sc = rc.start_checks()                                   # flags, FP environment, pinned chart caches (before any load)
    if sc: print('START CHECKS FAILED:', sc); sys.exit(1)
    if not (0 <= si < sn): print('bad shard', a.shard); sys.exit(1)
    msnap = rc.snapshot(a.manifest); lsnap = rc.snapshot(a.log)      # the replay parses exactly these bytes
    vh, vdet = binding('replay3.py', msnap['sha256'], a.family); bad = []
    def fail(msg):
        bad.append(msg)
        if len(bad) <= 50: print('FAIL:', msg, flush=True)
    try: man = json.loads(msnap['bytes'])
    except Exception: man = None
    if not isinstance(man, dict): man = dict(family=None, grid=None, roots=[], skipped=[], _bad=True)
    if set(man) != {'family', 'grid', 'roots', 'skipped'}: fail('manifest: keys are not exactly family/grid/roots/skipped'); man = dict(family=None, grid=None, roots=[], skipped=[])
    if man.get('family') != a.family: fail('manifest family mismatch')
    if json.dumps(man.get('grid')) != json.dumps(json.loads(json.dumps(run_fam.GRIDS[a.family]))): fail('manifest grid differs from run_fam.GRIDS')
    def cell_ok(c, n):
        return (isinstance(c, list) and len(c) == n and isinstance(c[0], str) and vec4(c[1]) and vec4(c[2])
                and (n == 3 or isinstance(c[3], str)))
    if not (isinstance(man['roots'], list) and all(cell_ok(c, 3) for c in man['roots'])): fail('manifest: malformed root entry'); man['roots'] = []
    if not (isinstance(man['skipped'], list) and all(cell_ok(c, 4) for c in man['skipped'])): fail('manifest: malformed skipped entry'); man['skipped'] = []
    # 0. grid / manifest consistency and skipped cells
    axes = [run_fam.axis_cells(*g) for g in run_fam.GRIDS[a.family]]
    import itertools
    expect = {}
    for idx in itertools.product(*[range(len(x)) for x in axes]):
        rid = a.family + '_' + '_'.join(map(str, idx))
        expect[rid] = ([axes[i][idx[i]][0] for i in range(4)], [axes[i][idx[i]][1] for i in range(4)])
    mroots = {}; seen = set()
    def is_cell(rid, p0, r):
        return (rid in expect and [float(x) for x in expect[rid][0]] == [float(x) for x in p0]   # RIGOR: exact: float() of floats; exact equality
                and [float(x) for x in expect[rid][1]] == [float(x) for x in r])   # RIGOR: exact: float() of floats; exact equality
    for rid, p0, r in man['roots']:
        if rid in seen: fail(f'duplicate manifest cell {rid}')
        seen.add(rid); mroots[rid] = (p0, r)
        if not is_cell(rid, p0, r): fail(f'manifest root {rid} is not the grid cell (p0 and r must both match)')
        try: exact_box(p0, r)
        except InexactBox as e: fail(f'manifest root {rid}: {e}')
    for rid, p0, r, why in man['skipped']:
        if rid in seen: fail(f'duplicate manifest cell {rid}')
        seen.add(rid)
        if not is_cell(rid, p0, r): fail(f'skipped {rid} is not the grid cell (p0 and r must both match)'); continue
        try:
            exact_box(p0, r)
            if not fam.check(why, np.array(p0, float), np.array(r, float)): fail(f'skipped cell {rid}: predicate {why!r} does not hold')   # RIGOR: exact: float64 arrays of Python floats (no rounding)
        except InexactBox as e: fail(f'skipped cell {rid}: {e}')
    if seen != set(expect): fail(f'manifest cells {len(seen)} != grid cells {len(expect)}')
    if a.partial:
        present = set()
        for line in lsnap['bytes'].decode().split('\n'):
            try: present.add(json.loads(line).get('root'))
            except Exception: pass
        mroots = {k: v for k, v in mroots.items() if k in present}
        print('PARTIAL replay (test only):', len(mroots), 'roots', flush=True)
    # 1-2. log parsing
    leaves = {}; seen_paths = set(); done = {}; nlines = 0; handoffs = {}; items_done = {}
    if True:
        for ln, line in enumerate(lsnap['bytes'].decode().split('\n')[:-1] if lsnap['bytes'].endswith(b'\n') else lsnap['bytes'].decode().split('\n'), 1):
            nlines += 1
            s = line.strip()
            if not s: fail(f'line {ln}: empty'); continue
            try: d = json.loads(s)
            except Exception: fail(f'line {ln}: malformed JSON'); continue
            if not isinstance(d, dict): fail(f'line {ln}: not an object'); continue
            why = check_schema(d)
            if why: fail(f'line {ln}: {why}'); continue
            root = d['root']
            if root not in mroots: fail(f'line {ln}: root {root} not in manifest roots'); continue
            if d['type'] == 'root_done':
                if root in done: fail(f'line {ln}: duplicate root_done {root}')
                done[root] = d; continue
            if d['type'] in ('handoff', 'item_done'):
                if not d['path'].startswith(root) or any(c not in '01k' for c in d['path'][len(root):]): fail(f'line {ln}: bad path')
                bag = handoffs if d['type'] == 'handoff' else items_done
                if (root, d['path']) in bag: fail(f'line {ln}: duplicate {d["type"]} {d["path"]}')
                bag[(root, d['path'])] = d; continue
            pk = (root, d['path'])
            if pk in seen_paths: fail(f'line {ln}: duplicate path')
            seen_paths.add(pk)
            if not d['path'].startswith(root) or any(c not in '01k' for c in d['path'][len(root):]): fail(f'line {ln}: bad path')
            if d['type'] == 'leaf':
                try: exact_box(d['p0'], d['r'])
                except InexactBox as e: fail(f'line {ln}: leaf box {e}'); continue
                d['_line'] = ln; d['_sha'] = rec_hash({k: v for k, v in d.items() if not k.startswith('_')})
                leaves.setdefault(root, []).append(d)
    for (rt, pth), d in items_done.items():
        if d['ok'] is not True: fail(f'root {rt} item {pth}: ok = {d["ok"]}')
    for (rt, pth) in handoffs:
        if (rt, pth) not in items_done: fail(f'root {rt}: handed-off item {pth} not done (pending)')
    for rid, (p0, r) in mroots.items():
        sp0, sr = dyadic_root(p0, r)
        if rid in done:
            rd = done[rid]
            if rd['ok'] is not True: fail(f'root {rid}: ok = {rd["ok"]}')
            if [float(x) for x in rd['box'].get('p0', [])] != [float(x) for x in sp0] or [float(x) for x in rd['box'].get('r', [])] != [float(x) for x in sr]:   # RIGOR: exact: float() of JSON floats; exact equality
                fail(f'root {rid}: recorded box differs from the snapped manifest box')
        elif (rid, rid) not in items_done:
            fail(f'root {rid}: pending (no root_done and no item_done of the root item)'); continue
        ok, why = tiles(dict(p0=sp0, r=sr), leaves.get(rid, []))
        if not ok: fail(f'root {rid}: tiling: {why}')
    # 4. re-verification
    from bo2 import bo_family
    from p2fix import p2_cert
    out = open(a.out, 'x') if a.out else None                # exclusive creation: never overwrites an output
    gcs = {}
    def gc_for(k):
        from gcert import GCert
        if k not in gcs: gcs[k] = GCert(run_fam.fem_problem(a.family)(k))
        return gcs[k]
    nchk = nok = 0; t0 = time.time(); shas = []
    for j, rid in enumerate(sorted(mroots)):
        if j % sn != si: continue
        for L in leaves.get(rid, []):
            p0 = np.array(L['p0'], float); r = np.array(L['r'], float)   # RIGOR: exact: float64 arrays of Python floats (no rounding)
            res_ok = False; info = ''
            try:
                if L['status'] == 'discard':
                    res_ok = fam.check(L['reason'], p0, r); info = L['reason']
                elif a.no_recert:
                    res_ok = True; info = 'not recertified (--no-recert)'
                elif L['method'] == 'BO2':
                    from bo2 import arm_variant
                    armA, armC = arm_variant(*fam.arms(p0, r), L.get('arms', 'polar'))
                    E, _ = bo_family(armA, armC, L['theta'], L['K'])
                    inn = fam.inners(p0, r)
                    if not (0 <= L['inner'] < len(inn)): raise ValueError('inner index out of range')
                    Ai, Ci = inn[L['inner']]
                    if (list(Ai) if Ai is not None else None) != L['Ain'] or list(Ci) != L['Cin']: raise ValueError('inner dart differs')
                    if not E > 0: raise ValueError('E <= 0')
                    res = p2_cert(Ai, Ci, E, L['k'])
                    res_ok = bool(res.get('ok')); info = f'E {E:.6g} ' + ('' if res_ok else str(res.get('why')))
                elif L['method'] == 'INC':
                    import inc
                    if fam.name != 'PP' or not getattr(fam, 'inc', False): raise ValueError('INC certificate outside family PP')
                    if not inc.verify_inclusion(p0, r, L['Aout'], L['Cout'], L['Ain'], L['Cin']): raise ValueError('inclusion check fails')
                    if not L['E'] > 0: raise ValueError('E <= 0')
                    E2, inf2 = inc.l1_point(L['Aout'], L['Cout'], L['k'])
                    if E2 is None or not L['E'] <= E2: raise ValueError(f'L1 at the outer dart not confirmed: {E2} {inf2.get("why")}')
                    res = p2_cert(L['Ain'], L['Cin'], L['E'], L['kp'])
                    res_ok = bool(res.get('ok')); info = f'E {L["E"]:.6g} (fresh {E2:.6g}) ' + ('' if res_ok else str(res.get('why')))
                elif L['method'] == 'FEM':
                    if not fam.fem: raise ValueError('FEM certificate in a family without FEM')
                    from gsigma import sigma_box
                    if hasattr(fam, 'fem_setup'):
                        st = fam.fem_setup(p0, r, L['k'])
                        if st is None: raise ValueError('FEM chart setup fails (validity not verified)')
                        pb, pc_, rc_ = st
                        if [float(x) for x in pc_] != L['pc'] or [float(x) for x in rc_] != L['rc']: raise ValueError('FARL box differs')   # RIGOR: exact: float() of floats; exact equality
                        from gcert import GCert
                        if ('setup', id(pb)) not in gcs: gcs[('setup', id(pb))] = (pb, GCert(pb))   # keep pb alive: ids stay unique
                        gc = gcs[('setup', id(pb))][1]; pc, rcb = np.array(pc_), np.array(rc_)   # RIGOR: exact: float64 arrays of Python floats (no rounding)
                    else:
                        cb = fam.cart_box(p0, r)
                        if cb is None or [float(x) for x in cb[0]] != L['pc'] or [float(x) for x in cb[1]] != L['rc']:   # RIGOR: exact: float() of floats; exact equality
                            raise ValueError('Cartesian enclosure differs')
                        pc, rcb = np.array(cb[0]), np.array(cb[1])   # RIGOR: exact: float64 arrays of Python floats (no rounding)
                        gc = gc_for(L['k'])
                    sol = gc.solve(pc, rcb)
                    sc, _ = sigma_box(gc, sol, pc, rcb)
                    res = gc.certify(pc, rcb, sc, sol=sol) if sc is not None else dict(ok=False, why='sigma')
                    res_ok = bool(res.get('ok')); info = '' if res_ok else str(res.get('why'))
            except Exception as e:
                res_ok = False; info = 'exception: ' + repr(e)[:200]
            try: fpenv.check()
            except fpenv.FPEnvError as e: res_ok = False; info = str(e)[:200]
            nchk += 1; nok += res_ok; shas.append(L['_sha'])
            if not res_ok: fail(f'root {rid} path {L["path"]} line {L["_line"]}: not re-verified: {info}')
            if out:
                out.write(json.dumps(dict(root=rid, path=L['path'], line=L['_line'], status=L['status'], method=L.get('method'),
                                          ok=res_ok, info=info, record_sha256=L['_sha'], binding_sha256=vh, family=a.family)) + '\n')
    summ = dict(log=a.log, family=a.family, shard=a.shard, partial=bool(a.partial), recert=not a.no_recert, lines=nlines, roots=len(mroots), skipped=len(man['skipped']),
                leaves_checked=nchk, ok=nok, failures=len(bad), binding_sha256=vh, manifest_sha256=vdet.get('manifest_sha256'),
                environment=vdet['environment'], leaves_sha256=hashlib.sha256(''.join(sorted(shas)).encode()).hexdigest(),
                log_sha256=lsnap['sha256'], time_s=round(time.time() - t0, 1))   # RIGOR: info: timing
    for w in rc.end_checks('replay3.py', msnap['sha256'], a.family, vh, [msnap, lsnap], vdet): fail('end of run: ' + w)
    summ['failures'] = len(bad); summ['late_loaded'] = dict(rc.LATE); summ['niceness'] = os.nice(0)   # (recorded, not bound)
    summ['binding_sources'] = vdet['sources']; summ['binding_refcache'] = vdet['refcache']
    print(json.dumps(summ), flush=True)
    if out: out.write(json.dumps(dict(summary=summ)) + '\n'); out.close()
    sys.exit(0 if not bad else 1)

if __name__ == '__main__':
    main()
