"""Pre-replay completeness gate of a re-certification directory (default certs4): a fast check (seconds, not hours) that a replay
may start.  For every family F (log prefix f = F.lower()) of the list:
  * the worker logs are exactly the files DIR/<f>.w<ID>.jsonl (merge_pool.py merges every DIR/<f>.*.jsonl but the merged log, so
    any other such file is reported), with worker IDs 1..N, N = DIR/RECERT_N (recorded by recert_all.sh at the first start;
    required unless require_lock is False, which is for test directories only);
  * every line is a complete JSON object with a string root (no torn or garbage line) of a known record type, every leaf has
    status cert, discard or FAIL;
  * every record names a manifest root; every manifest root has exactly ONE root_done record, and it has ok = true;
  * no FAIL leaf -- every FAIL leaf is listed, so that it can be decided on before any replay time is spent;
  * no (root, path) tree record (leaf / split / esc) occurs twice; every handoff has its item_done; no handoff or item_done occurs
    twice; every item_done has ok = true;
  * every claim DIR/claims/<F>/<root> names a manifest root and a worker 1..N; every manifest root is claimed; every claimed root is
    done (a claim without its root_done is an orphan: the root would never be finished).
check(dir, families, require_lock) -> dict(problems, fails, inputs {file name: sha256 of the bytes checked}, roots {F: count}, N).
The proof driver (replay_all.py) runs check() before any gate, merge or replay and aborts on any problem or FAIL leaf; the merge
inputs must then have exactly the checked sha256.  recert_all.sh runs it when its workers have finished.
Usage: certs_complete.py DIR [FAMILY ...] [--no-lock]   (exit 0 iff complete; DIR absolute or relative to the package)"""
import sys, os, json, hashlib, re
HERE = '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts'
ALL = ('PP', 'SP_spk', 'SP_col', 'FP', 'PF', 'FF', 'SF_spk', 'SF_col')
TYPES = ('leaf', 'split', 'esc', 'root_done', 'handoff', 'item_done')

def read_n(d):
    """N recorded in DIR/RECERT_N, or None if absent / malformed"""
    try:
        with open(os.path.join(d, 'RECERT_N')) as fh: s = fh.read().strip()
    except OSError: return None
    return int(s) if re.fullmatch(r'[1-9][0-9]?', s) else None

def check(d, families=ALL, require_lock=True):
    d = d if os.path.isabs(d) else os.path.join(HERE, d)
    problems = []; fails = []; inputs = {}; roots = {}
    if not os.path.isdir(d): return dict(problems=[f'{d}: no such directory'], fails=[], inputs={}, roots={}, N=None)
    N = read_n(d)
    if N is None and require_lock: problems.append('RECERT_N missing or malformed (recert_all.sh records the worker count at the first start)')
    listing = sorted(os.listdir(d))
    for F in families:
        if F not in ALL: problems.append(f'unknown family {F!r}'); continue
        f = F.lower()
        try:
            with open(os.path.join(HERE, 'manifests', f + '.json'), 'rb') as fh: man = json.loads(fh.read())
            mroots = {r[0] for r in man['roots']}
        except Exception as e: problems.append(f'{F}: manifest unreadable: {e!r}'); continue
        roots[F] = len(mroots)
        files = [x for x in listing if x.startswith(f + '.') and x.endswith('.jsonl') and x != f + '.merged.jsonl']
        done = {}; seen = {}; hand = {}; items = {}
        for fn in files:
            m = re.fullmatch(re.escape(f) + r'\.w([1-9][0-9]?)\.jsonl', fn)
            if not m: problems.append(f'{fn}: not a worker log name (merge_pool.py would merge it)'); continue
            if N is not None and int(m.group(1)) > N: problems.append(f'{fn}: worker ID above N = {N}')
            with open(os.path.join(d, fn), 'rb') as fh: raw = fh.read()            # the bytes checked (hashed for the merge)
            inputs[fn] = hashlib.sha256(raw).hexdigest()
            if raw and not raw.endswith(b'\n'): problems.append(f'{fn}: torn last line')
            for i, line in enumerate(raw.split(b'\n'), 1):
                if not line.strip(): continue
                try: rec = json.loads(line)
                except Exception: problems.append(f'{fn}:{i}: not JSON (torn or garbage line)'); continue
                if not isinstance(rec, dict) or not isinstance(rec.get('root'), str): problems.append(f'{fn}:{i}: not a record with a root'); continue
                rid = rec['root']; t = rec.get('type'); where = f'{fn}:{i}'
                if t not in TYPES: problems.append(f'{where}: unknown record type {t!r}'); continue
                if rid not in mroots: problems.append(f'{where}: root {rid} is not a manifest root'); continue
                if t == 'root_done': done.setdefault(rid, []).append((where, rec.get('ok'))); continue
                key = (rid, rec.get('path'))
                if t in ('handoff', 'item_done'):
                    bag = hand if t == 'handoff' else items
                    if key in bag: problems.append(f'{where}: duplicate {t} {key[1]} (also {bag[key][0]})'); continue
                    bag[key] = (where, rec.get('ok'))
                    if t == 'item_done' and rec.get('ok') is not True: problems.append(f'{where}: item_done {key[1]} has ok = {rec.get("ok")!r}')
                    continue
                if key in seen: problems.append(f'{where}: duplicate tree record {rid} {key[1]} (also {seen[key]})'); continue
                seen[key] = where
                if t == 'leaf':
                    st = rec.get('status')
                    if st == 'FAIL':
                        fails.append(dict(family=F, root=rid, path=key[1], p0=rec.get('p0'), r=rec.get('r'), why=str(rec.get('res'))[:300], at=where))
                    elif st not in ('cert', 'discard'): problems.append(f'{where}: leaf status {st!r}')
        for rid in sorted(mroots):
            ds = done.get(rid, [])
            if len(ds) != 1: problems.append(f'{F} root {rid}: {len(ds)} root_done records' + (f' ({", ".join(w for w, _ in ds)})' if ds else ''))
            elif ds[0][1] is not True: problems.append(f'{F} root {rid}: root_done has ok = {ds[0][1]!r} ({ds[0][0]})')
        for key in sorted(set(hand) - set(items)): problems.append(f'{F} root {key[0]}: handoff {key[1]} has no item_done ({hand[key][0]})')
        cdir = os.path.join(d, 'claims', F); claims = {}
        for c in (sorted(os.listdir(cdir)) if os.path.isdir(cdir) else []):
            if c.startswith('.'): continue                              # temporary files of the atomic claim protocol
            try:
                with open(os.path.join(cdir, c)) as fh: owner = fh.read().strip()
            except OSError as e: problems.append(f'{F} claim {c}: unreadable ({e})'); continue
            if c not in mroots: problems.append(f'{F} claim {c}: not a manifest root')
            if not re.fullmatch(r'[1-9][0-9]?', owner) or (N is not None and int(owner) > N): problems.append(f'{F} claim {c}: owner {owner!r} is not a worker 1..{N}')
            claims[c] = owner
        for rid in sorted(mroots):
            if rid not in claims: problems.append(f'{F} root {rid}: never claimed')
            elif rid not in done: problems.append(f'{F} root {rid}: claimed by worker {claims[rid]} but without root_done (orphaned or unfinished)')
    return dict(problems=problems, fails=fails, inputs=inputs, roots=roots, N=N)

def main():
    args = [x for x in sys.argv[1:] if not x.startswith('--')]
    if not args: print(__doc__); sys.exit(2)
    fams = tuple(args[1:]) or ALL
    res = check(args[0], fams, require_lock='--no-lock' not in sys.argv)
    for p in res['problems'][:200]: print('INCOMPLETE:', p)
    if len(res['problems']) > 200: print(f'INCOMPLETE: ... {len(res["problems"]) - 200} more')
    for x in res['fails']: print('FAIL LEAF:', json.dumps(x))
    ok = not res['problems'] and not res['fails']
    print(json.dumps(dict(complete=ok, problems=len(res['problems']), fail_leaves=len(res['fails']), roots=res['roots'], N=res['N'],
                          files=len(res['inputs']))))
    sys.exit(0 if ok else 1)

if __name__ == '__main__':
    main()
