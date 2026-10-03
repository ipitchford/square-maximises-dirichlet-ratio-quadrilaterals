"""Strict validation of the outputs of a sharded replay (adversarial reviews round 2 (SERIOUS 6, 7) and round 3 (items 4, 7)).
validate(outputs, log_snap, manifest_snap, family, binding_hex) checks EVERY record of the shard outputs, given as BYTE SNAPSHOTS
(name, bytes) taken once by the caller (the driver hashes and summarises exactly the same bytes), i = 0..n-1:
  * every line parses; every line but the last is a leaf receipt with exactly the keys
    {root, path, line, status, method, ok, info, record_sha256, binding_sha256, family}, ok = true, the expected binding and family;
  * the last line is a summary for shard i/n with partial = false, recert = true, failures = 0, leaves_checked = ok = the number of
    receipts in the file, the expected family, binding, log sha256 and manifest sha256, and leaves_sha256 = the sha256 of the
    file's sorted receipt record hashes;
  * across the shards: no duplicate receipt (log line), and the multiset of receipt record hashes and the set of receipt lines
    equal those RECOMPUTED from the bound log bytes (every leaf record of the log, hashed as by the replays); every receipt's
    root, path, status and method agree with its log record (info must be a string).
Returns a list of problems (empty = valid).  The proof driver (replay_all.py) always replays fresh and then applies validate()
to its own outputs; reuse of an older replay (this script's command line) is for development only and never yields a proof.
Usage: replay_ok.py TAG LOG SCRIPT MANIFEST FAMILY [--dir DIR]   (exit 0 iff a complete valid replay of LOG exists in DIR)"""
import sys, os, glob, json, hashlib
HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
RECEIPT_KEYS = {'root', 'path', 'line', 'status', 'method', 'ok', 'info', 'record_sha256', 'binding_sha256', 'family'}

def log_leaf_hashes(log_bytes):
    """{line number: (record sha256, root, path, status, method)} of every leaf record of the log, hashed exactly as by replay2 /
    replay3; method as the replays report it (legacy cert records without a method key are FEM certificates)"""
    from replay2 import rec_hash
    out = {}
    text = log_bytes.decode()
    lines = text.split('\n')
    if text.endswith('\n'): lines = lines[:-1]
    for ln, line in enumerate(lines, 1):
        s = line.strip()
        if not s: continue
        try: d = json.loads(s)
        except Exception: continue
        if isinstance(d, dict) and d.get('type') == 'leaf':
            out[ln] = (rec_hash({k: v for k, v in d.items() if not k.startswith('_')}), d.get('root'), d.get('path'), d.get('status'),
                       d.get('method', 'FEM' if d.get('status') == 'cert' else None))
    return out

def validate(outputs, log_snap, manifest_snap, family, binding_hex=None, expect_shards=None, recert=True):
    """outputs: list of (name, bytes) snapshots of the shard output files; returns a list of problems (empty = valid).
    binding_hex None: the binding is taken from the first summary and every shard and receipt must carry the same one (the caller
    checks that binding's recorded sources, chart caches and environment, which a launching process cannot recompute itself)"""
    bad = []
    if not outputs: return ['no replay output files']
    if binding_hex is None:
        try: binding_hex = json.loads(outputs[0][1].decode().rstrip('\n').split('\n')[-1])['summary']['binding_sha256']
        except Exception: return ['first shard has no readable summary']
    seen_lines = {}; meta = {}; n_expected = None; idx = set()
    for f, rawb in outputs:
        try: raw = rawb.decode()
        except Exception as e: bad.append(f'{f}: not UTF-8 ({e})'); continue
        if not raw.endswith('\n'): bad.append(f'{f}: torn last line'); continue
        lines = raw.split('\n')[:-1]
        if not lines: bad.append(f'{f}: empty'); continue
        try: summ = json.loads(lines[-1])['summary']
        except Exception: bad.append(f'{f}: last line is not a summary'); continue
        try: i, n = (int(x) for x in str(summ.get('shard')).split('/'))
        except Exception: bad.append(f'{f}: bad shard field'); continue
        if n_expected is None: n_expected = n
        if n != n_expected or not (0 <= i < n) or i in idx: bad.append(f'{f}: shard {i}/{n} inconsistent or duplicate'); continue
        idx.add(i)
        recs = []
        for k, line in enumerate(lines[:-1], 1):
            try: d = json.loads(line)
            except Exception: bad.append(f'{f}:{k}: not JSON'); continue
            if not isinstance(d, dict) or set(d) != RECEIPT_KEYS: bad.append(f'{f}:{k}: receipt keys'); continue
            if d['ok'] is not True: bad.append(f'{f}:{k}: receipt not ok'); continue
            if d['binding_sha256'] != binding_hex: bad.append(f'{f}:{k}: receipt binding differs'); continue
            if d['family'] != family: bad.append(f'{f}:{k}: receipt family {d["family"]!r} != {family!r}'); continue
            if not isinstance(d['line'], int): bad.append(f'{f}:{k}: receipt line'); continue
            if d['line'] in seen_lines: bad.append(f'{f}:{k}: duplicate receipt for log line {d["line"]}'); continue
            if not isinstance(d['info'], str): bad.append(f'{f}:{k}: receipt info is not a string'); continue
            seen_lines[d['line']] = d['record_sha256']; meta[d['line']] = (d['root'], d['path'], d['status'], d['method'])
            recs.append(d['record_sha256'])
        want = dict(partial=False, recert=recert, failures=0, family=family, binding_sha256=binding_hex, log_sha256=log_snap['sha256'],
                    manifest_sha256=manifest_snap['sha256'])
        for k, v in want.items():
            if summ.get(k) != v: bad.append(f'{f}: summary {k} = {summ.get(k)!r}, expected {v!r}')
        if not (summ.get('leaves_checked') == summ.get('ok') == len(recs)): bad.append(f'{f}: summary counts {summ.get("leaves_checked")}/{summ.get("ok")} vs {len(recs)} receipts')
        if summ.get('leaves_sha256') != hashlib.sha256(''.join(sorted(recs)).encode()).hexdigest(): bad.append(f'{f}: leaves_sha256 does not match its receipts')
    if n_expected is not None and idx != set(range(n_expected)): bad.append(f'shards present {sorted(idx)} != 0..{n_expected - 1}')
    if expect_shards is not None and n_expected != expect_shards: bad.append(f'{n_expected} shards, expected {expect_shards}')
    exp = log_leaf_hashes(log_snap['bytes'])
    if set(seen_lines) != set(exp): bad.append(f'receipts cover {len(seen_lines)} log lines, the log has {len(exp)} leaf records (sets differ)')
    else:
        if any(seen_lines[ln] != exp[ln][0] for ln in exp): bad.append('a receipt record hash differs from the recomputed hash of its log line')
        wrong = [ln for ln in exp if meta[ln] != exp[ln][1:]]
        if wrong: bad.append(f'{len(wrong)} receipts disagree with their log records in root / path / status / method (first: line {wrong[0]})')
    return bad

if __name__ == '__main__':
    import replay_common as rc
    from binding import binding
    tag, log, script, manifest, family = sys.argv[1:6]
    d = sys.argv[sys.argv.index('--dir') + 1] if '--dir' in sys.argv else os.path.join(HERE, 'replays')
    ls = rc.snapshot(os.path.join(HERE, log) if not os.path.isabs(log) else log)
    ms = rc.snapshot(os.path.join(HERE, manifest) if not os.path.isabs(manifest) else manifest)
    files = sorted(glob.glob(os.path.join(d, f'{tag}.[0-9]*.jsonl')))
    outs = [(f, open(f, 'rb').read()) for f in files]
    # the binding is taken from the summaries (its runtime part belongs to the replaying process, which this one cannot recompute);
    # its recorded source closure and chart caches must equal the current files
    problems = validate(outs, ls, ms, family, None)
    if not problems:
        det = binding(script, ms['sha256'], family)[1]
        summ = json.loads(outs[0][1].decode().rstrip('\n').split('\n')[-1])['summary']
        if summ.get('binding_sources') != det['sources']: problems.append('the replay was made with other source files than the current ones')
        if summ.get('binding_refcache') != det['refcache']: problems.append('the replay was made with other chart caches than the current ones')
    for p in problems[:20]: print('INVALID:', p)
    if problems: sys.exit(1)
    print(f'valid: {tag} ({len(files)} shards, every receipt checked against the bound log)')
