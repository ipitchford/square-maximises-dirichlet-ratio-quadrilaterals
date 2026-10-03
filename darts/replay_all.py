"""Aggregate fail-closed proof replay of the dart package (adversarial reviews round 2 (SERIOUS 2, 5, 7) and round 3).
Threat model (README §5): accidental corruption, stale or changed files, concurrent agents editing files, environment drift.
Deliberate tampering by the owner of the machine is out of scope.

Launch ONLY through replay_all.sh, which runs this file as
    env -i PATH=/usr/bin:/bin HOME=... <venv python> -I -S -B -X pycache_prefix=<new empty dir> <darts>/bootstrap.py <darts>/replay_all.py
Usage: replay_all.py --pin PINFILE [--shards N] [--mode proof|quick|test] [--certs DIR] [--run-id ID] [--only TAGS --gates all|none]
  --shards N   canonical integer 1..64 (default 6, the queue's CPU reservation): exactly N replay processes per log
  --mode       proof (default): every certificate recomputed; the ONLY mode that can print PASS and exit 0.
               quick: structure only, never PASS (exit 3 when structurally OK); test: machinery test on a subset (--only, --gates),
               never PASS (exit 3 when OK).  Exit 1: FAIL; exit 2: usage / launch error.
What a run does (every step recorded in the receipt):
  * reads the pin list ONCE (its sha256 is the hash of exactly the parsed bytes), writes a private read-only copy into the new
    run directory (children read that copy through DARTS_PIN), and checks: every pinned file has its pinned sha256, every
    transitive decision source of every stage is pinned, the chart cache directory equals the pinned refcache/* inventory;
  * checks its own interpreter flags (-I -S -B, empty pycache prefix, no -O) and the FP environment, lowers its own priority
    in-process to the TARGET niceness 10 (os.nice(max(0, 10 - current)): under a launcher that already applied nice 10, such as
    research-run, the run stays at 10, not 20; the children inherit it and report it in their summaries), and launches every
    child as <venv python, absolute> -I -S -B
    -X pycache_prefix=<run>/pycache bootstrap.py <script> with a sanitised environment (PATH=/usr/bin:/bin, HOME, TMPDIR, the BLAS
    thread variables, DARTS_PIN; nothing else);
  * completeness (seconds, before any gate, merge or replay): certs_complete.check() of the re-certification directory for the
    selected cover5 families -- every manifest root exactly one root_done with ok = true, no FAIL leaf (every FAIL leaf is listed
    on stdout and in the receipt), no torn / garbage / duplicate record, every handoff with its item_done, every claim of a
    worker 1..N (RECERT_N) and done; any problem ABORTS the run before any replay starts; the merge inputs must then have exactly
    the sha256 of the bytes this check read;
  * gates, merges and replays: gate outputs are captured in memory, hashed and written exclusively; a merged log is snapshotted
    once after its merge and its sha256 recorded; before every replay the chart caches are re-verified against the pin, the log
    and manifest are snapshotted once and must equal their pin entry (old logs, manifests) or the recorded merge output (family
    logs); the children create their outputs exclusively; after the children exit each output is snapshotted ONCE and exactly
    those bytes are validated (replay_ok.validate, incl. log-record metadata), summarised and hashed; a missing or unreadable
    summary FAILS the run;
  * end of run: the original pin file and the private copy, every pinned file, the chart caches, every consumed log and manifest
    and every shard output must still have their recorded sha256; a FRESH runtime inventory must equal the start one; the
    driver's loaded modules must come from unchanged bound files; the pycache directories must be empty;
  * receipt replays/receipts/<ID>.jsonl: created exclusively; one canonical JSON record per line with seq (strictly 0, 1, 2, ...),
    prev_sha256 (the previous record's self_sha256), time, kind (start, complete, gate, merge, replay, end, abort), data and self_sha256 (the
    sha256 of the canonical record without it); every record written completely (write loop) and fsynced, the directory fsynced;
    after the final record the receipt is closed, made read-only, re-read and fully validated (read_receipt) before the final
    status is printed together with the final digest (the self_sha256 of the end record), which is to be anchored externally."""
import argparse, hashlib, json, os, re, subprocess, sys, time, secrets
HERE = os.path.dirname(os.path.realpath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
GATES = [['lemma_sector.py'], ['lemma_longarm.py'], ['check_charts.py', '8'], ['rigor_lint.py'], ['test_containment.py'],
         ['hull_to_box.py', '56', '64', '0.0339', '0.12', '0.5', '0.875'], ['check_cbox_cover.py'], ['accounting.py']]
LEGACY = [('certs/crit_collar.jsonl', 'COLLAR', 'manifests/crit_collar.json', 'crit_collar'),
          ('certs/crit_spike_X1.jsonl', 'SPIKE_X1', 'manifests/crit_spike_X1.json', 'crit_spike_X1'),
          ('certs/main_cbox.jsonl', 'MAIN_CBOX', 'manifests/main_cbox.json', 'main_cbox'),
          ('certs/shallow_crit.jsonl', 'BOX', 'manifests/shallow_crit.json', 'shallow_crit'),
          ('certs2/kf_collar.jsonl', 'KF_COLLAR', 'manifests/kf_collar.json', 'kf_collar'),
          ('certs2/kf_spike.jsonl', 'KF_SPIKE', 'manifests/kf_spike.json', 'kf_spike')]
FAMILIES = ['PP', 'SP_col', 'SP_spk', 'PF', 'FP', 'FF', 'SF_col', 'SF_spk']
ENTRY_POINTS = ['replay_all.py', 'bootstrap.py', 'replay2.py', 'replay3.py', 'replay_ok.py', 'merge_pool.py'] + [g[0] for g in GATES]
RECORD_KEYS = {'seq', 'prev_sha256', 'time', 'kind', 'data', 'self_sha256'}
KINDS = {'start', 'complete', 'gate', 'merge', 'replay', 'end', 'abort'}
EXIT_PASS, EXIT_FAIL, EXIT_USAGE, EXIT_NOT_PROOF = 0, 1, 2, 3

def canon(obj): return json.dumps(obj, sort_keys=True, separators=(',', ':'), ensure_ascii=True)
def plain(obj): return json.loads(json.dumps(obj, default=str))
def sha_bytes(b): return hashlib.sha256(b).hexdigest()
def sha_file(p):
    with open(p, 'rb') as fh: return sha_bytes(fh.read())
def write_all(fd, data):
    view = memoryview(data)
    while view:
        n = os.write(fd, view)
        if n <= 0: raise OSError('write returned 0')
        view = view[n:]
def fsync_dir(d):
    fd = os.open(d, os.O_RDONLY)
    try: os.fsync(fd)
    finally: os.close(fd)
def write_exclusive(path, data, mode=0o444):
    fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    try: write_all(fd, data); os.fsync(fd)
    finally: os.close(fd)
    os.chmod(path, mode); fsync_dir(os.path.dirname(path))

def final_status(mode, failures):
    """only a proof-mode run without any failure is a PASS"""
    if failures: return 'FAIL'
    return {'proof': 'PASS', 'quick': 'STRUCTURE-ONLY (NOT A PROOF)', 'test': 'MACHINERY TEST OK (NOT A PROOF)'}[mode]
def exit_code(status):
    return EXIT_PASS if status == 'PASS' else (EXIT_FAIL if status == 'FAIL' else EXIT_NOT_PROOF)

class Receipt:
    def __init__(self, path):
        self.path = path; self.dir = os.path.dirname(path)
        os.makedirs(self.dir, exist_ok=True)
        self.fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY | os.O_APPEND, 0o644)
        fsync_dir(self.dir)
        self.seq = 0; self.prev = '0'*64
    def add(self, kind, **data):
        if kind not in KINDS: raise ValueError(kind)
        rec = dict(seq=self.seq, prev_sha256=self.prev, time=time.strftime('%Y-%m-%dT%H:%M:%S'), kind=kind, data=plain(data))
        rec['self_sha256'] = sha_bytes(canon(rec).encode())
        write_all(self.fd, (canon(rec) + '\n').encode()); os.fsync(self.fd)
        self.prev = rec['self_sha256']; self.seq += 1
        return self.prev
    def close(self):
        os.fsync(self.fd); os.close(self.fd); os.chmod(self.path, 0o444); fsync_dir(self.dir)

def read_receipt(raw):
    """(records, final digest, problems) of receipt bytes: canonical lines, strict sequence, chained and self hashes, schema"""
    problems = []; recs = []; prev = '0'*64
    if not raw.endswith(b'\n'): problems.append('receipt does not end with a newline (torn)')
    lines = raw.split(b'\n')[:-1] if raw.endswith(b'\n') else raw.split(b'\n')
    for i, line in enumerate(lines):
        try: rec = json.loads(line)
        except Exception: problems.append(f'record {i}: not JSON'); break
        if not isinstance(rec, dict) or set(rec) != RECORD_KEYS: problems.append(f'record {i}: keys'); break
        if canon(rec).encode() != line: problems.append(f'record {i}: not in canonical form'); break
        if rec['seq'] != i: problems.append(f'record {i}: seq {rec["seq"]}'); break
        if rec['prev_sha256'] != prev: problems.append(f'record {i}: prev_sha256 does not chain'); break
        if rec['kind'] not in KINDS or not isinstance(rec['data'], dict): problems.append(f'record {i}: kind / data'); break
        body = {k: v for k, v in rec.items() if k != 'self_sha256'}
        if sha_bytes(canon(body).encode()) != rec['self_sha256']: problems.append(f'record {i}: self_sha256 mismatch'); break
        prev = rec['self_sha256']; recs.append(rec)
    if recs and recs[0]['kind'] != 'start': problems.append('first record is not start')
    ends = [r for r in recs if r['kind'] in ('end', 'abort')]
    if len(ends) > 1 or (ends and ends[0] is not recs[-1]): problems.append('end / abort record is not unique and last')
    return recs, prev, problems

def summaries_of(outputs):
    """(summaries, problems) of shard output snapshots [(name, bytes)]: a missing or unreadable summary is a problem (FAIL)"""
    out = []; bad = []
    for name, b in outputs:
        try:
            s = json.loads(b.decode().rstrip('\n').split('\n')[-1])['summary']
            if not isinstance(s, dict): raise ValueError('summary is not an object')
            out.append(s)
        except Exception: bad.append(f'{name}: summary missing or unreadable'); out.append(None)
    return out, bad

def read_pin_bytes(raw):
    out = {}
    for line in raw.decode().split('\n'):
        if not line: continue
        m = re.fullmatch(r'([0-9a-f]{64})  (\S.*)', line)
        if not m or m.group(2) in out: raise ValueError(f'malformed or duplicate pin line: {line!r}')
        out[m.group(2)] = m.group(1)
    return out

def check_pinned(pins, which=None):
    bad = []
    for f, h in sorted(pins.items()):
        if which is not None and not which(f): continue
        p = os.path.join(HERE, f)
        if not os.path.isfile(p) or os.path.islink(p): bad.append(f'{f}: missing or not a regular file'); continue
        if sha_file(p) != h: bad.append(f'{f}: sha256 differs from the pin')
    return bad

def check_cache_inventory(pins):
    want = {f[len('refcache/'):] for f in pins if f.startswith('refcache/')}
    have = set(os.listdir(os.path.join(HERE, 'refcache'))) if os.path.isdir(os.path.join(HERE, 'refcache')) else set()
    bad = []
    if not want: bad.append('the pin lists no chart cache files')
    if have != want: bad.append(f'chart cache inventory differs from the pin: extra {sorted(have - want)}, missing {sorted(want - have)}')
    return bad + check_pinned(pins, lambda f: f.startswith('refcache/'))

def own_flags_problems():
    f = sys.flags; bad = []
    if f.optimize: bad.append('python -O: refused')
    if not (f.isolated and f.no_site and f.dont_write_bytecode and sys.pycache_prefix): bad.append('driver not launched by replay_all.sh (-I -S -B -X pycache_prefix)')
    return bad

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pin', required=True); ap.add_argument('--shards', default='6'); ap.add_argument('--mode', default='proof')
    ap.add_argument('--certs', default='certs4'); ap.add_argument('--run-id', default=None)
    ap.add_argument('--only', default=None); ap.add_argument('--gates', default='all')
    a = ap.parse_args()
    if not re.fullmatch(r'[1-9][0-9]?', a.shards) or not (1 <= int(a.shards) <= 64): print(f'--shards must be a canonical integer in 1..64, got {a.shards!r}'); sys.exit(EXIT_USAGE)
    if a.mode not in ('proof', 'quick', 'test'): print(f'--mode must be proof, quick or test, got {a.mode!r}'); sys.exit(EXIT_USAGE)
    if a.mode != 'test' and (a.only is not None or a.gates != 'all'): print('--only / --gates are allowed in --mode test only'); sys.exit(EXIT_USAGE)
    if a.gates not in ('all', 'none'): print('--gates must be all or none'); sys.exit(EXIT_USAGE)
    fl = own_flags_problems()
    if fl: print('LAUNCH ERROR:', fl); sys.exit(EXIT_USAGE)
    only = set(a.only.split(',')) if a.only else None
    N = int(a.shards); quick = a.mode == 'quick'
    pin_path = os.path.realpath(a.pin)
    try:
        with open(pin_path, 'rb') as fh: pin_raw = fh.read()              # the ONLY read of the pin file before the end check
        pins = read_pin_bytes(pin_raw)
    except (OSError, ValueError) as e:
        print(f'pin file {a.pin!r} unreadable or malformed: {e}'); sys.exit(EXIT_USAGE)
    if not pins: print('empty pin'); sys.exit(EXIT_USAGE)
    pin_sha = sha_bytes(pin_raw)
    run_id = a.run_id or (time.strftime('%Y%m%dT%H%M%S') + '_' + secrets.token_hex(4))
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', run_id): print('bad run id'); sys.exit(EXIT_USAGE)
    base = os.path.join(HERE, 'replays') if a.mode != 'test' else os.path.join(HERE, 'replays', 'test_runs')
    run_dir = os.path.join(base, f'run_{run_id}')
    os.makedirs(run_dir, exist_ok=False)
    for d in ('pycache', 'tmp'): os.makedirs(os.path.join(run_dir, d))
    pin_copy = os.path.join(run_dir, 'PIN.sha256'); write_exclusive(pin_copy, pin_raw)
    rc_ = Receipt(os.path.join(base, 'receipts', f'{run_id}.jsonl'))
    # niceness TARGET 10 (in-process, inherited by every child): under a launcher that already applied nice 10 (research-run)
    # the run stays at 10, not 20; a higher start value is kept (an unprivileged process cannot lower it)
    os.nice(max(0, 10 - os.nice(0)))   # RIGOR: exact: integer niceness arithmetic (no numerics)
    import fpenv, replay_common as rcm
    from binding import closure, static_inventory, runtime_identity, binding, loaded_modules_problems
    from replay_ok import validate
    failures = []
    def failed(stage, msg):
        failures.append(f'{stage}: {msg}'); print(f'FAIL {stage}: {msg}', flush=True)
    # ---- start
    for b in check_pinned(pins): failed('pin', b)
    needed = set(closure(ENTRY_POINTS)) | {'replay_all.sh', 'rigor_allowlist.json'}
    for f in sorted(needed - set(pins)): failed('pin', f'decision source {f} is not pinned')
    for b in check_cache_inventory(pins): failed('chart caches', b)
    try: fpenv.check(); fp_ok = True
    except fpenv.FPEnvError as e: fp_ok = False; failed('fpenv', str(e))
    inv0 = static_inventory(); ident0 = runtime_identity(inv0)
    exe = sys.executable
    if not (os.path.isabs(exe) and os.path.isfile(os.path.join(os.path.dirname(os.path.dirname(exe)), 'pyvenv.cfg'))): failed('launch', f'interpreter {exe} is not the absolute venv interpreter')
    child_env = {'PATH': '/usr/bin:/bin', 'HOME': os.environ.get('HOME', '/var/empty'), 'TMPDIR': os.path.join(run_dir, 'tmp'),
                 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1', 'DARTS_PIN': pin_copy}
    def child(script, *args):
        return [exe, '-I', '-S', '-B', '-X', 'pycache_prefix=' + os.path.join(run_dir, 'pycache'), os.path.join(HERE, 'bootstrap.py'),
                os.path.join(HERE, script)] + list(args)
    rc_.add('start', run_id=run_id, mode=a.mode, shards=N, certs=a.certs, argv=sys.argv, pin_file=pin_path, pin_sha256=pin_sha,
            pin_count=len(pins), pin_copy=os.path.relpath(pin_copy, HERE), driver_sha256=sha_file(os.path.realpath(__file__)),
            wrapper_sha256=sha_file(os.path.join(HERE, 'replay_all.sh')), runtime_identity=ident0, child_env=child_env,
            child_command=child('<script>'), niceness=os.nice(0), failures_so_far=list(failures))
    if failures:
        rc_.add('abort', reason='start checks failed', failures=failures); rc_.close()
        print('ABORT (start checks); receipt', os.path.relpath(rc_.path, HERE)); sys.exit(EXIT_FAIL)
    # ---- completeness of the re-certification directory (seconds): no gate, merge or replay starts on incomplete certificates
    import certs_complete
    fams_sel = [F for F in FAMILIES if only is None or F.lower() in only]
    t0 = time.time()
    cc = (certs_complete.check(os.path.join(HERE, a.certs), fams_sel, require_lock=(a.mode != 'test')) if fams_sel
          else dict(N=None, roots={}, inputs={}, problems=[], fails=[]))     # (test runs of legacy logs only)
    rc_.add('complete', certs=a.certs, families=fams_sel, N=cc['N'], roots=cc['roots'], inputs=cc['inputs'], problems=len(cc['problems']),
            first_problems=cc['problems'][:100], fail_leaves=len(cc['fails']), fail_leaf_list=cc['fails'][:500],
            seconds=round(time.time() - t0, 1))   # RIGOR: info: timing
    for x in cc['fails']: print('FAIL LEAF:', json.dumps(x), flush=True)
    for pr in cc['problems'][:50]: failed('complete', pr)
    if len(cc['problems']) > 50: failed('complete', f'... {len(cc["problems"]) - 50} more completeness problems')
    if cc['fails']: failed('complete', f'{len(cc["fails"])} FAIL leaves (listed above and in the receipt)')
    if failures:
        rc_.add('abort', reason='re-certification directory incomplete', failures=failures); rc_.close()
        print('ABORT (incomplete certificates: no replay started); receipt', os.path.relpath(rc_.path, HERE)); sys.exit(EXIT_FAIL)
    # ---- gates
    for g in (GATES if a.gates == 'all' else []):
        t0 = time.time()
        r = subprocess.run(child(*g), cwd=HERE, env=child_env, capture_output=True)
        outb = r.stdout + b'\n--- stderr ---\n' + r.stderr
        out = os.path.join(run_dir, 'gate_' + g[0].replace('.py', '') + '.out'); write_exclusive(out, outb)
        rc_.add('gate', command=g, exit=r.returncode, output=os.path.relpath(out, HERE), output_sha256=sha_bytes(outb),
                seconds=round(time.time() - t0, 1))   # RIGOR: info: timing
        if r.returncode != 0: failed('gate ' + g[0], f'exit {r.returncode}')
        print(f'gate {" ".join(g)}: exit {r.returncode}', flush=True)
    # ---- merges
    expected = {}                                                     # path -> sha256 the replays must consume
    replays = [(log, fam, man, 'replay2.py', tag) for log, fam, man, tag in LEGACY if only is None or tag in only]
    for log, fam, man, script, tag in replays: expected[log] = pins.get(log); expected[man] = pins.get(man)
    for F in FAMILIES:
        if only is not None and F.lower() not in only: continue
        L = F.lower(); merged = os.path.join(a.certs, f'{L}.merged.jsonl')
        inputs = sorted(f for f in os.listdir(os.path.join(HERE, a.certs)) if f.startswith(L + '.') and f.endswith('.jsonl') and not f.endswith('.merged.jsonl')) if os.path.isdir(os.path.join(HERE, a.certs)) else []
        in_sha = {f: sha_file(os.path.join(HERE, a.certs, f)) for f in inputs}
        if in_sha != {f: h for f, h in cc['inputs'].items() if f.startswith(L + '.')}: failed('merge ' + F, 'merge inputs differ from the files the completeness check read')
        r = subprocess.run(child('merge_pool.py', F, '--dir', a.certs), cwd=HERE, env=child_env, capture_output=True)
        mp = os.path.join(HERE, merged)
        msha = sha_file(mp) if (r.returncode == 0 and os.path.isfile(mp)) else None
        rc_.add('merge', family=F, exit=r.returncode, inputs=in_sha, merged=merged, merged_sha256=msha, output_sha256=sha_bytes(r.stdout + r.stderr))
        if msha is None: failed('merge ' + F, f'exit {r.returncode} or merged log missing'); continue
        man = f'manifests/{L}.json'
        expected[merged] = msha; expected[man] = pins.get(man)
        replays.append((merged, F, man, 'replay3.py', L))
    # ---- replays
    out_sha = {}
    for log, fam, man, script, tag in replays:
        stage = f'replay {tag}'; t0 = time.time()
        for b in check_cache_inventory(pins): failed(stage, 'chart caches: ' + b)
        ls = rcm.snapshot(os.path.join(HERE, log)); ms = rcm.snapshot(os.path.join(HERE, man))
        if expected.get(log) is None or ls['sha256'] != expected[log]: failed(stage, f'log {log}: sha256 {ls["sha256"][:12]} differs from its pin / merge entry'); continue
        if expected.get(man) is None or ms['sha256'] != expected[man]: failed(stage, f'manifest {man}: sha256 differs from its pin entry'); continue
        procs = []
        for i in range(N):
            f = os.path.join(run_dir, f'{tag}.{i}.jsonl')
            lf = open(os.path.join(run_dir, f'{tag}.{i}.log'), 'xb')
            procs.append((i, f, subprocess.Popen(child(script, log, fam, man, '--shard', f'{i}/{N}', '--out', f) + (['--no-recert'] if quick else []),
                                                 cwd=HERE, env=child_env, stdout=lf, stderr=subprocess.STDOUT), lf))
        codes = {}; outputs = []
        for i, f, p, lf in procs:
            codes[i] = p.wait(); lf.close()
            if codes[i] != 0: failed(stage, f'shard {i} exit {codes[i]}')
            try:
                with open(f, 'rb') as fh: b = fh.read()                    # the ONLY read of this output
                outputs.append((os.path.basename(f), b)); out_sha[f] = sha_bytes(b)
            except OSError: failed(stage, f'shard {i}: output missing')
        summaries, sbad = summaries_of(outputs)
        for b_ in sbad: failed(stage, b_)
        vh = summaries[0].get('binding_sha256') if summaries and isinstance(summaries[0], dict) else None
        problems = validate(outputs, ls, ms, fam, vh, expect_shards=N, recert=not quick) if vh else ['no binding in the first summary']
        # the binding a child computed (which this process cannot recompute: other loaded images, environment, pycache prefix) must
        # report exactly the pinned sources and chart caches, the sanitised child environment and the run's pycache prefix
        for k, sm in enumerate(summaries):
            if not isinstance(sm, dict): continue
            for src, h in (sm.get('binding_sources') or {}).items():
                if pins.get(src) != h: problems.append(f'shard {k}: source {src} differs from the pin')
            if not sm.get('binding_sources'): problems.append(f'shard {k}: no binding sources reported')
            if {('refcache/' + f): h for f, h in (sm.get('binding_refcache') or {}).items()} != {f: h for f, h in pins.items() if f.startswith('refcache/')}:
                problems.append(f'shard {k}: chart caches differ from the pin')
            env_ = sm.get('environment') or {}
            if env_.get('env') != {k2: v for k2, v in child_env.items() if k2.startswith(('OMP_', 'OPENBLAS_', 'VECLIB_', 'DARTS_'))}: problems.append(f'shard {k}: environment differs from the sanitised child environment')
            if env_.get('pycache_prefix') != os.path.join(run_dir, 'pycache'): problems.append(f'shard {k}: wrong pycache prefix')
            if (env_.get('fp') or {}).get('behaviour_failures'): problems.append(f'shard {k}: FP environment')
            if sm.get('family') != fam or sm.get('log_sha256') != expected[log] or sm.get('manifest_sha256') != expected[man]: problems.append(f'shard {k}: family / log / manifest identity')
        if len({(sm or {}).get('environment', {}).get('inventory_sha256') for sm in summaries}) != 1: problems.append('shards report different runtime inventories')
        for pr in problems[:20]: failed(stage, pr)
        if len(problems) > 20: failed(stage, f'... {len(problems) - 20} more validation problems')
        rc_.add('replay', tag=tag, script=script, family=fam, log=log, log_sha256=ls['sha256'], manifest=man, manifest_sha256=ms['sha256'],
                expected_binding=vh, shards=N, exit_codes=codes, outputs={n: sha_bytes(b) for n, b in outputs}, summaries=summaries,
                validation_problems=len(problems), seconds=round(time.time() - t0, 1))   # RIGOR: info: timing
        print(f'{stage}: shards ok {sum(1 for c in codes.values() if c == 0)}/{N}, validation problems {len(problems)}', flush=True)
        del ls, ms, outputs
    # ---- end of run
    end_bad = []
    for p in (pin_path, pin_copy):
        try:
            if sha_file(p) != pin_sha: end_bad.append(f'pin file {p} changed during the run')
        except OSError: end_bad.append(f'pin file {p} unreadable at the end')
    end_bad += check_pinned(pins) + check_cache_inventory(pins)
    for path, h in expected.items():
        if h is not None and os.path.isfile(os.path.join(HERE, path)) and sha_file(os.path.join(HERE, path)) != h: end_bad.append(f'{path} changed during the run')
    for f, h in out_sha.items():
        if sha_file(f) != h: end_bad.append(f'shard output {os.path.basename(f)} changed after it was validated')
    inv1 = static_inventory(paths=list(inv0))                       # the start inventory's files, re-read now (no memoisation)
    if inv1 != inv0: end_bad.append('runtime files changed during the run: ' + str(sorted(k for k in inv0 if inv0.get(k) != inv1.get(k))[:5]))
    srcs = {f: pins[f] for f in pins if f.endswith('.py') and '/' not in f}
    late = {}
    end_bad += loaded_modules_problems(inv0, srcs, late)
    for d in (os.path.join(run_dir, 'pycache'),):
        if any(fn for _, _, fn in os.walk(d)): end_bad.append(f'bytecode written under {d}')
    try: fpenv.check()
    except fpenv.FPEnvError as e: end_bad.append(str(e))
    for b in end_bad: failed('end', b)
    status = final_status(a.mode, failures)
    digest = rc_.add('end', end_state_unchanged=not end_bad, end_problems=end_bad, failures=failures, status=status, driver_late_loaded=late)
    rc_.close()
    with open(rc_.path, 'rb') as fh: rraw = fh.read()
    recs, final, rprob = read_receipt(rraw)
    if rprob or final != digest or recs[-1]['data'].get('status') != status:
        print('FAIL receipt re-read:', rprob[:3]); status = 'FAIL'
    print(f'REPLAY_ALL STATUS: {status}  receipt {os.path.relpath(rc_.path, HERE)}  run directory {os.path.relpath(run_dir, HERE)}')
    print(f'RECEIPT FINAL DIGEST (self_sha256 of the end record): {final}')
    print(f'RECEIPT FILE SHA256: {sha_bytes(rraw)}')
    sys.exit(exit_code(status))

if __name__ == '__main__':
    main()
