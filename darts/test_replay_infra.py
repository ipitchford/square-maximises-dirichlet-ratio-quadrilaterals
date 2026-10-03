"""Probes of the replay / receipt machinery (adversarial reviews round 2 (SERIOUS 2-7) and round 3 (items 1-7)).  Threat model:
accidental corruption, stale or changed files, concurrent agents editing files, environment drift.  Every probe must FAIL CLOSED
(or, for the positive probes, pass).  Scratch space: replays/test_runs/infra_<random>/ (inside the package).
  R2   SHARDS / QUICK validation; only a proof-mode run without failures is PASS and exits 0 (quick / test exit 3, FAIL exits 1)
  R2   wrong family; replay_ok: genuine replay accepted, wrong family / summary-only / dropped / forged receipts rejected; a log
       changed during a replay; FE_UPWARD / FE_DOWNWARD / FE_TOWARDZERO / FZ at start
  1    runtime inventory recomputed at every call (no memoisation); the loaded Python framework binary is bound; a stale .pyc with
       a valid header is executed by a normal interpreter but NOT by the proof flags (-B, empty pycache prefix); sys.flags /
       _xoptions / pycache prefix are part of the identity; python -O is refused; a late-loaded module outside the inventory and
       its installer RECORD, or differing from its RECORD, is reported
  2    a pin without chart caches aborts the driver; a cache file whose bytes differ from the pin is NOT unpickled; an extra cache
       file is reported
  3    the pin file changed during a run, a pinned log differing from its pin entry, a merged log changed after its merge: FAIL
  4    outputs are created exclusively; a shard output changed after validation is caught at the end; a missing summary FAILS
  5    receipts: a forged terminal status, reordered / extra-key / non-canonical / torn records are rejected; write loop under
       partial writes; the printed final digest is the end record's self hash
  6    a symlinked copy of the launcher next to a counterfeit driver still runs the canonical driver; a PATH shim 'nice' is never
       executed; a polluted caller environment (PYTHONPATH, PYTHONSTARTUP, other variables) does not reach the children
  7    lint: NumPy aliases, from-numpy names, ufunc reductions, np.all, .var() flagged; an annotation inside a string literal is
       ignored; Lemma 9.6 interval check rejects a violating configuration; FIZ / AH / FZ16 (where the CPU supports them) and a
       rounding-mode change inside a leaf are caught; receipt metadata differing from the log record is rejected
  4b   launch path under research-run (which execs its job as /usr/bin/nice -n 10 CMD): the replay driver and recert_all.sh run
       at the TARGET niceness max(10, launcher's) -- under nice -n 10, 10 and not 20 -- and so do all their children; default
       worker counts are 6; recert_all.sh works from any cwd under a PATH of failing shims and a polluted environment; it records
       N, the code identity and the families at the first start and refuses a resume with another N (including the default),
       other code; a worker refuses invalid claims and a second process with the same worker ID; a torn last line is cut off at
       the next start; certs_complete.check flags every kind of incompleteness (missing / duplicate / not-ok root_done, FAIL
       leaf, handoff without item_done, unclaimed or orphaned root, garbage line, stray file, worker ID above N, no RECERT_N); the
       driver aborts on an incomplete directory before any gate, merge or replay, and replays a complete re-certified family
Usage: test_replay_infra.py   (exit 0 iff every probe behaves as required; about 10 minutes)"""
import os, sys, json, shutil, subprocess, hashlib, threading, time, secrets, ctypes, importlib
HERE = os.path.dirname(os.path.realpath(__file__)); sys.path.insert(0, HERE)
PY = os.environ.get('PY') or sys.executable
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
RESULTS = []
def probe(name, cond, detail=''):
    RESULTS.append((name, bool(cond))); print(('ok   ' if cond else 'FAIL ') + name + (f'  ({detail})' if detail else ''), flush=True)
def run(cmd, env=None, **kw):
    return subprocess.run(cmd, cwd=HERE, env=env or ENV, capture_output=True, text=True, **kw)
SCR = os.path.join(HERE, 'replays', 'test_runs', 'infra_' + secrets.token_hex(4)); os.makedirs(SCR)
RSCR = os.path.relpath(SCR, HERE)
def make_pin(path, drop=None, alter=None):
    r = run([PY, 'make_pin.py', path]); assert r.returncode == 0, r.stdout + r.stderr
    lines = open(path).read().split('\n')[:-1]
    if drop: lines = [l for l in lines if not drop(l.split('  ', 1)[1])]
    if alter:
        lines = [(('0'*64 + '  ' + l.split('  ', 1)[1]) if alter(l.split('  ', 1)[1]) else l) for l in lines]
    open(path, 'w').write('\n'.join(lines) + '\n')
    return path
def driver(run_id, pin, only='crit_collar', certs=None, env_extra=None, wrapper=os.path.join(HERE, 'replay_all.sh'), path=None, prefix=()):
    env = dict(ENV, SHARDS='2', PIN=pin); env.update(env_extra or {})
    if path: env['PATH'] = path
    args = list(prefix) + ['bash', wrapper, '--mode', 'test', '--only', only, '--gates', 'none', '--run-id', run_id]
    if certs: args += ['--certs', certs]
    return subprocess.run(args, cwd=HERE, env=env, capture_output=True, text=True)
def receipt_path(run_id): return os.path.join(HERE, 'replays', 'test_runs', 'receipts', run_id + '.jsonl')
def wait_for(pred, timeout=300):
    t0 = time.time()
    while time.time() - t0 < timeout:
        try:
            if pred(): return True
        except Exception: pass
        time.sleep(0.2)
    return False
def receipt_kinds(run_id):
    try: return [json.loads(l)['kind'] + ':' + str(json.loads(l)['data'].get('tag', json.loads(l)['data'].get('family', ''))) for l in open(receipt_path(run_id))]
    except OSError: return []

def main():
    import replay_all as RA
    BASE_PIN = make_pin(os.path.join(SCR, 'pin.sha256'))
    # ---------------------------------------------------------------- round-2 probes (kept)
    for s in ('0', '-1', 'abc', '07', '65'):
        r = run(['bash', 'replay_all.sh'], env=dict(ENV, SHARDS=s)); probe(f'R2 SHARDS={s!r} rejected', r.returncode == 2 and 'SHARDS must be' in r.stdout)
    for q in ('0', 'false', 'yes'):
        r = run(['bash', 'replay_all.sh'], env=dict(ENV, QUICK=q)); probe(f'R2 QUICK={q!r} rejected', r.returncode == 2 and 'QUICK must be' in r.stdout)
    probe('R2/7 exit codes: 0 only for a proof PASS', RA.exit_code(RA.final_status('proof', [])) == 0 and RA.exit_code(RA.final_status('quick', [])) == 3
          and RA.exit_code(RA.final_status('test', [])) == 3 and RA.exit_code(RA.final_status('proof', ['x'])) == 1 and 'PASS' not in RA.final_status('quick', []))
    out = os.path.join(SCR, 'ok'); os.makedirs(out)
    r = run([PY, 'replay2.py', 'certs/crit_collar.jsonl', 'COLLAR', 'manifests/crit_collar.json', '--out', os.path.join(out, 'crit_collar.0.jsonl')])
    probe('genuine replay of crit_collar passes', r.returncode == 0, r.stdout[-200:] if r.returncode else '')
    r = run([PY, 'replay2.py', 'certs/crit_collar.jsonl', 'SPIKE_X1', 'manifests/crit_collar.json', '--out', os.path.join(SCR, 'wf.jsonl')])
    probe('R2 wrong-family replay fails', r.returncode != 0 and 'family' in r.stdout)
    ok_cmd = [PY, 'replay_ok.py', 'crit_collar', 'certs/crit_collar.jsonl', 'replay2.py', 'manifests/crit_collar.json']
    probe('R2 replay_ok accepts the genuine replay', run(ok_cmd + ['COLLAR', '--dir', out]).returncode == 0)
    probe('R2 replay_ok rejects a wrong expected family', run(ok_cmd + ['SPIKE_X1', '--dir', out]).returncode != 0)
    lines = open(os.path.join(out, 'crit_collar.0.jsonl')).read().split('\n')[:-1]
    def forged(name, new_lines):
        d = os.path.join(SCR, 'syn_' + name); os.makedirs(d)
        open(os.path.join(d, 'crit_collar.0.jsonl'), 'w').write('\n'.join(new_lines) + '\n')
        return run(ok_cmd + ['COLLAR', '--dir', d]).returncode
    s0 = json.loads(lines[-1]); s0['summary']['leaves_checked'] = s0['summary']['ok'] = 123; s0['summary'].pop('leaves_sha256')
    probe('R2 summary-only output rejected', forged('summary_only', [json.dumps(s0)]) != 0)
    probe('R2 dropped receipt rejected', forged('dropped', lines[:2] + lines[3:]) != 0)
    rec = json.loads(lines[0]); rec['record_sha256'] = '0'*64
    probe('R2 forged record hash rejected', forged('forged', [json.dumps(rec)] + lines[1:]) != 0)
    rec = json.loads(lines[0]); rec['status'] = 'discard' if rec['status'] == 'cert' else 'cert'
    probe('7 receipt metadata differing from its log record rejected', forged('meta', [json.dumps(rec)] + lines[1:]) != 0)
    logc = os.path.join(SCR, 'crit_collar_copy.jsonl'); shutil.copyfile(os.path.join(HERE, 'certs/crit_collar.jsonl'), logc)
    w = os.path.join(SCR, 'mutate.py')
    open(w, 'w').write(f'''import sys; sys.path.insert(0, {HERE!r})
import replay2
orig = replay2.tiles
def tiles(*a, **k):
    with open({logc!r}, 'a') as fh: fh.write('\\n')
    return orig(*a, **k)
replay2.tiles = tiles
sys.argv = ['replay2.py', {logc!r}, 'COLLAR', {os.path.join(HERE, 'manifests/crit_collar.json')!r}]
replay2.main()
''')
    r = run([PY, w]); probe('R2 log changed during the replay -> failure', r.returncode != 0 and 'changed during the replay' in r.stdout)
    import fpenv
    for mode, name in ((0x00400000, 'FE_UPWARD'), (0x00800000, 'FE_DOWNWARD'), (0x00C00000, 'FE_TOWARDZERO')):
        w = os.path.join(SCR, f'fp_{name}.py')
        open(w, 'w').write(f'''import sys; sys.path.insert(0, {HERE!r})
import fpenv; fpenv._lib().fesetround({mode})
import replay2
sys.argv = ['replay2.py', {os.path.join(HERE, 'certs/crit_collar.jsonl')!r}, 'COLLAR', {os.path.join(HERE, 'manifests/crit_collar.json')!r}]
replay2.main()
''')
        r = run([PY, w]); probe(f'R2 replay refuses to start under {name}', r.returncode != 0 and 'START CHECKS FAILED' in r.stdout)
    # FPCR bits: FZ, FZ16, AH, FIZ (AH and FIZ exist only with FEAT_AFP; an unsupported bit does not stick and cannot be tested)
    for bit, name in ((24, 'FZ'), (19, 'FZ16'), (1, 'AH'), (0, 'FIZ')):
        w = os.path.join(SCR, f'fpcr_{name}.py')
        open(w, 'w').write(f'''import sys, ctypes; sys.path.insert(0, {HERE!r})
import fpenv
env = fpenv._FenvArm64Darwin(); lib = fpenv._lib(); lib.fegetenv(ctypes.byref(env)); env.fpcr |= (1 << {bit}); lib.fesetenv(ctypes.byref(env))
lib.fegetenv(ctypes.byref(env))
if not (env.fpcr >> {bit}) & 1: print('UNSUPPORTED'); sys.exit(0)
try: fpenv.check(); print('NOT DETECTED')
except fpenv.FPEnvError as e: print('DETECTED', e)
''')
        r = run([PY, w])
        probe(f'7 FPCR.{name} detected (or unsupported by this CPU)', 'DETECTED' in r.stdout or 'UNSUPPORTED' in r.stdout, r.stdout.strip()[:60])
    w = os.path.join(SCR, 'fp_midleaf.py')
    open(w, 'w').write(f'''import sys; sys.path.insert(0, {HERE!r})
import fpenv, replay2, regions
orig = regions.check_reason
calls = [0]
def cr(*a, **k):
    calls[0] += 1
    if calls[0] == 1: fpenv._lib().fesetround(0x00400000)       # a persistent mode change inside the first leaf
    return orig(*a, **k)
regions.check_reason = cr
sys.argv = ['replay2.py', {os.path.join(HERE, 'certs/main_cbox.jsonl')!r}, 'MAIN_CBOX', {os.path.join(HERE, 'manifests/main_cbox.json')!r}, '--no-recert']
replay2.main()
''')
    r = run([PY, w]); probe('7 a rounding-mode change inside a leaf is caught by the per-leaf guard', r.returncode != 0 and 'round-to-nearest' in r.stdout)
    # ---------------------------------------------------------------- item 1: runtime inventory
    import binding
    t = os.path.join(SCR, 'invfile.bin'); open(t, 'wb').write(b'a')
    i1 = binding.static_inventory(paths=[t]); open(t, 'wb').write(b'b'); i2 = binding.static_inventory(paths=[t])
    probe('1 runtime inventory recomputed at every call (no memoisation)', i1[t] != i2[t])
    paths = binding.inventory_paths()
    probe('1 the loaded Python framework binary is bound', any(p.endswith('/Python') and 'Python.framework' in p for p in paths))
    ident = binding.runtime_identity()
    probe('1 sys.flags, _xoptions and the pycache prefix are part of the runtime identity', all(k in ident for k in ('flags', 'xoptions', 'pycache_prefix')))
    pm = os.path.join(SCR, 'pycmod'); os.makedirs(os.path.join(pm, '__pycache__'))
    src = os.path.join(pm, 'probe_mod.py'); open(src, 'w').write('X = 1\n')
    import importlib.util, importlib._bootstrap_external as be, marshal
    st = os.stat(src)
    code = compile('X = 2\n', src, 'exec')
    data = be._code_to_timestamp_pyc(code, st.st_mtime, st.st_size)
    open(importlib.util.cache_from_source(src), 'wb').write(data)
    r1 = run([PY, '-c', f'import sys; sys.path.insert(0, {pm!r}); import probe_mod; print(probe_mod.X)'])
    empty = os.path.join(SCR, 'pyc_empty'); os.makedirs(empty)
    r2 = run([PY, '-I', '-S', '-B', '-X', f'pycache_prefix={empty}', '-c', f'import sys; sys.path.insert(0, {pm!r}); import probe_mod; print(probe_mod.X)'])
    probe('1 a stale .pyc with a valid header: executed normally, ignored under the proof flags', r1.stdout.strip() == '2' and r2.stdout.strip() == '1', f'normal {r1.stdout.strip()}, proof {r2.stdout.strip()}')
    r = run([PY, '-O', '-I', '-S', '-B', '-X', f'pycache_prefix={empty}', os.path.join(HERE, 'bootstrap.py'), os.path.join(HERE, 'replay2.py'),
             'certs/crit_collar.jsonl', 'COLLAR', 'manifests/crit_collar.json'], env=dict(ENV, DARTS_PIN=BASE_PIN))
    probe('1 python -O is refused', r.returncode != 0 and 'refused' in r.stdout)
    fake = os.path.join(SCR, 'fake_late.py'); open(fake, 'w').write('Y = 1\n')
    import types
    m = types.ModuleType('fake_late_mod'); m.__spec__ = importlib.util.spec_from_file_location('fake_late_mod', fake); sys.modules['fake_late_mod'] = m
    pr = binding.loaded_modules_problems(binding.static_inventory(), {})
    probe('1 a late-loaded module outside the inventory and every RECORD is reported', any('fake_late.py' in p and 'RECORD' in p for p in pr))
    orig_rh = binding.record_hashes
    binding.record_hashes = lambda: {os.path.realpath(fake): ('fakedist', 'wronghash')}
    pr = binding.loaded_modules_problems(binding.static_inventory(), {})
    binding.record_hashes = orig_rh; del sys.modules['fake_late_mod']
    probe('1 a late-loaded module differing from its installer RECORD is reported', any('differs from its installer RECORD' in p for p in pr))
    # ---------------------------------------------------------------- item 2: chart caches
    nopin = make_pin(os.path.join(SCR, 'pin_nocache.sha256'), drop=lambda f: f.startswith('refcache/'))
    r = driver('p2a_' + secrets.token_hex(3), nopin)
    probe('2 a pin without the chart caches aborts the driver', r.returncode == 1 and 'the pin lists no chart cache files' in r.stdout)
    bad_cache = sorted(f for f in os.listdir(os.path.join(HERE, 'refcache')) if f.startswith('chart_'))[0]
    altpin = make_pin(os.path.join(SCR, 'pin_altcache.sha256'), alter=lambda f: f == 'refcache/' + bad_cache)
    w = os.path.join(SCR, 'cache_probe.py')
    open(w, 'w').write(f'''import sys, pickle; sys.path.insert(0, {HERE!r})
calls = []
orig = pickle.loads
pickle.loads = lambda *a, **k: (calls.append(1), orig(*a, **k))[1]
import chart
print('INVENTORY', chart.verify_cache_inventory())
try: chart.load_pinned(chart.os.path.join(chart.CACHE, {bad_cache!r})); print('LOADED')
except ValueError as e: print('REFUSED', e)
print('UNPICKLE_CALLS', len(calls))
''')
    r = run([PY, w], env=dict(ENV, DARTS_PIN=altpin))
    probe('2 a cache file whose bytes differ from the pin is not unpickled', 'REFUSED' in r.stdout and 'UNPICKLE_CALLS 0' in r.stdout and 'differs from the pin' in r.stdout)
    droppin = make_pin(os.path.join(SCR, 'pin_dropcache.sha256'), drop=lambda f: f == 'refcache/' + bad_cache)
    r = run([PY, w], env=dict(ENV, DARTS_PIN=droppin))
    probe('2 an extra (unpinned) cache file is reported and not loaded', 'extra' in r.stdout and bad_cache in r.stdout and 'UNPICKLE_CALLS 0' in r.stdout)
    # ---------------------------------------------------------------- item 3: pin and inputs
    pin3 = os.path.join(SCR, 'pin_changing.sha256'); shutil.copyfile(BASE_PIN, pin3)
    rid = 'p3a_' + secrets.token_hex(3)
    th = threading.Thread(target=lambda: wait_for(lambda: len(receipt_kinds(rid)) >= 1) and open(pin3, 'a').write('\n'))
    th.start(); r = driver(rid, pin3, only='crit_collar,ff', certs='recert_tests/pipe1'); th.join()
    probe('3 the pin file changed during the run -> FAIL', r.returncode == 1 and 'changed during the run' in r.stdout)
    altlog = make_pin(os.path.join(SCR, 'pin_altlog.sha256'), alter=lambda f: f == 'certs/crit_collar.jsonl')
    r = driver('p3b_' + secrets.token_hex(3), altlog)
    probe('3 a pinned log differing from its pin entry aborts', r.returncode == 1 and 'certs/crit_collar.jsonl: sha256 differs from the pin' in r.stdout)
    cdir = os.path.join(SCR, 'certs_ff'); os.makedirs(cdir)
    for f in os.listdir(os.path.join(HERE, 'recert_tests', 'pipe1')):
        if f.startswith('ff.w') and f.endswith('.jsonl'): shutil.copyfile(os.path.join(HERE, 'recert_tests', 'pipe1', f), os.path.join(cdir, f))
    shutil.copytree(os.path.join(HERE, 'recert_tests', 'pipe1', 'claims'), os.path.join(cdir, 'claims'))     # (complete: passes certs_complete)
    rid = 'p3c_' + secrets.token_hex(3)
    th = threading.Thread(target=lambda: wait_for(lambda: any(k.startswith('merge') for k in receipt_kinds(rid))) and open(os.path.join(cdir, 'ff.merged.jsonl'), 'a').write('\n'))
    th.start(); r = driver(rid, BASE_PIN, only='ff', certs=os.path.relpath(cdir, HERE)); th.join()
    probe('3 a merged log changed after its merge -> FAIL', r.returncode == 1 and ('merge entry' in r.stdout or 'changed during the run' in r.stdout or 'log / manifest identity' in r.stdout))
    # ---------------------------------------------------------------- item 4: shard outputs
    ex = os.path.join(SCR, 'exists.jsonl'); open(ex, 'w').write('x\n')
    r = run([PY, 'replay2.py', 'certs/crit_collar.jsonl', 'COLLAR', 'manifests/crit_collar.json', '--out', ex])
    probe('4 an existing output file is never overwritten (exclusive creation)', r.returncode != 0 and 'FileExistsError' in r.stderr and open(ex).read() == 'x\n')
    rid = 'p4b_' + secrets.token_hex(3)
    run_dir = os.path.join(HERE, 'replays', 'test_runs', f'run_{rid}')
    def tamper():
        if wait_for(lambda: 'replay:crit_collar' in receipt_kinds(rid)):
            with open(os.path.join(run_dir, 'crit_collar.0.jsonl'), 'a') as fh: fh.write('\n')
    th = threading.Thread(target=tamper); th.start(); r = driver(rid, BASE_PIN, only='crit_collar,ff', certs='recert_tests/pipe1'); th.join()
    probe('4 a shard output changed after its validation is caught at the end', r.returncode == 1 and 'changed after it was validated' in r.stdout)
    summ, bad = RA.summaries_of([('x.jsonl', lines[0].encode() + b'\n')])
    probe('4 a missing summary is a failure', summ == [None] and bad)
    # ---------------------------------------------------------------- items 5, 6: receipts and bootstrap (one clean machinery run)
    rid = 'p5_' + secrets.token_hex(3)
    r = driver(rid, BASE_PIN, env_extra={'PYTHONPATH': SCR, 'PYTHONSTARTUP': os.path.join(SCR, 'nope.py'), 'EVIL_VAR': '1'})
    raw = open(receipt_path(rid), 'rb').read()
    recs, digest, prob = RA.read_receipt(raw)
    probe('6 a polluted caller environment does not reach the children (machinery run passes)', r.returncode == 3 and not prob and recs[-1]['data']['status'].startswith('MACHINERY'), r.stdout.strip().split('\n')[-3][:100] if r.stdout else r.stderr[-200:])
    probe('5 the printed final digest is the end record\'s self hash', f'RECEIPT FINAL DIGEST (self_sha256 of the end record): {digest}' in r.stdout)
    start = recs[0]['data']
    probe('6 children run with the fixed sanitised environment', set(start['child_env']) == {'PATH', 'HOME', 'TMPDIR', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'DARTS_PIN'}
          and start['child_env']['PATH'] == '/usr/bin:/bin' and start['child_command'][1:5] == ['-I', '-S', '-B', '-X'])
    lines_r = raw.split(b'\n')[:-1]
    forged_end = lines_r[-1].replace(b'MACHINERY TEST OK (NOT A PROOF)', b'PASS')
    probe('5 a forged terminal status is rejected', bool(RA.read_receipt(b'\n'.join(lines_r[:-1] + [forged_end]) + b'\n')[2]))
    probe('5 reordered records are rejected', bool(RA.read_receipt(b'\n'.join([lines_r[1], lines_r[0]] + lines_r[2:]) + b'\n')[2]))
    probe('5 a non-canonical record is rejected', bool(RA.read_receipt(b'\n'.join([lines_r[0].replace(b'":', b'": ', 1)] + lines_r[1:]) + b'\n')[2]))
    probe('5 a torn receipt is rejected', bool(RA.read_receipt(raw[:-1])[2]))
    x = json.loads(lines_r[0]); x['extra'] = 1
    probe('5 a record with an extra key is rejected', bool(RA.read_receipt(json.dumps(x, sort_keys=True, separators=(',', ':')).encode() + b'\n')[2]))
    orig_write = RA.os.write
    RA.os.write = lambda fd, data: orig_write(fd, bytes(data[:3]))
    try:
        rp = os.path.join(SCR, 'partial', 'r.jsonl'); rr = RA.Receipt(rp); rr.add('start', a=1, b='x'*100); rr.add('end', status='FAIL'); rr.close()
    finally: RA.os.write = orig_write
    probe('5 the receipt write loop survives partial writes', not RA.read_receipt(open(rp, 'rb').read())[2])
    fk = os.path.join(SCR, 'fake'); os.makedirs(fk)
    os.symlink(os.path.join(HERE, 'replay_all.sh'), os.path.join(fk, 'replay_all.sh'))
    open(os.path.join(fk, 'replay_all.py'), 'w').write('print("COUNTERFEIT")\n')
    r = subprocess.run(['bash', os.path.join(fk, 'replay_all.sh')], cwd=fk, env=dict(ENV, SHARDS='2', PIN=os.path.join(SCR, 'no_such_pin')), capture_output=True, text=True)
    probe('6 a symlinked launcher next to a counterfeit driver runs the canonical driver', 'COUNTERFEIT' not in r.stdout and 'unreadable or malformed' in r.stdout)
    sb = os.path.join(SCR, 'shimbin'); os.makedirs(sb); marker = os.path.join(SCR, 'nice_was_run')
    open(os.path.join(sb, 'nice'), 'w').write(f'#!/bin/sh\ntouch {marker}\nexec /usr/bin/nice "$@"\n'); os.chmod(os.path.join(sb, 'nice'), 0o755)
    r = driver('p6b_' + secrets.token_hex(3), BASE_PIN, path=sb + ':' + ENV.get('PATH', '/usr/bin:/bin'))
    probe('6 a PATH shim named nice is never executed', r.returncode == 3 and not os.path.exists(marker) and "'nice'" not in open(os.path.join(HERE, 'replay_all.py')).read())
    # ---------------------------------------------------------------- item 7: lint and Lemma 9.6
    import rigor_lint as RL, ast
    srcl = ('import numpy as xp\nfrom numpy import max as npmax, sum\nimport numpy\n'
            'a = npmax(x)\nb = sum(y)\nc = xp.sum(z)\nd = numpy.add.reduce(w)\ne = xp.all(v)\nf = q.var()\n'
            's = "# RIGOR: exact: inside a string"; g = float(s)\n')
    kinds = [k for _, k in RL.flagged(ast.parse(srcl))]
    probe('7 lint flags NumPy aliases, from-numpy names, ufunc reductions, np.all and .var()', all(any(t in k for k in kinds) for t in
          ('npmax (from numpy)', 'sum (from numpy)', 'np.sum', 'ufunc.reduce', 'np.all', '.var()')), str(kinds))
    probe('7 an annotation inside a string literal is not an annotation', 11 not in RL.annotations(srcl))
    from test_containment import farl_conditions_interval
    from fractions import Fraction as Fr
    good = farl_conditions_interval(Fr(5), Fr(2), Fr(1, 4), Fr(1), Fr(3), [(Fr(1, 2), Fr(-1))])
    badc = farl_conditions_interval(Fr(5), Fr(2), Fr(1, 4), Fr(1), Fr(3), [(Fr(1, 2), Fr(3))])
    probe('7 the Lemma 9.6 interval check accepts a valid and rejects a violating configuration', all(good[:3]) and not badc[1])
    launch_path_probes(RA, BASE_PIN)
    bad = [n for n, c in RESULTS if not c]
    print(f'test_replay_infra: {"PASS" if not bad else "FAIL"} ({len(RESULTS)} probes, {len(bad)} failed); scratch {RSCR}')
    sys.exit(0 if not bad else 1)

def launch_path_probes(RA, BASE_PIN):
    import certs_complete as CC, fcntl
    base = os.nice(0); want10 = max(10, base + 10); want0 = max(10, base)
    # ---- replay driver: target niceness
    rid = 'n10_' + secrets.token_hex(3)
    r = driver(rid, BASE_PIN, prefix=('/usr/bin/nice', '-n', '10'))
    recs = RA.read_receipt(open(receipt_path(rid), 'rb').read())[0] if os.path.exists(receipt_path(rid)) else []
    nd = recs[0]['data'].get('niceness') if recs else None
    nc = [s_.get('niceness') for rr in recs if rr['kind'] == 'replay' for s_ in rr['data']['summaries']]
    probe(f'4b driver launched under nice -n 10 (from {base}): driver and every child at niceness {want10}, not {base + 20}',
          r.returncode == 3 and nd == want10 and nc and all(x == want10 for x in nc), f'driver {nd}, children {nc}')
    rid = 'n0_' + secrets.token_hex(3)
    r = driver(rid, BASE_PIN)
    recs = RA.read_receipt(open(receipt_path(rid), 'rb').read())[0] if os.path.exists(receipt_path(rid)) else []
    nd = recs[0]['data'].get('niceness') if recs else None
    nc = [s_.get('niceness') for rr in recs if rr['kind'] == 'replay' for s_ in rr['data']['summaries']]
    probe(f'4b driver launched at niceness {base}: driver and every child at {want0}', r.returncode == 3 and nd == want0 and nc and all(x == want0 for x in nc), f'driver {nd}, children {nc}')
    sh = open(os.path.join(HERE, 'replay_all.sh')).read(); rsh = open(os.path.join(HERE, 'recert_all.sh')).read()
    probe('4b default worker counts are 6 (replay_all.sh SHARDS, replay_all.py --shards, recert_all.sh N)',
          'SH="${SHARDS:-6}"' in sh and "ap.add_argument('--shards', default='6')" in open(os.path.join(HERE, 'replay_all.py')).read() and 'N="${1:-6}"' in rsh)
    # ---- recert_all.sh on a test directory (family FF, 42 roots, 540 leaves) under nice -n 10, from cwd /, with a PATH of shims
    rdir = os.path.join(SCR, 'rc_ff'); shim = os.path.join(SCR, 'shim_rc'); os.makedirs(shim)
    for tool in ('nice', 'ps', 'mkdir', 'env', 'python', 'python3', 'cat', 'date', 'tr', 'mktemp', 'rmdir'):
        fn = os.path.join(shim, tool); open(fn, 'w').write(f'#!/bin/sh\n: > {os.path.join(SCR, "shim_used_" + tool)}\nexit 97\n'); os.chmod(fn, 0o755)
    env_rc = dict(PATH=shim + ':/usr/bin:/bin', HOME=os.environ.get('HOME', '/var/empty'), RECERT_DIR=rdir, RECERT_FAMILIES='FF',
                  PYTHONPATH=SCR, PYTHONSTARTUP=os.path.join(SCR, 'nope.py'), PYTHONHOME='/nonexistent', EVIL_VAR='1')
    def recert_all(*args, prefix=()):
        return subprocess.run(list(prefix) + ['/bin/bash', os.path.join(HERE, 'recert_all.sh')] + list(args), cwd='/', env=env_rc, capture_output=True, text=True, timeout=3600)
    def starts():
        out = []
        for w in (1, 2):
            try: out += [json.loads(l) for l in open(os.path.join(rdir, f'recert.w{w}.out')) if l.startswith('{"event": "start"')]
            except OSError: pass
        return out
    t0 = time.time(); r = recert_all('2', prefix=('/usr/bin/nice', '-n', '10')); dt = time.time() - t0
    st = starts()
    probe(f'4b recert_all.sh 2 under nice -n 10 (cwd /, shim PATH, polluted environment): complete, every worker at niceness {want10}',
          r.returncode == 0 and '"complete": true' in r.stdout and len(st) == 2 and all(x['niceness'] == want10 and x['isolated'] == 1 and x['no_site'] == 1 for x in st)
          and not any(f.startswith('shim_used_') for f in os.listdir(SCR)), f'exit {r.returncode}, {dt:.0f} s, niceness {[x["niceness"] for x in st]}, {r.stdout.strip()[-160:]}')
    r = recert_all('2')
    st = [ev for w_ in (1, 2) for ev in [e for e in starts() if e.get('worker') == w_][-1:]]   # the LAST start of each worker (the resume's)
    probe(f'4b resume with the same N from niceness {base}: workers at {want0}, nothing redone, complete', r.returncode == 0 and len(st) == 2 and all(x['niceness'] == want0 for x in st), f'exit {r.returncode}')
    r = recert_all('3')
    probe('4b resume with another N is refused (exit 2)', r.returncode == 2 and 'refusing to resume' in r.stdout and "RECERT_N records '2'" in r.stdout)
    r = recert_all()
    probe('4b resume with the default N (6) on a directory started with 2 is refused', r.returncode == 2 and "this start has '6'" in r.stdout)
    code_f = os.path.join(rdir, 'RECERT_CODE'); code = open(code_f).read()
    os.chmod(code_f, 0o644); open(code_f, 'w').write('0'*64 + '\n')
    r = recert_all('2'); open(code_f, 'w').write(code)
    probe('4b resume with other re-certification code is refused', r.returncode == 2 and 'RECERT_CODE' in r.stdout and 'refusing' in r.stdout)
    r = run([PY, 'recert.py', 'FF', '3'], env=dict(ENV, RECERT_DIR=rdir))
    probe('4b a worker ID above the recorded N is refused', r.returncode != 0 and 'is not in 1..2' in r.stderr)
    croot = sorted(os.listdir(os.path.join(rdir, 'claims', 'FF')))[0]; cfn = os.path.join(rdir, 'claims', 'FF', croot); owner = open(cfn).read()
    open(cfn, 'w').write('7')
    r = run([PY, 'recert.py', 'FF', '1'], env=dict(ENV, RECERT_DIR=rdir)); cc = CC.check(rdir, ['FF'])
    open(cfn, 'w').write(owner)
    probe('4b a claim naming no worker 1..N stops the worker and the completeness check', r.returncode != 0 and 'invalid claims' in r.stderr
          and any("owner '7'" in p for p in cc['problems']))
    lk = os.open(os.path.join(rdir, 'locks', 'w1.lock'), os.O_RDWR); fcntl.flock(lk, fcntl.LOCK_EX)
    r = run([PY, 'recert.py', 'FF', '1'], env=dict(ENV, RECERT_DIR=rdir)); os.close(lk)
    probe('4b a second process with the same worker ID is refused', r.returncode != 0 and 'already running' in r.stderr)
    w1 = os.path.join(rdir, 'ff.w1.jsonl'); size = os.path.getsize(w1)
    with open(w1, 'ab') as fh: fh.write(b'{"type": "leaf", "root": "FF_')
    torn = any('torn last line' in p for p in CC.check(rdir, ['FF'])['problems'])
    r = recert_all('2')
    probe('4b a torn last line is reported, then cut off at the next start (directory complete again)', torn and r.returncode == 0 and os.path.getsize(w1) == size
          and 'torn_tail_removed' in open(os.path.join(rdir, 'recert.w1.out')).read())
    # ---- completeness check: one defect per copy of the complete directory
    lines1 = open(w1).read().split('\n')[:-1]
    i_done = [i for i, l in enumerate(lines1) if '"root_done"' in l]; i_leaf = [i for i, l in enumerate(lines1) if '"type": "leaf"' in l]
    def variant(name, edit_lines=None, fs=None):
        d = os.path.join(SCR, 'cc_' + name); shutil.copytree(rdir, d, ignore=shutil.ignore_patterns('locks'))
        if edit_lines:
            ls = edit_lines(list(lines1)); open(os.path.join(d, 'ff.w1.jsonl'), 'w').write('\n'.join(ls) + '\n')
        if fs: fs(d)
        return d, CC.check(d, ['FF'])
    def setok(l, v): d = json.loads(l); d['ok'] = v; return json.dumps(d)
    def setfail(l): d = json.loads(l); d['status'] = 'FAIL'; d['res'] = 'probe'; d.pop('method', None); return json.dumps(d)
    root0 = json.loads(lines1[i_done[0]])['root']
    cases = [('missing_root_done', dict(edit_lines=lambda ls: ls[:i_done[0]] + ls[i_done[0] + 1:]), '0 root_done records'),
             ('not_ok', dict(edit_lines=lambda ls: ls[:i_done[0]] + [setok(ls[i_done[0]], False)] + ls[i_done[0] + 1:]), 'root_done has ok = False'),
             ('dup_root_done', dict(edit_lines=lambda ls: ls + [ls[i_done[0]]]), '2 root_done records'),
             ('handoff', dict(edit_lines=lambda ls: ls + [json.dumps(dict(type='handoff', root=root0, path=root0 + '0', p0=[0, 0, 0, 0], r=[0, 0, 0, 0], depth=1, k=4))]), 'has no item_done'),
             ('garbage', dict(edit_lines=lambda ls: ls[:5] + ['{"type": "leaf", "ro'] + ls[5:]), 'not JSON'),
             ('dup_tree', dict(edit_lines=lambda ls: ls + [ls[i_leaf[0]]]), 'duplicate tree record'),
             ('unclaimed', dict(fs=lambda d: os.remove(os.path.join(d, 'claims', 'FF', root0))), 'never claimed'),
             ('stray_file', dict(fs=lambda d: open(os.path.join(d, 'ff.0.jsonl'), 'w').write('')), 'not a worker log name'),
             ('id_above_n', dict(fs=lambda d: os.rename(os.path.join(d, 'ff.w2.jsonl'), os.path.join(d, 'ff.w5.jsonl'))), 'worker ID above N'),
             ('no_lock', dict(fs=lambda d: os.remove(os.path.join(d, 'RECERT_N'))), 'RECERT_N missing')]
    for name, kw, want in cases:
        d, cc = variant(name, **kw)
        probe(f'4b completeness check flags: {name}', any(want in p for p in cc['problems']), (cc['problems'][:1] or ['-'])[0][:90])
    dfail, cc = variant('fail_leaf', edit_lines=lambda ls: ls[:i_leaf[0]] + [setfail(ls[i_leaf[0]])] + ls[i_leaf[0] + 1:])
    probe('4b completeness check lists a FAIL leaf', len(cc['fails']) == 1 and cc['fails'][0]['path'] == json.loads(lines1[i_leaf[0]])['path'])
    probe('4b completeness check passes the complete directory', CC.check(rdir, ['FF']) ['problems'] == [] and not CC.check(rdir, ['FF'])['fails'])
    # ---- driver: an incomplete directory aborts before any gate / merge / replay; the complete one replays (FF recert -> merge -> replay3)
    rid = 'cc_' + secrets.token_hex(3)
    r = driver(rid, BASE_PIN, only='ff', certs=os.path.relpath(dfail, HERE))
    kinds = [k.split(':')[0] for k in receipt_kinds(rid)]
    probe('4b driver aborts on a FAIL leaf before any gate, merge or replay, listing it', r.returncode == 1 and kinds == ['start', 'complete', 'abort'] and 'FAIL LEAF:' in r.stdout, str(kinds))
    rid = 'ccok_' + secrets.token_hex(3)
    r = driver(rid, BASE_PIN, only='ff', certs=os.path.relpath(rdir, HERE))
    kinds = [k.split(':')[0] for k in receipt_kinds(rid)]
    probe('4b driver replays the complete re-certified family (recert -> completeness -> merge -> replay3, 2 shards)',
          r.returncode == 3 and kinds == ['start', 'complete', 'merge', 'replay', 'end'], str(kinds) + ' ' + r.stdout.strip()[-120:])

if __name__ == '__main__':
    main()
