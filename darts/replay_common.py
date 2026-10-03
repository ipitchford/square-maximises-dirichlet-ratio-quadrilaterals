"""Shared fail-closed plumbing of replay2.py / replay3.py (adversarial review round 2, SERIOUS 3, 4, 6).
  snapshot(path)        reads a file ONCE; the replay parses exactly these bytes and records their sha256 (no re-read for the hash)
  changed(snap)         re-reads the path at the end: a different sha256 or file identity (device, inode, size, mtime) is a failure
  start_checks()        at the start of a replay: interpreter flags of a proof process (isolated, -S, -B, empty pycache prefix, no -O),
                        FP environment, and in a proof run (DARTS_PIN) the chart-cache inventory equal to the pin BEFORE any chart
                        is loaded
  end_checks(...)       at the end of a replay: binding recomputed from a FRESH file-system inventory and equal to the start binding
                        (sources, chart cache, manifest hash, family, runtime identity, FP environment); every snapshot unchanged;
                        FP environment still round-to-nearest with gradual underflow; every loaded module and image comes from an
                        unchanged file of the start inventory or a pinned local source and read no cached bytecode; the chart-cache
                        inventory still equal to the pin; no production-only module was loaded (this verifies the RIGOR-MODULE
                        annotations of rigor_lint.py at run time)."""
import os, sys, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
# production-only modules: rigor_lint.py exempts them from the per-site rules ONLY because no replay ever executes them
EXEMPT_MODULES = ('cover2', 'cover3', 'cover4', 'cover5', 'bo2f', 'manifests_legacy', 'pool', 'recert', 'run_spike_all')

def snapshot(path):
    with open(path, 'rb') as fh:
        st = os.fstat(fh.fileno()); raw = fh.read()
    return dict(path=path, bytes=raw, sha256=hashlib.sha256(raw).hexdigest(),
                ident=(st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns))

def changed(snap):
    try:
        with open(snap['path'], 'rb') as fh:
            st = os.fstat(fh.fileno()); raw = fh.read()
    except OSError as e:
        return f'{snap["path"]}: unreadable at the end ({e})'
    if hashlib.sha256(raw).hexdigest() != snap['sha256']: return f'{snap["path"]}: content changed during the replay'
    if (st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns) != snap['ident']: return f'{snap["path"]}: file identity changed during the replay'
    return None

def start_checks():
    """problems that forbid starting a replay"""
    import fpenv
    from binding import proof_flags_problems
    bad = proof_flags_problems()
    try: fpenv.check()
    except fpenv.FPEnvError as e: bad.append(str(e))
    if os.environ.get('DARTS_PIN'):
        import chart
        bad += chart.verify_cache_inventory()
    return bad

LATE = {}           # modules / images first loaded during the run, verified against their installer RECORD (reported in summaries)
def end_checks(script, manifest_sha, family, vh_start, snaps, start_det):
    import fpenv
    from binding import binding, loaded_modules_problems, static_inventory
    bad = []
    try: fpenv.check()
    except fpenv.FPEnvError as e: bad.append(str(e))
    inv_end = static_inventory(paths=list(start_det['inventory']))     # the same files, re-read now (no memoisation)
    vh_end, det_end = binding(script, manifest_sha, family, inventory=inv_end)
    if vh_end != vh_start:
        bad.append('binding changed during the replay (sources, chart cache, runtime or environment)')
        a, b = start_det['inventory'], det_end['inventory']
        ch = sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))
        if ch: bad.append(f'runtime files changed: {ch[:5]}')
    for s in snaps:
        w = changed(s)
        if w: bad.append(w)
    srcs = dict(start_det['sources'])
    srcs['bootstrap.py'] = __import__('binding').sha_file(os.path.join(HERE, 'bootstrap.py'))
    bad += loaded_modules_problems(start_det['inventory'], srcs, LATE)
    if os.environ.get('DARTS_PIN'):
        import chart
        bad += chart.verify_cache_inventory()
    for m in EXEMPT_MODULES:
        if m in sys.modules: bad.append(f'production-only module {m} was loaded by the replay')
    return bad
