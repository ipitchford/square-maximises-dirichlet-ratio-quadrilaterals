"""Binding of replay results (adversarial reviews of 2026-10-02, rounds 1-3).  Threat model (README §5): accidental corruption,
stale or changed files, concurrent agents editing files and environment drift; deliberate tampering by the machine owner is out of
scope.  A replay result is valid only for the exact inputs and runtime it was computed with.

binding(script, manifest_sha256, family) returns (sha256, details) over
  * every TRANSITIVE local source of the replay script (AST import closure), byte for byte;
  * every cached chart file (refcache/*), byte for byte (the proof driver additionally requires them to equal the reviewed pin
    before any of them is unpickled, see chart.py);
  * the sha256 of the manifest bytes that the replay PARSED and the FAMILY name;
  * the runtime identity runtime_identity(): a FRESH (never memoised) inventory static_inventory() of
      - every Mach-O image loaded in this process that exists as a file (dyld: the interpreter, the Python framework binary,
        every extension module and its bundled libraries),
      - the interpreter executable and the standard library tree (sources, extension modules),
      - every installed file of python-flint, numpy, scipy, scikit-fem, mpmath, and every *.pth file and every top-level module
        (*.py: sitecustomize, usercustomize, targets of *.pth import hooks such as virtualenv's _virtualenv.py) of the site
        directories,
    and sys.flags, sys._xoptions, sys.pycache_prefix, the NumPy build / BLAS configuration, the FLINT version and arb context,
    the OS, the floating-point environment (fpenv.state()) and the numerics environment variables.
loaded_modules_problems(start_inventory, sources): every module loaded in the process (at the time of the call) must come from a
file of the start inventory or a pinned local source with an unchanged sha256, and must not have read a cached .pyc (the proof
runs use an empty pycache prefix and -B: all bytecode is compiled from the bound sources).
proof_flags_problems(): in a proof run (DARTS_PIN set) the interpreter must run isolated (-I), without the site module (-S: no
*.pth startup hooks), with -B and an empty pycache prefix (bootstrap.py); -O is always refused (asserts must be live)."""
import ast, os, sys, hashlib, json, platform, sysconfig, ctypes
HERE = os.path.dirname(os.path.realpath(__file__))
PACKAGES = ('python-flint', 'numpy', 'scipy', 'scikit-fem', 'mpmath')
ENV_PREFIXES = ('OMP_', 'OPENBLAS_', 'VECLIB_', 'MKL_', 'FLINT', 'TM_', 'NPY_', 'PYTHON', 'GOTO', 'BLIS_', 'ACCELERATE', 'DARTS_')

def local_imports(f):
    t = ast.parse(open(os.path.join(HERE, f)).read())
    out = set()
    for n in ast.walk(t):
        if isinstance(n, ast.Import):
            for a in n.names: out.add(a.name.split('.')[0])
        elif isinstance(n, ast.ImportFrom) and n.module and n.level == 0:
            out.add(n.module.split('.')[0])
    return sorted(m + '.py' for m in out if os.path.exists(os.path.join(HERE, m + '.py')))

def closure(entries):
    seen = set(); st = list(entries)
    while st:
        f = st.pop()
        if f in seen: continue
        seen.add(f); st.extend(local_imports(f))
    return sorted(seen)

def sha_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''): h.update(chunk)
    return h.hexdigest()

def loaded_images():
    """paths of the Mach-O images loaded in this process (dyld), or of the mapped files on Linux"""
    if sys.platform == 'darwin':
        lib = ctypes.CDLL('/usr/lib/libSystem.B.dylib')
        lib._dyld_image_count.restype = ctypes.c_uint32
        lib._dyld_get_image_name.restype = ctypes.c_char_p; lib._dyld_get_image_name.argtypes = [ctypes.c_uint32]
        return sorted({lib._dyld_get_image_name(i).decode() for i in range(lib._dyld_image_count())})
    try:
        return sorted({l.split()[-1] for l in open('/proc/self/maps') if len(l.split()) >= 6 and l.split()[-1].startswith('/')})
    except OSError:
        return []

def _walk(root, suffixes=None):
    for dp, dn, fn in os.walk(root):
        dn[:] = [d for d in dn if d != '__pycache__']
        for f in fn:
            if suffixes is None or f.endswith(suffixes): yield os.path.join(dp, f)

SYSTEM_PREFIXES = ('/usr/lib/', '/System/')     # OS-owned images (dyld shared cache): bound through the OS build (uname)
def inventory_paths():
    """the set of runtime files to bind: interpreter, every non-OS image loaded so far, the standard library tree, every file of
    the numerical packages, the *.pth files and every top-level *.py module of the site directories (sitecustomize, usercustomize
    and the targets of *.pth import hooks, e.g. virtualenv's _virtualenv.py, which no installer RECORD lists)"""
    import importlib.metadata as md
    paths = set()
    paths.add(os.path.realpath(sys.executable))
    for p in loaded_images():
        if os.path.isfile(p) and not p.startswith(SYSTEM_PREFIXES): paths.add(os.path.realpath(p))
    sp = sysconfig.get_paths()
    for p in _walk(os.path.realpath(sp['stdlib']), ('.py', '.so')): paths.add(os.path.realpath(p))
    for pkg in PACKAGES:
        try: dist = md.distribution(pkg)
        except md.PackageNotFoundError: continue
        for f in dist.files or []:
            p = dist.locate_file(f)
            if os.path.isfile(p) and '__pycache__' not in str(f): paths.add(os.path.realpath(p))
    sitedirs = {sp['purelib'], sp['platlib']}
    try:
        import bootstrap; sitedirs.add(bootstrap.venv_site_packages())
    except SystemExit: pass
    for d in sitedirs:
        if os.path.isdir(d):
            for f in os.listdir(d):
                if f.endswith(('.pth', '.py')) and os.path.isfile(os.path.join(d, f)): paths.add(os.path.realpath(os.path.join(d, f)))
    return paths

def static_inventory(paths=None):
    """{realpath: sha256} of inventory_paths() (or of the given paths: the same key set re-hashed at the end of a run), read from the
    file system at EVERY call -- nothing is memoised"""
    if paths is None: paths = inventory_paths()
    out = {}
    for p in sorted(paths):
        try: out[p] = sha_file(p)
        except OSError: out[p] = None
    return out

def record_hashes():
    """{realpath: (distribution, urlsafe-b64 sha256)} from the installers' RECORD files of every distribution on sys.path"""
    import importlib.metadata as md
    out = {}
    for dist in md.distributions():
        for f in dist.files or []:
            if f.hash is not None and f.hash.mode == 'sha256':
                out[os.path.realpath(dist.locate_file(f))] = (dist.metadata['Name'], f.hash.value)
    return out
def _b64_sha256(path):
    import base64
    return base64.urlsafe_b64encode(hashlib.sha256(open(path, 'rb').read()).digest()).rstrip(b'=').decode()

def runtime_identity(inventory=None):
    import numpy as np, flint
    import fpenv, tm
    import importlib.metadata as md
    try: npcfg = np.show_config(mode='dicts')
    except Exception as e: npcfg = {'error': repr(e)}
    vers = {}
    for p in PACKAGES:
        try: vers[p] = md.version(p)
        except Exception: vers[p] = 'missing'
    inv = static_inventory() if inventory is None else inventory
    return dict(
        TM_ORDER=os.environ.get('TM_ORDER', '<unset>'), tm_order_effective=tm.ORDER, arb_prec=flint.ctx.prec,
        flint=dict(python_flint=getattr(flint, '__version__', None), FLINT=getattr(flint, '__FLINT_VERSION__', None),
                   ctx=repr(flint.ctx), threads=getattr(flint.ctx, 'threads', None)),
        python=sys.version, implementation=platform.python_implementation(), executable=sys.executable,
        flags=repr(sys.flags), xoptions=dict(sys._xoptions), pycache_prefix=sys.pycache_prefix,
        packages=vers, numpy_config=npcfg,
        os=dict(platform=platform.platform(), uname=platform.uname()._asdict(), mac_ver=platform.mac_ver()[0]),
        fp=fpenv.state(),
        env={k: v for k, v in sorted(os.environ.items()) if k.startswith(ENV_PREFIXES)},
        inventory_sha256=hashlib.sha256(json.dumps(inv, sort_keys=True).encode()).hexdigest(), inventory_count=len(inv))

def binding(script, manifest_sha256=None, family=None, inventory=None):
    """inventory: a static_inventory() computed by the caller in this same call sequence (fresh); None: computed here"""
    h = hashlib.sha256(); det = {}
    srcs = closure([script])
    det['sources'] = {f: sha_file(os.path.join(HERE, f)) for f in srcs}
    for f in srcs: h.update(b'SRC'); h.update(f.encode()); h.update(det['sources'][f].encode())
    rc = os.path.join(HERE, 'refcache'); det['refcache'] = {}
    if os.path.isdir(rc):
        for f in sorted(os.listdir(rc)):
            hb = sha_file(os.path.join(rc, f)); det['refcache'][f] = hb
            h.update(b'CACHE'); h.update(f.encode()); h.update(hb.encode())
    det['manifest_sha256'] = manifest_sha256; det['family'] = family
    h.update(b'MANIFEST'); h.update(str(manifest_sha256).encode()); h.update(b'FAMILY'); h.update(str(family).encode())
    inv = static_inventory() if inventory is None else inventory
    env = runtime_identity(inv); det['environment'] = env; det['inventory'] = inv
    h.update(b'RUNTIME'); h.update(json.dumps(env, sort_keys=True, default=str).encode())
    return h.hexdigest(), det

def loaded_modules_problems(inventory, sources, late=None):
    """every loaded module / image comes from an unchanged file of the START inventory, a pinned local source, an OS image, or --
    for a module first imported during the run from an installed distribution (e.g. numpy.f2py's optional charset_normalizer) --
    a file whose sha256 equals its installer RECORD entry (recorded in `late`); no cached bytecode was read"""
    bad = []; rec = None
    def late_ok(rp, what):
        nonlocal rec
        if rec is None: rec = record_hashes()
        r = rec.get(rp)
        if r is None: bad.append(f'{what}: {rp} is neither in the bound inventory nor in an installer RECORD'); return
        if _b64_sha256(rp) != r[1]: bad.append(f'{what}: {rp} differs from its installer RECORD ({r[0]})'); return
        if late is not None: late[rp] = r[0]
    for name, m in list(sys.modules.items()):
        spec = getattr(m, '__spec__', None)
        origin = getattr(spec, 'origin', None) if spec is not None else getattr(m, '__file__', None)
        if not origin or origin in ('built-in', 'frozen') or not os.path.isabs(str(origin)): continue
        rp = os.path.realpath(origin)
        if os.path.dirname(rp) == HERE and os.path.basename(rp) in sources:
            if sha_file(rp) != sources[os.path.basename(rp)]: bad.append(f'local module {name}: source changed')
        elif rp in inventory:
            if sha_file(rp) != inventory[rp]: bad.append(f'module {name}: file {rp} changed')
        else:
            late_ok(rp, f'module {name}')
        cached = getattr(spec, 'cached', None) if spec is not None else None
        if cached and sys.pycache_prefix and not os.path.realpath(cached).startswith(os.path.realpath(sys.pycache_prefix)):
            bad.append(f'module {name}: cached bytecode {cached} outside the empty pycache prefix')
    for p in loaded_images():
        if not os.path.isfile(p) or p.startswith(SYSTEM_PREFIXES): continue
        rp = os.path.realpath(p)
        if rp in inventory:
            if sha_file(rp) != inventory[rp]: bad.append(f'loaded image {p} changed')
        else:
            late_ok(rp, 'loaded image')
    if sys.pycache_prefix and os.path.isdir(sys.pycache_prefix):
        for dp, dn, fn in os.walk(sys.pycache_prefix):
            if fn: bad.append(f'bytecode was written under the pycache prefix ({dp})'); break
    return bad

def proof_flags_problems():
    """the interpreter flags required of every proof process (DARTS_PIN set): isolated, -B, empty pycache prefix, no -O"""
    bad = []
    f = sys.flags
    if f.optimize: bad.append('python -O / -OO: asserts are stripped; refused')
    if os.environ.get('DARTS_PIN'):
        if not (f.isolated and f.ignore_environment and f.no_user_site and f.safe_path): bad.append('not running isolated (-I)')
        if not f.no_site: bad.append('the site module ran (-S missing): *.pth startup hooks were executed')
        if not f.dont_write_bytecode: bad.append('bytecode writing enabled (-B missing)')
        if not sys.pycache_prefix: bad.append('no pycache prefix (-X pycache_prefix=<empty dir> missing)')
    return bad

if __name__ == '__main__':
    hx, det = binding(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None, sys.argv[3] if len(sys.argv) > 3 else None)
    env = det['environment']
    print(hx); print(json.dumps(dict(sources=len(det['sources']), family=det['family'], manifest_sha256=det['manifest_sha256'],
                                     inventory=env['inventory_count'], flags=env['flags'], fp=env['fp']['fpcr'], os=env['os']['platform']), indent=1))
