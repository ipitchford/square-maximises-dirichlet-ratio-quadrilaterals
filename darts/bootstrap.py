"""Proof-run bootstrap (adversarial review round 3, item 6).  Every proof process is started as
    <venv python, absolute> -I -S -B -X pycache_prefix=<empty run directory> <darts>/bootstrap.py <darts>/SCRIPT ARGS...
-I: isolated (no PYTHON* variables, no user site, no script directory on sys.path); -S: the site module is not run, so no *.pth
startup hooks and no sitecustomize / usercustomize; -B and the empty pycache prefix: no bytecode is read or written, every module
is compiled from its bound source.  Without the site module the virtual environment is not activated, so this file puts its own
directory and the venv's site-packages (<venv>/lib/pythonX.Y/site-packages, the venv being the directory of sys.executable's
pyvenv.cfg) on sys.path, and runs SCRIPT (which must lie in this directory) as __main__."""
import sys, os, runpy
HERE = os.path.dirname(os.path.realpath(__file__))
def venv_site_packages():
    root = os.path.dirname(os.path.dirname(os.path.abspath(sys.executable)))
    if not os.path.isfile(os.path.join(root, 'pyvenv.cfg')): raise SystemExit(f'bootstrap: {sys.executable} is not a venv interpreter')
    sp = os.path.join(root, 'lib', f'python{sys.version_info.major}.{sys.version_info.minor}', 'site-packages')
    if not os.path.isdir(sp): raise SystemExit(f'bootstrap: no site-packages at {sp}')
    return sp
if __name__ == '__main__':
    if len(sys.argv) < 2: raise SystemExit('usage: bootstrap.py SCRIPT [ARGS...]')
    script = os.path.realpath(sys.argv[1])
    if os.path.dirname(script) != HERE or not os.path.isfile(script): raise SystemExit(f'bootstrap: {sys.argv[1]!r} is not a script of {HERE}')
    sys.argv = [script] + sys.argv[2:]
    sys.path[:0] = [HERE]
    sys.path.append(venv_site_packages())
    # load the standard-library extension modules (and their external libraries, e.g. libmpdec for _decimal) that the verifier
    # uses, before any runtime inventory is taken
    import decimal, fractions, zlib, hashlib, ctypes, json, pickle, struct, math   # noqa: F401
    runpy.run_path(script, run_name='__main__')
