"""Floating-point environment guard (adversarial review round 2, SERIOUS 4).
Every float error model of the verifier (ldl.py, gcert.py, gsigma.py: gamma_n with u = 2^-53) assumes IEEE binary64 with
round-to-nearest-even and gradual underflow.  check() raises FPEnvError unless
  (a) behavioural tests, evaluated at RUN time on non-constant operands (no compile-time folding), in Python floats AND in
      NumPy ufunc loops, show round-to-nearest-even (ties to even, both signs, below / above half an ulp) and gradual underflow
      (no flush-to-zero of results, no flush-to-zero / denormals-are-zero of inputs);
  (b) on Darwin arm64, the control register FPCR read through fegetenv() has RMode = 00 (nearest), FZ = 0 (bit 24), FIZ = 0
      (bit 0), AH = 0 (bit 1), FZ16 = 0 (bit 19), and fegetround() returns FE_TONEAREST.
The replays call check() at start, after every leaf, and at the end; state() is part of the replay binding."""
import ctypes, ctypes.util, platform, sys

class FPEnvError(ArithmeticError): pass

_libc = None
def _lib():
    global _libc
    if _libc is None:
        name = ctypes.util.find_library('c') or ('libSystem.dylib' if sys.platform == 'darwin' else None)
        _libc = ctypes.CDLL(name) if name else None
    return _libc

class _FenvArm64Darwin(ctypes.Structure):
    _fields_ = [('fpsr', ctypes.c_ulonglong), ('fpcr', ctypes.c_ulonglong)]

def _fpcr():
    """the AArch64 FPCR on Darwin arm64 (via fegetenv), else None"""
    if not (sys.platform == 'darwin' and platform.machine() == 'arm64'): return None
    lib = _lib(); env = _FenvArm64Darwin()
    if lib.fegetenv(ctypes.byref(env)) != 0: raise FPEnvError('fegetenv failed')
    return int(env.fpcr)

def _behaviour():
    """run-time tests; returns a list of failures (empty if round-to-nearest-even with gradual underflow)"""
    import numpy as np
    one = float('1'); h = float.fromhex('0x1p-53'); q = float.fromhex('0x1p-54'); tq = float.fromhex('0x3p-54')   # RIGOR: exact: exact binary constants parsed at run time (the test operands)
    tiny = float.fromhex('0x1p-1074'); mn = float.fromhex('0x1p-1022'); half = float('0.5')   # RIGOR: exact: exact binary constants parsed at run time
    bad = []
    if one + h != one: bad.append('1 + 2^-53 is not 1 (ties-to-even / upward)')
    if -one - h != -one: bad.append('-1 - 2^-53 is not -1 (ties-to-even / downward)')
    if one + tq != float.fromhex('0x1.0000000000001p+0'): bad.append('1 + 3*2^-54 does not round up (toward zero / downward)')
    if -one - tq != float.fromhex('-0x1.0000000000001p+0'): bad.append('-1 - 3*2^-54 does not round down (toward zero / upward)')
    if one + q != one: bad.append('1 + 2^-54 is not 1 (upward)')
    if float.fromhex('0x1.0000000000001p+0') + h != float.fromhex('0x1.0000000000002p+0'): bad.append('odd + half ulp does not round to even')
    if mn*half != float.fromhex('0x1p-1023'): bad.append('2^-1022 * 0.5 is flushed (FTZ)')
    if tiny + tiny != float.fromhex('0x1p-1073'): bad.append('subnormal inputs are flushed (DAZ / FIZ)')
    if not (tiny/float('2') == 0.0 and tiny*float('0.75') == tiny): bad.append('subnormal rounding is not to nearest')   # RIGOR: exact: exact constants; the test itself
    a = np.array([one, -one, one, mn, tiny]); b = np.array([h, -h, tq, 0.0, tiny])   # RIGOR: exact: the NumPy test operands
    s = a + b
    if not (s[0] == one and s[1] == -one and s[2] == float.fromhex('0x1.0000000000001p+0') and s[4] == float.fromhex('0x1p-1073')):
        bad.append('NumPy ufunc rounding differs from round-to-nearest / gradual underflow')
    if np.multiply(np.array([mn]), np.array([half]))[0] != float.fromhex('0x1p-1023'): bad.append('NumPy flushes subnormal results')   # RIGOR: exact: the NumPy test itself
    return bad

def state():
    """a description of the floating-point environment (bound into every replay)"""
    pc = _fpcr()
    lib = _lib()
    rnd = int(lib.fegetround()) if lib is not None else None
    return dict(fegetround=rnd, fpcr=pc, behaviour_failures=_behaviour(), float_info=repr(sys.float_info))

def check():
    bad = _behaviour()
    pc = _fpcr()
    if pc is not None:
        if (pc >> 22) & 3: bad.append(f'FPCR.RMode = {(pc >> 22) & 3} (not nearest)')
        if (pc >> 24) & 1: bad.append('FPCR.FZ = 1 (flush to zero)')
        if pc & 1: bad.append('FPCR.FIZ = 1 (flush inputs to zero)')
        if (pc >> 1) & 1: bad.append('FPCR.AH = 1 (alternate handling)')
        if (pc >> 19) & 1: bad.append('FPCR.FZ16 = 1')
        lib = _lib()
        if int(lib.fegetround()) != 0: bad.append(f'fegetround() = {int(lib.fegetround()):#x}, not FE_TONEAREST')
    if bad: raise FPEnvError('floating-point environment is not IEEE round-to-nearest with gradual underflow: ' + '; '.join(bad))
    return True

if __name__ == '__main__':
    print(state()); print('check:', check())
