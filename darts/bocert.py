"""BO box certificate: lambda_1 >= E from the rigorous 1-D bound (bo1d.py, Lemma 11.1) on the outer domain of the box, and
lambda_2 < (5/2) E from the rigorous P2 Rayleigh-Ritz bound over the box (p2up.py).  Hence xi < 5/2 on the whole box.
  bd boxes   : params (a_x, a_y, c_x, c_y); P2 domain = the dart itself (BD chart); slicing axes: 9 candidate directions
  spike boxes: params (t, -, c_x, c_y); darts A = a (1, t) with a >= X0 (SPIKE2 P2 domain = dart with A = X0 (1, t)) or any a
               (COLLAR, P2 domain = CBD); vertical slicing with the spike tail."""
import math, numpy as np
from bo1d import bo_box, bo_spike
from p2up import solve_p2, p2_test

def cand_dirs(p):   # RIGOR-FUNCTION: heuristic: candidate slicing directions (any direction is valid)
    ax, ay, cx, cy = p
    bA = math.atan2(ay, ax); bC = math.atan2(cy, cx)
    return [0.0, math.pi/2, math.atan2(cy - ay, cx - ax), math.atan2(cy - ay, cx - ax) + math.pi/2, bA, bC,
            math.atan2(ay, ax - 1), math.atan2(cy, cx - 1), 0.5*(bA + bC)]

def bo_lower_bd(p0, r, K=120, dirs=None):
    p0 = np.asarray(p0, float); r = np.asarray(r, float)   # RIGOR: exact: float64 arrays of Python floats (no rounding)
    Abox = ((p0[0] - r[0], p0[1] - r[1]), (p0[0] + r[0], p0[1] + r[1]))
    Cbox = ((p0[2] - r[2], p0[3] - r[3]), (p0[2] + r[2], p0[3] + r[3]))
    best = (0.0, None)
    for th in (dirs if dirs is not None else cand_dirs(p0)):
        E, _ = bo_box(Abox, Cbox, th, K)
        if E > best[0]: best = (E, th)
    return best

def bo_cert_bd(gc, p0, r, K=120, sol=None):
    """gc: GCert of a BD-chart problem (only its P2 part is used)"""
    if sol is None: sol = solve_p2(gc, p0, r)
    lam2 = float(sol['lamp2'][1])   # RIGOR: heuristic: float pre-check only (it can only make the certificate fail)
    E, th = bo_lower_bd(p0, r, K)
    out = dict(method='BO', E=E, theta=th, K=K, lam2p2=lam2, xi_bo=lam2/E if E > 0 else float('inf'))   # RIGOR: info: logged ratio
    if not (E > 0 and lam2 < 2.5*E): out.update(ok=False, why='float: lam2 >= 2.5 E'); return out
    res = p2_test(gc, p0, r, sol, E)
    out.update(res)
    return out

def bo_cert_spike(gc, p0, r, K=160, sol=None):
    """gc: GCert of the SPIKE (COLLAR) or SPIKE2 problem (only its P2 part is used); the BO bound covers every a > 0"""
    p0 = np.asarray(p0, float); r = np.asarray(r, float)   # RIGOR: exact: float64 arrays of Python floats (no rounding)
    t_hi = p0[0] + r[0]
    Cbox = ((p0[2] - r[2], p0[3] - r[3]), (p0[2] + r[2], p0[3] + r[3]))
    if sol is None: sol = solve_p2(gc, p0, r)
    lam2 = float(sol['lamp2'][1])   # RIGOR: heuristic: float pre-check only (it can only make the certificate fail)
    E, info = bo_spike(t_hi, Cbox, K)
    out = dict(method='BO', E=E, theta=None, K=K, lam2p2=lam2, xi_bo=lam2/E if E > 0 else float('inf'))   # RIGOR: info: logged ratio
    if not (E > 0 and lam2 < 2.5*E): out.update(ok=False, why='float: lam2 >= 2.5 E'); return out
    res = p2_test(gc, p0, r, sol, E)
    out.update(res)
    return out
