"""Pilot of the rigorous box cover on representative core regions (BD chart).  Usage: pilot.py NAME"""
import sys, time, json, numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts')
from cover import Cover
REG = {
 'mod':     ((1.43, 0.82, 1.43, -0.82), (0.10, 0.10, 0.10, 0.10)),     # symmetric moderate dart (equilateral T, beta 0.3)
 'mod2':    ((1.72, 0.69, 2.41, -1.03), (0.10, 0.08, 0.12, 0.10)),     # (90,45,0.3,0.5)
 'crit':    ((2.40, 0.10, 0.55, -0.865), (0.10, 0.02, 0.05, 0.05)),    # spike, CBD near equilateral (xi ~ 2.33)
 'shallow': ((1.05, 0.61, 1.05, -0.61), (0.05, 0.05, 0.05, 0.05)),     # near Prop 5.1 (beta ~ 0.05)
}
name = sys.argv[1]; kmin = int(sys.argv[2]) if len(sys.argv) > 2 else 4
p0, r = REG[name]
cv = Cover(f'/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts/pilot_runs/pilot_{name}.jsonl', kmin=kmin, kmax=6)
t = time.time(); ok = cv.process_root(name, p0, r)
st = cv.stats; st['wall_s'] = time.time() - t; st['root_volume'] = float(np.prod(2*np.array(r)))
print(json.dumps(dict(name=name, ok=ok, **st)), flush=True)
