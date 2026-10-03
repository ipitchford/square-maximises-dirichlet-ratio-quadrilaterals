"""Audit (coordinator, 2026-10-01): did any logged FARL FEM certificate pass through the NaN path of the old FARL tail bound
(v = (1 - k0 u^2)**2 + (k0 u)**2 on balls; ball**2 is NaN for a ball centred at 0, and max() dropped it)?  For every FEM leaf of
the FARL families (FP, PF, and PP boxes with one long arm) the old 64-cell loop is re-evaluated and NaN cells are counted.
Usage: audit_farl_nan.py"""
import sys, os, glob, json, re, numpy as np
from fractions import Fraction
HERE = '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts'
sys.path.insert(0, HERE)
from flint import arb
import fam5
tot = farl = nan_leaves = nan_cells = 0
for fam_name in ('FP', 'PF', 'PP'):
    fam = fam5.FAMILIES[fam_name]()
    for f in sorted(glob.glob(os.path.join(HERE, 'certs3', fam_name.lower() + '.*.jsonl'))):
        if f.endswith('.merged.jsonl'): continue
        for line in open(f):
            if '"FEM"' not in line: continue
            d = json.loads(line)
            if d.get('method') != 'FEM': continue
            tot += 1
            st = fam.fem_setup(np.array(d['p0']), np.array(d['r']), d['k']) if hasattr(fam, 'fem_setup') else None
            if st is None or not st[0].name.startswith('FARL'): continue
            farl += 1
            k0 = Fraction(float(re.search(r'_k0([0-9.e-]+)$', st[0].name).group(1)))
            kk = arb(k0.numerator)/arb(k0.denominator)
            p0, r = d['pc'], d['rc']; bad = 0
            for i in range(64):
                lo = p0[0] - r[0] + 2*r[0]*i/64; hi = p0[0] - r[0] + 2*r[0]*(i + 1)/64
                uu = arb(lo).union(arb(hi))
                v = (1 - kk*uu*uu)**2 + (kk*uu)**2              # the OLD formula
                if not v.is_finite(): bad += 1
            if bad: nan_leaves += 1; nan_cells += bad
print(json.dumps(dict(fem_leaves=tot, farl_leaves=farl, farl_leaves_with_nan_cells=nan_leaves, nan_cells=nan_cells)))
