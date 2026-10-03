"""Root lists ('manifests') of the four runs of the first session, regenerated from the root definitions of their run scripts
(run_crit2.py COLLAR/SPIKE with XS = 1, run_main_cbox.py with AYMAX = 1.5479, run_box.py shallow box with 4 roots per axis),
for the fail-closed replay2.py.  Usage: manifests_legacy.py OUTDIR"""
# RIGOR-MODULE: heuristic: definitions of the legacy root lists (any floats); accounting.py checks the manifests against them and the coverage exactly
import sys, os, json, itertools, numpy as np
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from regions import T1, T2
CX = (0.40, 0.60); CY = (-0.966, -0.766); nC = 2
cxs = np.linspace(CX[0], CX[1], nC + 1); cys = np.linspace(CY[0], CY[1], nC + 1)
def collar():
    return [(f'COL_{i}{j}', [T1, 0.0, (cxs[i]+cxs[i+1])/2, (cys[j]+cys[j+1])/2], [0.0, 0.0, (cxs[1]-cxs[0])/2, (cys[1]-cys[0])/2])
            for i in range(nC) for j in range(nC)]
def spike():
    ts = [T1]
    while ts[-1] < T2: ts.append(min(ts[-1]*1.25, T2))
    return [(f'SPK_{a}_{i}{j}', [(ts[a]+ts[a+1])/2, 0.0, (cxs[i]+cxs[i+1])/2, (cys[j]+cys[j+1])/2],
             [(ts[a+1]-ts[a])/2, 0.0, (cxs[1]-cxs[0])/2, (cys[1]-cys[0])/2]) for a in range(len(ts) - 1) for i in range(nC) for j in range(nC)]
def main_cbox(AYMAX=1.5479):
    roots = []; ays = [0.25]
    while ays[-1] < AYMAX: ays.append(min(ays[-1]*1.2, AYMAX))
    for a in range(len(ays) - 1):
        y0, y1 = ays[a], ays[a+1]; xmax = y1/0.3 + 0.05; xs = [1.0]
        while xs[-1] < xmax: xs.append(xs[-1] + 0.2*max(1.0, y1))
        for b in range(len(xs) - 1):
            roots.append((f'M_{a}_{b}', [(xs[b]+xs[b+1])/2, (y0+y1)/2, (CX[0]+CX[1])/2, (CY[0]+CY[1])/2],
                          [(xs[b+1]-xs[b])/2, (y1-y0)/2, (CX[1]-CX[0])/2, (CY[1]-CY[0])/2]))
    return roots
def shallow(lim=(0.965, 1.262, 0.085, 0.796, 0.479, 1.213, -1.091, -0.483), nr=4):
    edges = [np.linspace(lim[2*i], lim[2*i+1], nr + 1) for i in range(4)]
    out = []
    for idx in itertools.product(range(nr), repeat=4):
        p0 = [(edges[i][idx[i]] + edges[i][idx[i]+1])/2 for i in range(4)]; r = [(edges[i][1] - edges[i][0])/2 for i in range(4)]
        out.append(('B_' + '_'.join(map(str, idx)), [float(x) for x in p0], [float(x) for x in r]))
    return out
if __name__ == '__main__':
    d = sys.argv[1]
    for name, fam, roots in (('crit_collar', 'COLLAR', collar()), ('crit_spike_X1', 'SPIKE_X1', spike()), ('main_cbox', 'MAIN_CBOX', main_cbox()),
                             ('shallow_crit', 'BOX', shallow())):
        json.dump(dict(family=fam, roots=[(a, [float(x) for x in b], [float(x) for x in c]) for a, b, c in roots]), open(f'{d}/{name}.json', 'w'))
        print(name, len(roots))

# ---------------------------------------------------------------------------------------------- nominal cells (for accounting.py)
# (root id, lo[4], hi[4]): the closed cell each root was meant to cover, with the SAME float edges as the generators above
# (adjacent cells share their float edge exactly); inactive coordinates have lo = hi.
def collar_cells():
    return [(f'COL_{i}{j}', [T1, 0.0, cxs[i], cys[j]], [T1, 0.0, cxs[i + 1], cys[j + 1]]) for i in range(nC) for j in range(nC)]
def spike_cells():
    ts = [T1]
    while ts[-1] < T2: ts.append(min(ts[-1]*1.25, T2))
    return [(f'SPK_{a}_{i}{j}', [ts[a], 0.0, cxs[i], cys[j]], [ts[a + 1], 0.0, cxs[i + 1], cys[j + 1]])
            for a in range(len(ts) - 1) for i in range(nC) for j in range(nC)]
def main_cbox_cells(AYMAX=1.5479):
    out = []; ays = [0.25]
    while ays[-1] < AYMAX: ays.append(min(ays[-1]*1.2, AYMAX))
    for a in range(len(ays) - 1):
        y0, y1 = ays[a], ays[a+1]; xmax = y1/0.3 + 0.05; xs = [1.0]
        while xs[-1] < xmax: xs.append(xs[-1] + 0.2*max(1.0, y1))
        for b in range(len(xs) - 1):
            out.append((f'M_{a}_{b}', [xs[b], y0, CX[0], CY[0]], [xs[b + 1], y1, CX[1], CY[1]]))
    return out
def shallow_cells(lim=(0.965, 1.262, 0.085, 0.796, 0.479, 1.213, -1.091, -0.483), nr=4):
    edges = [np.linspace(lim[2*i], lim[2*i+1], nr + 1) for i in range(4)]
    return [('B_' + '_'.join(map(str, idx)), [float(edges[i][idx[i]]) for i in range(4)], [float(edges[i][idx[i] + 1]) for i in range(4)])
            for idx in itertools.product(range(nr), repeat=4)]
