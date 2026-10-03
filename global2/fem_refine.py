"""FEM refinement for the near-square rows where theta2 < FEM(k=6) lambda2: show FEM(k) decreases below theta2."""
import sys, csv
sys.path.insert(0, '../explore')
from quadfem import eigs
rows = list(csv.DictReader(open('validation_points.csv')))
out = open('validation_fem_refine.log', 'w')
seen = set()
for r in rows:
    if r['tag'] not in ('nearsq', 'random') or r['n'] != '12': continue
    V = [(float(r[f'v{j}x']), float(r[f'v{j}y'])) for j in range(1,5)]
    if tuple(V) in seen: continue
    seen.add(tuple(V))
    t2 = float(r['theta2']); t1 = float(r['theta1'])
    ls = {k: eigs(V, k=k, nev=3) for k in (6, 7, 8)}
    s = f"{r['tag']} p=({float(r['a']):.3f},{float(r['b']):.3f},{float(r['c']):.3f},{float(r['d']):.3f}) theta1={t1:.9f} FEM l1 k6/7/8 = {ls[6][0]:.9f} {ls[7][0]:.9f} {ls[8][0]:.9f} | theta2={t2:.9f} FEM l2 k6/7/8 = {ls[6][1]:.9f} {ls[7][1]:.9f} {ls[8][1]:.9f}  theta2<=FEM8? {t2 <= ls[8][1]}"
    print(s); out.write(s + '\n')
