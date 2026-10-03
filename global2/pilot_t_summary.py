"""Summary of a transport pilot directory: per-cell leaves, certificates, capped cells, extrapolation."""
import json, sys, os, re, glob, statistics as st
d = sys.argv[1]; ncells = int(sys.argv[2]); cap = int(sys.argv[3])
rows = []
for f in sorted(glob.glob(os.path.join(d, 'cell_*.jsonl')), key=lambda x: int(re.findall(r'cell_(\d+)', x)[0])):
    S = {}; npass = 0; t = 0.0
    lines = open(f).read().splitlines()
    for l in lines:
        r = json.loads(l); s = r['status']; S[s] = S.get(s, 0) + 1
        if s in ('pilot_pass', 'certified'): npass += 1; t += r.get('t', 0)
    rows.append((f, len(lines), npass, S, len(lines) >= cap, t))
K = len(rows); capped = sum(1 for r in rows if r[4])
passes = [r[2] for r in rows]; proc = [r[1] for r in rows]
print(f"cells sampled {K}/{ncells}; capped {capped}; mean processed {st.mean(proc):.0f}; mean certificates {st.mean(passes):.0f} (median {st.median(passes):.0f}, max {max(passes)})")
print(f"extrapolated: processed {ncells*st.mean(proc):.3g}, certificates {ncells*st.mean(passes):.3g} (capped cells counted at the cap: lower bound)")
tot_t = sum(r[5] for r in rows); tot_p = sum(passes)
if tot_p: print(f"mean Float64 time per certificate (incl. its share of nothing else): {tot_t/tot_p:.3f} s")
