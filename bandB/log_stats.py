"""log_stats.py — statistics of the band cover logs for BAND_README.md (reads the JSONL; no trust needed: the
replay re-derives everything).  Author: Anonymous.  usage: python log_stats.py runs/cover_chart1.jsonl runs/cover_chart2.jsonl"""
import sys, json, re
from fractions import Fraction as Fr
import numpy as np
for f in sys.argv[1:]:
    cnt = {}; margins = []; depth = []; nS = []; nR = []; secs = []; summ = None; vol = {}
    for ln in open(f):
        d = json.loads(ln)
        if "summary" in d: summ = d["summary"]; continue
        if "subtree" in d: continue
        st = d["status"]; cnt[st] = cnt.get(st, 0) + 1
        v = 1
        for a, b in zip(d["lo"], d["hi"]): v *= Fr(b) - Fr(a)
        vol[st] = vol.get(st, 0) + v
        depth.append(len(d["path"]))
        if st == "CERT":
            margins.append(float(d["margin"])); nS.append(d["nS"]); nR.append(d["nR"]); secs.append(d["sec"])
    tot = sum(vol.values())
    print(f"== {f}")
    print("  leaves:", sum(cnt.values()), " counts:", dict(sorted(cnt.items())))
    print("  volume fractions:", {k: f"{float(v / tot):.4f}" for k, v in sorted(vol.items())})
    m = np.array(margins)
    print(f"  CERT margins R2/(2.5E): max {m.max():.6f}, quantiles 50/90/99%: {np.percentile(m, 50):.4f} {np.percentile(m, 90):.4f} {np.percentile(m, 99):.4f}; min {m.min():.4f}")
    print(f"  depth: min {min(depth)} median {int(np.median(depth))} max {max(depth)};  nS {min(nS)}-{max(nS)}, nR {min(nR)}-{max(nR)};  cert seconds total {sum(secs):.1f}, mean {np.mean(secs)*1e3:.1f} ms")
    if summ: print("  summary line:", json.dumps(summ))
