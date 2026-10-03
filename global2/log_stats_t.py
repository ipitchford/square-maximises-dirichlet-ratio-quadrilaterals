"""Statistics of a transport cover log (JSONL from cover_t.jl): leaves by status, certificate details, volumes, times."""
import json, sys, collections, statistics, re
f = sys.argv[1]
st = collections.Counter(); vol = collections.Counter(); C = collections.defaultdict(collections.Counter)
kap = []; ts = []; tF = []; vols = []; worst = None; tleaf = 0.0; fails = []; ctree = collections.Counter(); cdepth = []
for ln in open(f):
    try:
        r = json.loads(ln)
    except Exception:
        continue
    s = r['status']; st[s] += 1
    if s == 'split':
        continue
    v = 1.0
    for i in range(4): v *= r['hi'][i] - r['lo'][i]
    vol[s] += v; tleaf += r.get('t', 0.0)
    if s in ('FAIL', 'ERROR'): fails.append(r)
    if s == 'COMPOUND':
        for m in re.findall(r'E[0-9]T?', r['tree']): ctree[m] += 1
        d = 0; dm = 0
        for ch in r['tree']:
            d += ch == '('; d -= ch == ')'; dm = max(dm, d)
        cdepth.append(dm)
    if s == 'certified':
        for k in ('mode', 'bound', 'fo', 'n', 'shift', 'which'):
            C[k][str(r.get(k))] += 1
        b = r['bound']; C['bound family'][b.split(':')[0] if ':' in b else ('const/boot' if b in ('const', 'boot') else re.sub(r'\d+$', '', b))] += 1
        k = r['kappa_cert']; kap.append(k); ts.append(r['t']); tF.append(r.get('tF', 0.0)); vols.append(v)
        if worst is None or (k, r['kappa_f64']) > (worst['kappa_cert'], worst['kappa_f64']): worst = r
print("records by status:", dict(st))
print("leaves:", sum(n for s, n in st.items() if s != 'split'), " splits:", st['split'])
print("volume by status:", {s: float('%.4e' % v) for s, v in vol.items()})
if kap:
    print("certified %d:" % len(kap))
    for k in ('mode', 'bound family', 'bound', 'fo', 'n', 'shift', 'which'):
        print("   %-12s %s" % (k, dict(C[k].most_common())))
    print("logged kappa_cert (ladder value; refine with margin_refine_t.jl): max %.6f, #(=2.5 fallback) %d, median %.4f, min %.4f"
          % (max(kap), sum(1 for k in kap if k >= 2.5), statistics.median(kap), min(kap)))
    qs = statistics.quantiles(kap, n=20)
    print("   kappa_cert 5/25/50/75/95%%: %.3f %.3f %.3f %.3f %.3f" % (qs[0], qs[4], qs[9], qs[14], qs[18]))
    print("box volume: median %.3e, min %.3e, max %.3e; total certified volume %.4e" % (statistics.median(vols), min(vols), max(vols), sum(vols)))
    print("time per certified box (one thread, incl. screen): mean %.3f s (Float64 part %.3f s); median %.3f s; sum over all leaves %.0f thread-s"
          % (statistics.mean(ts), statistics.mean(tF), statistics.median(ts), tleaf))
    print("largest logged kappa_cert box:", worst['path'], "kappa_cert", worst['kappa_cert'], "kappa_f64", worst['kappa_f64'], "bound", worst['bound'], worst['mode'])
if ctree:
    print("COMPOUND sub-leaves by condition:", dict(ctree.most_common()), " max depth", max(cdepth), " mean depth %.2f" % statistics.mean(cdepth))
print("FAIL/ERROR leaves:", len(fails))
for r in fails[:20]:
    print("  ", r['status'], r['path'], r.get('why', ''))
