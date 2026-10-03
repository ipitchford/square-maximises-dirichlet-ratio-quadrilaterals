"""corrupt_log.py — make a deliberately corrupted copy of a cover log (fail-closed tests of replay.jl).
Author: Anonymous.  usage: python corrupt_log.py SRC DST CASE"""
import sys, re, json
from fractions import Fraction as Fr
src, dst, case = sys.argv[1:4]
L = open(src).read().splitlines()
leaf = [i for i, l in enumerate(L) if l.startswith('{"chart"')]
cert = [i for i in leaf if '"status":"CERT"' in L[i]]
lema = [i for i in leaf if '"status":"LEMMA_A"' in L[i]]
mark = [i for i, l in enumerate(L) if l.startswith('{"subtree"')]
summ = [i for i, l in enumerate(L) if l.startswith('{"summary"')]
path = lambda l: re.search(r'"path":"([A-Da-d]*)"', l).group(1)
def bump_summary(dcert=0, dleaves=0, dstat=None):
    s = L[summ[0]]
    m = re.search(r'"leaves":(\d+)', s); s = s.replace(m.group(0), f'"leaves":{int(m.group(1)) + dleaves}', 1)
    m = re.search(r'"CERT":(\d+)', s); s = s.replace(m.group(0), f'"CERT":{int(m.group(1)) + dcert}', 1)
    if dstat:
        k, dv = dstat; m = re.search(rf'"{k}":(\d+)', s); s = s.replace(m.group(0), f'"{k}":{int(m.group(1)) + dv}', 1)
    L[summ[0]] = s
def bump_marker(p, dv):
    for i in mark:
        m = re.match(r'\{"subtree":"([A-Da-d]*)","done":true,"leaves":(\d+),', L[i])
        if p.startswith(m.group(1)):
            L[i] = L[i].replace(f'"leaves":{m.group(2)},', f'"leaves":{int(m.group(2)) + dv},', 1); return
if case == "ok": pass
elif case == "unknown_status": L[lema[0]] = L[lema[0]].replace('"status":"LEMMA_A"', '"status":"LEMMA_Z"')
elif case == "malformed": L[cert[0]] = L[cert[0]][:-1]
elif case == "extra_field": L[cert[0]] = L[cert[0]][:-1] + ',"extra":1}'
elif case == "dup_line": L.insert(cert[0] + 1, L[cert[0]])
elif case == "dup_path": L.insert(cert[0] + 1, re.sub(r'"sec":[0-9.e-]+', '"sec":9.9999', L[cert[0]]))
elif case == "missing_leaf": del L[cert[3]]
elif case == "missing_marker": del L[mark[5]]
elif case == "fail_floor":
    l = L[cert[2]]; box = re.search(r'"lo":\[[^\]]*\],"hi":\[[^\]]*\]', l).group(0)
    L[cert[2]] = '{"chart":%s,"path":"%s","status":"FAIL_FLOOR",%s,"why":"forged","dim":1}' % (re.search(r'"chart":(\d)', l).group(1), path(l), box)
elif case == "wrong_box":        # canonical but wrong value: lo[0] + 2^-30
    l = L[cert[4]]; m = re.search(r'"lo":\["(\d+)/(\d+)"', l)
    q = Fr(int(m.group(1)), int(m.group(2))) + Fr(1, 2**30)
    L[cert[4]] = l.replace(m.group(0), f'"lo":["{q.numerator}/{q.denominator}"', 1)
elif case == "noncanonical":
    l = L[cert[4]]; m = re.search(r'"lo":\["(\d+)/(\d+)"', l)
    L[cert[4]] = l.replace(m.group(0), f'"lo":["{2 * int(m.group(1))}/{2 * int(m.group(2))}"', 1)
elif case == "bad_rule": L[lema[0]] = re.sub(r'7/200', '7/100', L[lema[0]], count=1)
elif case == "false_exclusion":   # a CERT leaf re-labelled as a whitelisted but false Lemma-A exclusion; counts kept consistent
    l = L[cert[6]]; box = re.search(r'"lo":\[[^\]]*\],"hi":\[[^\]]*\]', l).group(0); ch = re.search(r'"chart":(\d)', l).group(1)
    L[cert[6]] = '{"chart":%s,"path":"%s","status":"LEMMA_A","rule":"c0 <= 7/200 c2",%s}' % (ch, path(l), box)
    bump_summary(dcert=-1, dstat=("LEMMA_A", +1))
elif case == "forged_parent_cert":  # two sibling CERT leaves replaced by ONE CERT record of their parent; markers/summary consistent
    P = {path(L[i]): i for i in cert}
    for p, i in P.items():
        last = p[-1]
        if last.islower() and (p[:-1] + last.upper()) in P and len(p) > 9:
            j = P[p[:-1] + last.upper()]; break
    lo = re.search(r'"lo":\[[^\]]*\]', L[i]).group(0); hi = re.search(r'"hi":\[[^\]]*\]', L[j]).group(0)
    new = re.sub(r'"path":"[A-Da-d]*"', f'"path":"{p[:-1]}"', L[i], count=1)
    new = re.sub(r'"hi":\[[^\]]*\]', hi, new, count=1)
    L[i] = new
    bump_marker(p, -1); bump_summary(dcert=-1, dleaves=-1)   # before the deletion: indices still valid
    del L[j]
elif case == "summary_mismatch": bump_summary(dcert=+1)
elif case == "summary_extra_field": L[summ[0]] = L[summ[0]][:-2] + ',"extra":1}}'
elif case == "summary_dup_key":
    m = re.search(r'"CERT":(\d+)', L[summ[0]]); L[summ[0]] = L[summ[0]].replace(m.group(0), m.group(0) + "," + m.group(0), 1)
elif case == "summary_bad_subtrees": L[summ[0]] = re.sub(r'"subtrees":\d+', '"subtrees":255', L[summ[0]], count=1)
elif case == "summary_bad_depth": L[summ[0]] = re.sub(r'"max_depth":(\d+)', lambda m: '"max_depth":%d' % (int(m.group(1)) + 1), L[summ[0]], count=1)
elif case == "missing_summary": del L[summ[0]]
elif case == "mixed_chart": L[cert[1]] = re.sub(r'"chart":(\d)', lambda m: '"chart":%d' % (3 - int(m.group(1))), L[cert[1]], count=1)
elif case == "trailing_garbage": L.append("garbage")
elif case == "empty": L = []
else: raise SystemExit("unknown case " + case)
open(dst, "w").write("\n".join(L) + ("\n" if L else ""))
