"""Stage 2: Taylor coefficients (arb) of K_kl, M_kl, LL_kl for all mode pairs; pickled to disk."""
import sys, time, pickle
from fractions import Fraction as Fr
from flint import arb, arb_mat
from engine2 import *

def run(N, m, mLL, rho0, out):
    modes = [(i,j) for i in range(1,N+1) for j in range(1,N+1)]
    base = [(1,1),(2,1),(1,2)]
    T = build_tensors(m, rho0, mLL=mLL)
    kmax = 2*N+2; amax = 0
    for tm in list(T['TK'].values()) + list(T['TLL'].values()) + [T['J']]:
        amax = max(amax, tm.poly.degrees()[0], tm.poly.degrees()[1])
    Ic, Is = trig_tables(int(amax)+2, kmax)
    def order(k, l, mm):
        nb = (k in base) + (l in base)
        return mm - (2 - nb)
    # contractors per tensor and truncation order
    cache = {}
    def contractor(name, key, mm):
        ck = (name, key, mm)
        if ck not in cache:
            tm = {'TK': T['TK'], 'TLL': T['TLL'], 'J': {None: T['J']}}[name][key]
            cache[ck] = Contractor(tm, mm)
        return cache[ck]
    t0 = time.time(); res = {'K':{}, 'M':{}, 'LL':{}, 'modes':modes, 'base':base, 'm':m, 'mLL':mLL, 'rho0':rho0, 'N':N}
    def addto(dst, src, rem):
        d, r = dst
        for g, val in src.items(): d[g] = d.get(g, arb(0)) + val
        return (d, r + rem)
    npairs = 0
    for ik, k in enumerate(modes):
        for l in modes[ik:]:
            mk = order(k, l, m); mL = order(k, l, mLL)
            # K
            acc = ({}, arb(0))
            for (al,be), (dk, dl) in {(0,0):(0,0), (1,1):(1,1)}.items():
                acc = addto(acc, *contractor('TK', (al,be), mk).contract(Ic, Is, k, dk, l, dl))
            o1 = contractor('TK', (0,1), mk).contract(Ic, Is, k, 0, l, 1); o2 = contractor('TK', (0,1), mk).contract(Ic, Is, k, 1, l, 0)
            acc = addto(acc, *o1); acc = addto(acc, *o2)
            res['K'][(k,l)] = acc
            # M
            res['M'][(k,l)] = contractor('J', None, mk).contract(Ic, Is, k, None, l, None)
            # LL
            acc = ({}, arb(0))
            for i in range(5):
                for j in range(5):
                    key = (min(i,j), max(i,j))
                    acc = addto(acc, *contractor('TLL', key, mL).contract(Ic, Is, k, i, l, j))
            res['LL'][(k,l)] = acc
            npairs += 1
            if npairs % 200 == 0: print(f"  {npairs} pairs, {time.time()-t0:.0f}s", flush=True)
    print(f"assembled {npairs} pairs in {time.time()-t0:.0f}s", flush=True)
    # convert arb to strings for pickling (arb not picklable): store (mid, rad) as str
    def ser(d):
        return {g: (str(v.mid()), str(v.rad())) for g, v in d.items()}
    out_res = {kk: {pair: (ser(v[0]), (str(v[1].mid()), str(v[1].rad()))) for pair, v in res[kk].items()} for kk in ('K','M','LL')}
    out_res.update(modes=modes, base=base, m=m, mLL=mLL, rho0=str(rho0), N=N)
    pickle.dump(out_res, open(out, 'wb'))
    print("saved", out)

if __name__ == "__main__":
    N = int(sys.argv[1]); m = int(sys.argv[2]); mLL = int(sys.argv[3]); rho0 = Fr(sys.argv[4]); out = sys.argv[5]
    run(N, m, mLL, rho0, out)
