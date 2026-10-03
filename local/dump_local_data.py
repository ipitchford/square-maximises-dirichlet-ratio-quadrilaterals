"""Readable, exact records of the local certificate's data (referee report 2026-10-02, item 4.3):
  corrections.json  - every correction coefficient C_{jki} of the transported test functions, as an exact dyadic rational;
  remainders.json   - the Taylor-model remainder radii of every assembled K, M, LL entry used by the 3x3 pencil and of det S and
                      tr S for every shift level beta (upper bounds, as exact rationals and decimals).
The per-box sign enclosures are in the cover output (schur_cover2.json: for each of the 960 direction boxes and each s-interval the
shift level, the criterion used and the exact rational upper bound of det S/s^2 or tr S/s^2).
Usage: dump_local_data.py COEFFS.pkl SCHUR.pkl OUTDIR"""
import sys, json, pickle
from fractions import Fraction as Fr
from flint import arb
from localcert import load, corrections

def up(x):
    m, e = arb(x).upper().man_exp(); return Fr(int(m)) * Fr(2)**int(e)

def main(coeffs_fn, schur_fn, out):
    D = load(coeffs_fn)
    C = corrections(D)
    ent = [dict(j=list(j), k=list(k), i="abcd"[i], value=str(v)) for (j, k, i), v in sorted(C.items())]
    json.dump(dict(description="correction coefficients C_{jki} (exact dyadic rationals, denominator 2^40) of the test functions "
                   "w_j = (g_j + sum_k sum_i C_{jki} p_i g_k) o Phi_p^{-1}, j in {(1,1),(2,1),(1,2)}",
                   N=D['N'], count=len(ent), entries=ent), open(f"{out}/corrections.json", "w"), indent=0)
    base = D['base']
    asm = {}
    for name in ("K", "M", "LL"):
        rows = {}
        for a in base:
            for b in base:
                key = (a, b) if (a, b) in D[name] else (b, a)
                if key not in D[name]: continue
                cd, r = D[name][key]
                rows[f"{a}-{b}"] = dict(remainder_upper=str(up(r)), remainder_upper_decimal=float(arb(r).upper()), coefficients=len(cd))
        allr = [up(v[1]) for v in D[name].values()]
        asm[name] = dict(base_entries=rows, max_remainder_over_all_entries=float(max(allr)), entries=len(allr))
    S = pickle.load(open(schur_fn, "rb"))
    sch = {}
    for beta, v in sorted(S['res'].items(), key=lambda kv: Fr(kv[0])):
        sch[beta] = {}
        for q in ("detS", "trS"):
            cd, (rm, rr), m = v[q]
            r = arb(rm) + arb(rr) * arb("[+/- 1]")
            sch[beta][q] = dict(remainder_upper=str(up(r)), remainder_upper_decimal=float(r.upper()), order=m, coefficients=len(cd))
    json.dump(dict(description="Taylor-model remainder radii r (the model is |f - sum c_g p^g| <= r |p|_2^(m+1) on |p|_2 <= rho0)",
                   N=D['N'], m=D['m'], mLL=D['mLL'], rho0=str(D['rho0']), assembly=asm, schur=sch, schur_meta=S['meta']),
              open(f"{out}/remainders.json", "w"), indent=1)
    print(f"DUMPED corrections ({len(ent)}) and remainders ({len(sch)} shift levels) to {out}")

if __name__ == "__main__":
    main(*sys.argv[1:4])
