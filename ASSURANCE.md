# Assurance boundary

| Dimension | State | Evidence and limits |
|---|---|---|
| Availability | checked at publication | The GitHub release and the Zenodo record, compared byte for byte with the release SHA-256 sums. The result is recorded on the Evidence Press page and in its machine records, not in this frozen file |
| Internal replay | passed | Exhaustive fail-closed replays of every certified cover, with hash-bound receipts. The archive as frozen was replayed from a fresh extraction ([REPLAY_RECEIPT.md](REPLAY_RECEIPT.md)) |
| Separately written checks (same producer) | partial | Coverage and every exclusion of the shell and transport covers (`verify/xcover.py`) and of the band (`verify/xbandcover.py`). Certificates of random samples: 60 of 31 426 transport certificates (`verify/xcert.py`); 400 band certificates (`verify/xband.py`). The local-step records (`verify/xlocal.py`). Dart coverage by 160 000 sampled darts (`verify/xdart.py`). The receipt state machine (`verify/xreceipt.py`). The dart certificates and the Taylor-model enclosures of the local step have no second implementation |
| Environment reproducibility | partial | Versions in [ENVIRONMENT.txt](ENVIRONMENT.txt), [requirements.txt](requirements.txt) and `environment/julia/` (Project and Manifest). The dart replay binds a runtime inventory of the Python installation. The replays run for hours; CI runs only the structural checks |
| Independent rerun | not established | No unaffiliated party has replayed the archive |
| Independent reimplementation | not established | The separately written checks above were written by the same producer workflow |
| Formal verification | not established | No proof assistant is used |
| Specialist review | not established | Reviews were adversarial model reviews, plus an external referee report of unrecorded authorship supplied by the publisher |
| Editorial peer review | not established | No journal or venue has reviewed the work |

**Load-bearing trust (trusted base).**
1. The written arguments: Lemmas A, B, 3.1, 3.2, Corollary 3.3, Lemma 4.1 and §4.2, Lemma $B'$ and its Sturm rule, and Supplement S.
2. The correctness of the verifiers:
   - `local/` (Taylor-model assembly and Schur complement);
   - `global2/` (Julia: assembly, Taylor models, criterion, transport, replay);
   - `bandB/` (Julia: one-dimensional certificates);
   - `darts/` (Python: finite-element, one-dimensional and inclusion certificates, replay drivers).
3. Arb (FLINT 3.6.0) through python-flint 0.9.0 and Arblib.jl 1.8.1; Julia 1.12.7; CPython 3.13.5.
4. SymPy 1.14.0 for the exact lemmas of §3. scikit-fem 12.0.2 supplies mesh connectivity only; node coordinates are exact dyadic
   rationals.
5. Floating point (NumPy, SciPy, BLAS) only chooses test functions, fluxes and splits. Every inequality that enters the proof is
   decided in ball or exact arithmetic.

**Known historical defects, all repaired before the replays reported here (paper Appendix C).**
- Shell and transport replays: early versions did not reject unknown statuses. They were replaced by fail-closed versions and re-run.
- Lemma $B'$: an early certificate allowed a zero inside a Sturm cell. The band was re-certified with the cell rule and replayed.
- Dart pipeline: certificate and driver defects found in five adversarial review rounds; the families were re-certified and replayed.
- Separately written band checker: a floating-point series remainder in a rare branch was replaced by ball arithmetic, and the
  sample was re-run.
- Local cover: the shift-level choice compared levels with a nearest-rounded value. It now uses the exact rational lower endpoint
  (smallest margin 0.00117 before the fix, so no certificate was affected), and the cover was re-run.
- Dart verifier: absolute paths replaced by location-derived ones (portability patch), followed by a fresh proof replay.
