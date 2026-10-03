# The square maximises the ratio of the first two Dirichlet eigenvalues among quadrilaterals

This is an **unrefereed candidate**, version `0.1.0-candidate`, dated 2026-10-03. DOI: [10.5281/zenodo.23109314](https://doi.org/10.5281/zenodo.23109314).
The creator is Anonymous and the publisher is Evidence Press.

- Paper: [paper/PAPER.pdf](paper/PAPER.pdf) (authoritative). Source: [paper/PAPER.md](paper/PAPER.md).
- Supplement S, on non-convex quadrilaterals: [paper/SUPPLEMENT_S.pdf](paper/SUPPLEMENT_S.pdf). It is typeset from [darts/LEMMAS.md](darts/LEMMAS.md), which governs.
- Agent-readable index: [AI_INDEX.md](AI_INDEX.md). Status: [STATUS.md](STATUS.md). Assurance boundary: [ASSURANCE.md](ASSURANCE.md).
- Replay: [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md) and [replay_archive.sh](replay_archive.sh).
- Response to the external referee report of 2 October 2026: [REFEREE_RESPONSE.md](REFEREE_RESPONSE.md).

## Result

**Theorem 1.** Let $Q\subset\mathbb R^2$ be a simple quadrilateral, convex or not, and let $\lambda_1(Q)<\lambda_2(Q)$ be the first two
eigenvalues of the Dirichlet Laplacian on $Q$. Then $\lambda_2(Q)/\lambda_1(Q)\le 5/2$, with equality if and only if $Q$ is a square.

- This is the case $k=4$ of the polygonal Payne–Pólya–Weinberger conjecture. Siudeja (2010, Conjecture 1.2) and, on numerical
  evidence, Antunes and Freitas (2008) conjectured that the regular $k$-gon maximises the ratio among $k$-gons.
- The case $k=3$ is known: Siudeja (2010) proved it for acute triangles, and Arbon, Mannan, Psenka and Ragavan (2022) for all triangles.
- For quadrilaterals, Arbon (2022, arXiv:2207.05814) proved that the square is a strict local maximiser in a small explicit
  neighbourhood. Paper §1 identifies two gaps in that argument; the local statement follows from Theorem 2 here.
- The bounded search found no earlier global result for quadrilaterals ([PRIOR_ART.md](PRIOR_ART.md)). That is a statement about
  the search, not a priority claim.

## How the proof is organised

The proof is computer-assisted. The convex case uses the parametrisation of Endo and Osting.

| Part | Region | Method | Certified objects | Replay entry point |
|---|---|---|---|---|
| Lemmas A, B (§2) | near-triangles; width/diameter $\le 1/10$ | analytic | none | none (constants in exact arithmetic) |
| Theorem 2 (§3) | $0<\lvert p\rvert_\infty\le 1/20$ | inertia criterion for a $3\times3$ Rayleigh–Ritz pencil; Kato–Temple bound as Taylor models | 960 direction boxes | `local/replay_local_v2.sh` |
| Theorem 3 (§4) | shell $0.049<\lvert p\rvert_\infty\le0.225$ (with Theorem 2: $0<\lvert p\rvert_\infty\le0.225$) | Rayleigh–Ritz; Lehmann–Goerisch with an exact flux | 8 216 certificates (8 217 leaves) | `global2/replay.jl` |
| Transport cover (§4) | rest of the convex fundamental domain | test functions and flux transported by a piecewise-affine map and a Piola transform | 31 426 certificates; 18 756 exclusion and compound leaves | `global2/replay_t.jl` |
| Lemma $B'$ (§4.5) | $1/10<$ width/diameter $\le 2/5$ | one-dimensional lower bound certified by Sturm oscillation; weighted one-dimensional Ritz bound | 44 268 certificates; 8 949 exclusions | `bandB/replay.jl` |
| Theorem 5 (§5, Supplement S) | non-convex quadrilaterals (darts) | analytic lemmas; certified covers (finite elements, one-dimensional bounds, inclusion) | 462 849 leaves in 14 logs | `darts/replay_all.sh` |

Floating point is used only to choose test functions, fluxes, correction coefficients and box splits. Every inequality that enters
the proof is decided in ball arithmetic (Arb) or exact rational arithmetic.

## Status and limitations

- Nothing here has been formally verified. No unaffiliated party has reproduced it, and no specialist or journal has reviewed it.
- An external referee report of 2 October 2026, supplied by the publisher, recommended major revisions. Every item was addressed
  ([REFEREE_RESPONSE.md](REFEREE_RESPONSE.md)). The report's text is not redistributed.
- During development, a second AI model (OpenAI) reviewed the parts adversarially (paper, Appendix C).
- The dart certificates have no second, independently written implementation. The publisher decided not to commission one for this
  version. The coverage of the dart decomposition is checked by an exact tiling and accounting in the replay, and by sampling.
- The trusted base is listed in [ASSURANCE.md](ASSURANCE.md).

## Replay

```bash
./replay_archive.sh check        # minutes: manifests, published-record checks, receipt verification, negative controls
./replay_archive.sh local        # about 10 min, 4 GB: Theorem 2 from the exact assembly
./replay_archive.sh band         # about 10 min on 3 threads: Lemma B'
./replay_archive.sh shell        # about 50 min on 3 threads: Theorem 3
./replay_archive.sh transport    # about 2 h on 3 threads: the transport cover
./replay_archive.sh darts        # hours (about 5.5 h with SHARDS=2, 1.7 GB per shard): Theorem 5
```

Requirements:
- Python 3.13 with python-flint 0.9.0, NumPy, SciPy, scikit-fem 12.0.2, mpmath and SymPy ([requirements.txt](requirements.txt)).
- Julia 1.12.7 with Arblib.jl 1.8.1 ([environment/julia/](environment/julia/)).
- Set `PY` and `JULIA` to the interpreters.

Recorded versions are in [ENVIRONMENT.txt](ENVIRONMENT.txt). The results of the fresh-extraction replay of this archive are in
[REPLAY_RECEIPT.md](REPLAY_RECEIPT.md).

## Contents

| Path | Content |
|---|---|
| `paper/` | the paper (Markdown source, PDF), Supplement S, figures and their generator, the reviewed draft (`history/`) |
| `local/` | Theorem 2: exact assembly, Schur-complement Taylor models, the direction cover, the exact lemmas, the replay and its records (`replay_v2/`: exact correction coefficients, remainder radii, the 960 boxes with exact bounds, receipt) |
| `global2/` | Theorem 3 and the transport cover: Julia verifier, final logs (`runs/`), replays and receipts, specifications |
| `bandB/` | Lemma $B'$: certificates, the two chart logs, replay and fail-closed tests |
| `darts/` | Theorem 5: `LEMMAS.md` (Supplement S), verifier, manifests, final logs (`certs/`, `certs2/`, `certs4/`), pinned inputs (`certs3/` witnesses, `refcache/`), pin (`review/PINNED.sha256`), receipts and run records |
| `verify/` | separately written checkers: coverage and exclusions (`xcover.py`, `xbandcover.py`), sampled certificates (`xcert.py`, `xband.py`), local-step records (`xlocal.py`), dart coverage sampling (`xdart.py`), receipt state machine (`xreceipt.py`) and its mutation test |
| `referee_response/` | records made while answering the referee: the exact Arbon calculation, re-runs of the separately written checks, the verified bibliography |
| `CLAIMS.json`, `AI_INDEX.md` | machine-readable claims, dependencies and evidence pointers |
| `PRIOR_ART.md`, `CITATION_AUDIT.md`, `SOURCES.md` | bounded prior-art search, citation checks, sources |
| `PROVENANCE.md`, `ENVIRONMENT.txt`, `LICENSES.md`, `MANIFEST.sha256` | provenance, environment, licences, hashes of every file |
| `fresh_replay/` | logs and records of the fresh-extraction replay (`REPLAY_RECEIPT.md`) |
| `LARGE_FILES.md`, `tools/check_git_manifest.py` | the bulk files kept only in the release archive (not in the Git repository), and a checker for a Git checkout or an extraction |

Excluded from the archive: superseded configurations, exploratory and pilot runs, pre-repair logs not needed by any replay, and
regenerable binary intermediates. The list, with reasons, is in [STATUS.md](STATUS.md).

## Load-bearing points for a reviewer

1. Lemma 3.1 (inertia criterion), Lemma 3.2 and Corollary 3.3 (exact first-order structure and the coefficient deletions).
2. Lemma 4.1 and the two flux constructions of §4.2, including the normal continuity of the transported Piola field.
3. The Sturm cell rule (Appendix A, P8) used by Lemma $B'$ and the one-dimensional dart bounds.
4. Supplement S: mixed-boundary bracketing (Lemmas 9.1–9.4, 9.6), tails (Lemmas 6.1, 11.1–11.3), the exact chart accounting (§13).
5. The fail-closed replay semantics and the receipt state machine (`verify/xreceipt.py`).

## Licence

Prose, claims, data and diagrams: CC0-1.0 ([LICENSE](LICENSE)). Code: MIT ([LICENSE-CODE](LICENSE-CODE)). Component map: [LICENSES.md](LICENSES.md).
