# Response to the referee report of 2 October 2026

**Manuscript:** *The square maximises the ratio of the first two Dirichlet eigenvalues among quadrilaterals* (Anonymous).
**Version reviewed:** draft of 2 October 2026, 12 pages, PDF SHA-256 `e537647f2f2245254010d1bd06113a7f02f5a242f0dd09b6f1cd762001092caf`
(the same file as the producer's `paper/PAPER_draft.pdf` of that date).
**Recommendation received:** major revisions. **Referee's audit bundle:** `quadrilateral_referee_audit_2026-10-02.zip`, whose checks were
read and whose findings are addressed below. The report's text and the bundle are not redistributed; each item is summarised here.

Every request is listed with its disposition and the location of the evidence. "Paper" refers to the revised manuscript
(`paper/PAPER.pdf`), "Supplement S" to `paper/SUPPLEMENT_S.pdf` (a typeset rendering of `darts/LEMMAS.md`, which governs).

## Major comments

### 4.1 Essential proof objects were missing

**Disposition: done.** The research archive now contains the complete supplementary proofs, all certificate inputs and final logs, the
verification sources with their manifests and dependency versions, the replay receipts and one entry point per component:

| component | entry point | inputs and outputs | receipt |
|---|---|---|---|
| Theorem 2 (local step) | `local/replay_local_v2.sh` | `local/replay_v2/` (exact correction coefficients, Taylor-model remainders, 960 direction boxes with exact sign bounds) | `local/replay_v2/RECEIPT.json` |
| Theorem 3 and the transport cover | `global2/replay.jl`, `global2/replay_t.jl` | `global2/runs/` | `global2/RECEIPT_transport.md`, `global2/transport_replay_clean_RECEIPT.json`, `global2/shell0225_replay_failclosed_full.log` |
| Lemma $B'$ | `bandB/run_replay.sh` | `bandB/runs/recert_chart1.jsonl`, `recert_chart2.jsonl` | `bandB/replay.log` |
| Theorem 5 (darts) | `darts/replay_all.sh` | `darts/certs4/`, `darts/certs/`, `darts/certs2/`, `darts/manifests/`, `darts/refcache/` | the fresh-extraction proof receipt in `darts/replays/receipts/` (see `REPLAY_RECEIPT.md`), verified by `verify/xreceipt.py` |

The archive distinguishes the final certificates from exploratory material: superseded configurations and pre-repair logs are excluded
(the README lists what was excluded and why). Checking the final certificates does not require regenerating the search.

**On a clean replay of the frozen archive.** Preparing the frozen archive exposed a portability defect. Nine pinned files of the dart
verifier hard-coded the producer's package path, so a copy could not have been replayed (the runtime binding check is fail-closed).
The paths were made location-relative, with no predicate changed (`darts/review/PORTABILITY_PATCH_2026-10-02.md`), and the verifier
was re-pinned and re-tested. The archive's component files (code, inputs, final logs, pins) were then frozen, extracted into a new
directory, and every component was replayed there through the local job queue: the local step (960 direction boxes, cube $1/20$; the exact correction coefficients, remainder radii and all 960 box records are reproduced byte for byte), Lemma $B'$ (44 268 certificates re-run, both charts tiled exactly), the shell (8 216 certificates, 8 217 leaves), the transport cover (31 426 certificates and all 18 756 exclusion and compound leaves) and the darts (all 462 849 leaves of the 14 logs, 2 shards per log); every replay passed with no failure. The release archive contains these
component files byte for byte. `replay_archive.sh` gives the same commands to any reader.

An independent rerun by an unaffiliated party has not taken place. The producer cannot supply that observation. The archive is
organised so that such a rerun can be made from a fresh extraction, and the assurance table of the paper records the dimension as
not established.

### 4.2 Non-convex proofs and the angular hypotheses of the thin-wedge lemma

**Disposition: done.**
- Supplement S contains the complete proofs, including the mixed-boundary bracketing (S 1.2, Lemmas 9.1–9.4, 9.6), the tail estimates
  (Lemmas 6.1, 11.1–11.3) and the exact chart accounting (§13).
- Lemma 5.1. The referee's example ($b_A = 30^\circ$, $b_C = 120^\circ$, $\rho_A=\rho_C=4$) is a valid dart that satisfies the tangent
  condition. The supplement's statement of Lemma 5.1 already required both angles to be acute and gave its own obtuse example
  ($b_A = 10^\circ$, $b_C = 160^\circ$); the paper's summary omitted the qualification. The rigorous predicate
  (`darts/regions.py`, `P_lemma51`) checks $a_x>0$ and $c_x>0$ on the whole box, and **no record of any final log uses Lemma 5.1**
  (0 occurrences in all 14 logs); the decomposition of §5.4 uses Lemma 5.3, whose hypothesis $\theta_B\le0.1772$ forces both angles to be
  acute. The separately written coverage checker also requires acute angles. Paper §5.2 now states the hypothesis and this fact.

### 4.3 The local neighbourhood

**Disposition: done.**
- Published (all in `local/replay_v2/`): the 1 180 correction coefficients $C_{jki}$ as exact dyadic rationals (`corrections.json`); the
  remainder radii of every assembled entry and of $\det S$, $\operatorname{tr}S$ for all 17 shift levels (`remainders.json`); the 960 direction
  boxes with, for every $s$-interval, the shift level, the criterion and the exact rational upper bound of $\det S/s^2$ or $\operatorname{tr}S/s^2$
  (`schur_cover2.json`).
- Coefficient deletion. An audit of the final chain found five deletion sites, all justified by exact identities; none relies on an
  enclosure containing zero. The step from the exact identities for the base matrices (Lemma 3.2, SymPy) to the corrected pencil was
  previously argued only in code comments; it is now Corollary 3.3 with a proof (Green's identity at $p=0$). The paper now separates the
  justification (the identity) from the consistency check (the enclosure contains $0$).
- The exposition now describes one configuration throughout ($N=10$, $m=12$, $m_{LL}=14$, $\rho_0=1/10$, cube $1/20$); the $1/25$
  configuration appears only as a separate second certificate. Two sentences quoting exploratory runs made before a code repair were
  removed.
- While preparing the published records we found that the choice of shift level compared the levels with a nearest-rounded floating-point
  value of a rigorous lower bound. The smallest margin in the certificate was $0.00117$, so no certificate was affected, but the comparison
  now uses the exact rational lower endpoint and the replay was re-run (2 October 2026, `local/replay_local_v2.sh` under the shared
  job queue): the code matched its manifest; the exact lemmas were re-verified; all 960 direction boxes were certified with 0 failures,
  cube half-width exactly $1/20$, the same box count as before; receipt `local/replay_v2/RECEIPT.json`, SHA-256
  `9ebbaa46d1d600154c10f34d6bd47b6f77fc181c8cae48de4edd3f6240feed10`, binding the hashes of the code and of every output.
- A separately written checker, `verify/xlocal.py`, now checks the published records in exact rational arithmetic: exact tiling of the
  eight faces by the 960 boxes, coverage of $s\in(0,1/20]$, admissibility of every shift level (circumscribed rectangle or Krahn–Szegő,
  re-derived from the vertex formulas) and the sign of every logged bound. It passes, and it rejects five corrupted copies of the record
  (`referee_response/xlocal.log`, `xlocal_negative_controls.log`). It does not recompute the Taylor-model enclosures.

### 4.4 A mathematical specification of the certificate primitives

**Disposition: done.** Appendix A specifies inputs, hypotheses and conclusions of 15 primitives, with the files that enforce them. In
particular:
- the two flux constructions are now distinguished in §4.2: the shell cover uses a global polynomial field with
  $\operatorname{div}\tau = -u_1-\Delta u_1/\mu_0$; the transport cover uses the Piola field with $\operatorname{div}\hat\tau = -\hat u_1$ and a single polynomial
  $\hat\psi$, and §4.2 gives the normal-continuity argument (both affine maps agree on the diagonal, $\hat\tau\cdot\hat n = 0$ there, $\hat\psi$ is continuous,
  and the Piola map preserves normal flux);
- the Sturm certificate (P8) states the cell rule that excludes interior zeros: the enclosure of $V_k-E$ must exclude $0$, $\omega_k\ell_k<\pi$
  must be proved on oscillatory cells, and the solution must be positive at every node. An audit confirmed that both implementations enforce
  this. It also found that the separately written checker computed the series remainder of a rarely used near-zero branch in floating point.
  That remainder is now computed in ball arithmetic (fail-closed unless $|z|<1$ is proved), and the separately written check of the band sample was re-run: exact tiling of both charts, all 8 949 exclusions re-decided, and 400 of 400 sampled
  certificates reproduced (382 with the default mesh, 18 after mesh refinement, the same as before the repair);
- the unbounded tails, the geometric containment conditions for inner and outer darts and the cut conditions are P11–P13.

### 4.5 Separating replay, independent implementation and sampling

**Disposition: done.** The paper's assurance table (§6) records, for each component, the written argument, exact coverage, exhaustive replay
and separately written checks.

Re-examining the records for this table showed that the draft's sentence on the transported-certificate sample was wrong. It said "52 of 60,
8 over budget, none fails"; the records gave 51 verified, 8 over the work budget, and 1 that had failed with an earlier version of the
separately written checker (before its exclusion fix) and had not been retried. All 9 unverified boxes were re-run with the current checker
and a five-fold work budget (200 units): all 9 verified, using 10 to 103 units. The separately written implementation therefore
re-certifies all 60 sampled boxes; the paper says which 9 needed the larger budget.

The dart certificates have no second implementation; this is stated in the paper. The publisher decided not to commission one for this
version.

## "Should fix" recommendations

- **Independent checking of the non-convex spectral bounds, tails and chart boundaries:** not done for this version (publisher's decision);
  recorded as an open assurance dimension. The specification (Appendix A, P10–P14) is given so that such a check can be made.
- **The eight transported boxes:** done; see 4.5.
- **A theorem-to-certificate dependency table and the trusted components:** done (Appendix B; the trusted base is listed in `ASSURANCE.md`).
- **The Arbon comparison by an explicit calculation:** done. For the trapezoid with parallel sides $1\pm t$ the four single-vertex
  first-order matrices are given exactly; they sum to zero, so the true first-order term vanishes, while the rule of (3.16) gives
  $-\sqrt{16384+729\pi^4}/(18\pi^2)\approx-1.664$ per unit $t$, or $-3/2$ when split into two simultaneous two-vertex moves; for a single
  two-vertex move the true slope is $-3/4$ and the rule gives $\approx-0.832$ (`referee_response/arbon_perturbation.py`; its vertex sums
  reproduce the matrices of Lemma 3.2 and the referee's).

## "Could improve"

- Diagrams: Figure 1 (the BD chart of a dart, drawn to scale) and Figure 2 (the transported test functions and Piola flux) are added to the paper.
- General certification propositions: Appendix A.
- A quantitative stability statement: not attempted.

## Minor comments

- **Bibliography:** completed. Every entry was checked against Crossref, MathSciNet, zbMATH, arXiv or the publisher, and every DOI
  resolves to a matching record (`referee_response/REFERENCES_VERIFIED.md`, which also lists what could only be checked through
  secondary sources). Arbon, Mannan, Psenka and Ragavan is cited as *J. Spectr. Theory* 12(2) (2022) 515–533, doi:10.4171/JST/409.
  The check made further corrections: the Payne–Pólya–Weinberger journal is *J. Math. and Phys.* 35 (1956) 289–298; the polygonal
  conjecture is attributed to Siudeja (2010, Conjecture 1.2, stated in the 2007 preprint) as well as to Antunes and Freitas (2008);
  the Krahn–Szegő inequality is cited to Krahn (1926), Hong (1954) and Pólya (1955); the constant $0.1893$ is cited to Liu (2015);
  Lehmann–Goerisch, Andrews–Clutterbuck, Hersch, Allegretto–Piepenbrink, Crouzeix–Raviart, Hadamard's formula (Henrot; Henrot–Pierre)
  and the software are now cited where used; Endo, Liu and Mariano (2026), previously listed but not cited, is cited in §4.2.
- **Parameter consistency:** one configuration throughout (see 4.3). $A_0,A_1,A_2$ are defined before Lemma 4.1. The two uses of $\rho_0$ are
  separated (the flux constant is now $\mu_0$), the $D_4$ rotation is $\sigma_r$, and Lemma B's boundary functions are $\overline y,\underline y$; a
  notation note lists the remaining deliberate reuses ($\beta$ in §§3–4 and §5).
- **Topology and terminology:** $\mathcal R$ is described as contained in the compact a-priori box. "Simple quadrilateral" is defined (interior
  angles different from $\pi$), and the degenerate case is a triangle, dispatched by the triangle theorem.
- **Claims:** the dart supremum $7/3$ is labelled a numerical observation, not proved and not used.
- **Presentation:** the review chronology is in Appendix C. Placeholders are removed, and the review rounds for the darts are counted as five.

## Changes to code and records made during the response

| file | change | re-run |
|---|---|---|
| `verify/xband.py` | rigorous ball-arithmetic remainder in the near-zero Sturm branch | `verify/xbandcover.py` sample, 400 boxes |
| `verify/xsample4.py` (new) | parallel driver with a configurable work budget for `verify/xcert.py` | the 9 unverified sampled boxes |
| `local/cover2.py` | exact rational shift-level bound; exact upper bounds recorded per box | `local/replay_local_v2.sh` |
| `local/dump_local_data.py` (new) | exact correction coefficients and remainder radii | the same replay |
| `local/replay_local_v2.sh` | manifest check first; hash-bound receipt at the end | the same replay |
| dart verifier (9 pinned files) and launcher | absolute package paths replaced by location-derived paths; interpreter from `PY`; re-pinned (`58bd96bf…`) | rigour lint, 82 infrastructure probes, fresh-extraction quick and proof replays |
| `verify/xreceipt.py`, `verify/test_xreceipt.py`, `verify/xdart.py` | location-derived paths; `--darts`, `--home`, `--recorded-darts`; mutations needing 2 shards report n/a on 1-shard receipts | quick receipt verified; mutation test |
| `verify/xlocal_negative_controls.py` (new) | five corruptions of the local-step record must be rejected | run |
| `darts/LEMMAS.md` | six wording corrections, a precise statement in §7, a citation year (`darts/review/LEMMAS_WORDING_2026-10-02.diff`) | Supplement S re-typeset and re-checked token by token |
| `verify/xlocal.py` (new) | separately written check of the published local-step records | run on the replay's records; 5 negative controls |

No certificate, predicate or log of the convex covers, the band or the darts was changed by the response. Two changes touch code that
a proof replay executes. The first is the exact shift-level comparison in the local cover; the cover was re-run, with the same 960
boxes. The second is the dart portability patch, after which the dart proof replay was run again.
