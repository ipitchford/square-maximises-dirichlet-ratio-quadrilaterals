# Novelty gate — square maximises λ₂/λ₁ among (convex) quadrilaterals (Dirichlet Laplacian)

Search dates: 2026-09-30 and 2026-10-01. Corpus boundary: arXiv (abstract pages + export API queries), Semantic Scholar
citation graph (citations of arXiv:2207.05814, 2009.00927, 0707.3631), web search (2 rounds), full text of the cited sources.

## Frozen claim
Theorem 1: for every convex quadrilateral Q, λ₂(Q)/λ₁(Q) ≤ 5/2 with equality only for squares. (Darts: separate, see paper §5.)
Contribution units: (1) the statement (global, convex quadrilaterals); (2) a complete computer-assisted proof; (3) two analytic
reduction lemmas (near-triangle homothety lemma; thin-domain inscribed-rectangle lemma); (4) the local method (Schur complement of
transported square modes as a multivariate Taylor model); (5) an identification of two gaps in Arbon 2022's local argument.

## Inspected sources (load-bearing)
- Arbon 2022, arXiv:2207.05814 (full text): local maximality of the square (Thm 1.3), global case open; §3.2 says the parameter search
  "at this time" is not done. Gaps found in eq. (3.16)–(3.18) (paper §1).
- Arbon–Mannan–Psenka–Ragavan, J. Spectral Theory 2020 (arXiv:2009.00927) and Siudeja 2010 (arXiv:0707.3631): triangle case (cited).
- Endo–Osting, arXiv:2609.02784 (Sept 2026): NEUMANN analogue on convex quadrilaterals; same parametrisation and box-cover architecture.
  Not a collision (different eigenvalue problem; no λ₁ lower bound needed there).
- Endo–Liu–Mariano, arXiv:2605.04331 (May 2026): triangles, second-order shape derivatives, Lehmann–Goerisch lower bounds (method prior art).
- Citations of 2207.05814 (Semantic Scholar, 2026-09-30): only Bañuelos–Mariano 2024 (Pólya functional) and Endo–Liu 2022 — no follow-up on the quadrilateral ratio.
- arXiv listing queries (quadrilateral ∧ eigenvalue ∧ Dirichlet; "fundamental ratio" ∧ quadrilateral; Ashbaugh ∧ Benguria ∧ polygon;
  "ratio of the first two" ∧ Dirichlet ∧ polygon; quadrilateral ∧ PPW): no 2024-06…2026-10 paper on the global quadrilateral Dirichlet ratio.
- Web search hits checked and excluded: arXiv:2603.29015 (obstacles in a square), arXiv:2607.01135 (Ashbaugh–Benguria reciprocal-gap
  conjecture, all domains), arXiv:2609.28602 (Krejčiřík survey; no mention of the polygonal ratio problem).

## Outcome by contribution unit
| unit | outcome |
|---|---|
| statement (convex quadrilaterals) | CLEAR TO PROVE within the recorded boundary (not found as proved) |
| complete proof | CLEAR (bounded search) |
| Lemma A (near-triangle homothety) | BOUNDED UNCERTAINTY — elementary; likely folklore in spirit; present as a lemma with proof, no priority claim |
| Lemma B (thin domains, rectangle of length θD and height (1−θ)w) | BOUNDED UNCERTAINTY — elementary; no priority claim |
| local method | method prior art: Endo–Osting (second-order test at the square, Neumann), Kato–Temple/Lehmann–Goerisch classics; our multivariate-Taylor Schur complement is a variant, no priority claim beyond "we use" |
| Arbon gaps | our own reading of the published preprint; verified by exact first-order computation (sympy) and FEM |

Permitted language: "we are not aware of a proof of the global statement; the local statement is due to Arbon (2022)"; "not found in
the bounded search described here". Forbidden: "first proof", "open problem since 2008" (no source states it as open with that date).

---

## Retrospective addendum (2 October 2026): the statement for all simple quadrilaterals

The report above (30 September – 1 October 2026) froze the convex statement; the non-convex case (darts) was then proved separately
(paper §5, Supplement S). This addendum records the intake check for the published statement, Theorem 1 for all simple quadrilaterals. It
is written after the proof and is labelled retrospective; it does not change the dated record above.

**Searches (2 October 2026).** arXiv API (`export.arxiv.org`, most recent first; full log `referee_response/arxiv_search_2026-10-02.log`):
`quadrilateral AND Dirichlet AND eigenvalue` (6 results: Endo–Osting 2026 (Neumann), Bañuelos–Mariano 2024 (Pólya functional),
Dirac eigenvalues 2024, Arbon 2022, polygonal Faber–Krahn 2022, a finite-element paper 2017); `"non-convex" AND quadrilateral AND
eigenvalue` (0); `"fundamental ratio" AND quadrilateral` (1: Arbon 2022); `Payne AND Weinberger AND polygon` (0); `"ratio of the first two"
AND Dirichlet AND quadrilateral` (1: Arbon 2022); `quadrilateral AND "second eigenvalue" AND Laplacian` (0).

**External check.** The referee report of 2 October 2026 (supplied by the publisher; its text is not redistributed) included a bounded literature check (live web search,
primary sources on arXiv and EMS Press) and "did not locate a superseding global quadrilateral-ratio proof".

**Outcome.**
| unit | outcome |
|---|---|
| statement for all simple quadrilaterals (convex and non-convex) | CLEAR within the recorded boundary (not found as proved) |
| non-convex methods (Neumann bracketing of arms, one-dimensional bounds with tails, inclusion certificates) | method prior art in general form (domain monotonicity, Dirichlet–Neumann bracketing, Born–Oppenheimer-type reductions, Liu's verified finite-element bounds); no priority claim for the methods |

**Permitted language:** "we are not aware of a proof of the global statement for convex or non-convex quadrilaterals; the local statement
is due to Arbon (2022); not found in the bounded searches described here". **Forbidden:** "first proof", "resolves a long-standing
open problem", any priority date.

**Bibliographic corrections found by the citation audit of 2 October 2026** (`CITATION_AUDIT.md`). The dated record above is left as
written; these supersede it. Arbon, Mannan, Psenka and Ragavan appeared in *J. Spectr. Theory* 12(2) (2022) 515–533; 2020 is the
arXiv year. The polygonal ratio conjecture is stated by Siudeja (arXiv:0707.3631 v1, 24 July 2007, Conjecture 1.2, (1.13); published
2010) as well as by Antunes and Freitas (2008); attribution to Antunes and Freitas alone, used in early drafts, was incomplete.
