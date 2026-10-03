# Verified bibliography — square maximises λ2/λ1 among quadrilaterals

Checked 2026-10-02. Sources consulted per item: Crossref REST API (`https://api.crossref.org/works/<DOI>`), MathSciNet MRef
(`https://mathscinet.ams.org/mathscinet-mref`), zbMATH Open API (`https://api.zbmath.org`, incl. JFM for pre-1940 items),
arXiv API and PDFs, publisher PDFs where open, GitHub/PyPI/Zenodo APIs for software. Every DOI below resolved in Crossref
(or Zenodo) to a record whose authors, title, volume and pages match; none was inferred.

Legend. **VERIFIED**: bibliographic fields match ≥1 publisher/indexing record (usually 2–3). **LIKELY**: a content claim supported
by ≥2 independent secondary sources but not checked in the primary text. **UNVERIFIED**: not confirmed; reason given.

## Corrections to the current draft reference list

| Draft entry | Problem | Correction |
|---|---|---|
| Arbon, Mannan, Psenka, Ragavan (2020) | Wrong year; no volume/pages | J. Spectr. Theory 12(2) (**2022**) 515–533; 2020 is the arXiv preprint 2009.00927. Also fix in-text "(J. Spectral Theory, 2020)" and "Arbon, Mannan, Psenka, Ragavan 2020". |
| Siudeja (2010) | Issue and DOI missing | 59(3), 1097–1120, doi:10.1512/iumj.2010.59.3744. Proves only the **acute** case (Thm 1.1, ineq. (1.9)). Pages 1097–1120 are right; zbMATH's "1087–1120" is a zbMATH error. |
| Antunes, Freitas (2008) | Incomplete; attribution not checkable | J. Phys. A: Math. Theor. 41(5), 055201 (19 pp.), doi:10.1088/1751-8113/41/5/055201. Siudeja's arXiv:0707.3631v1 (24 Jul 2007), published as Siudeja (2010), already states the polygonal ratio conjecture (Conj. 1.2, (1.13)). Item 6 has the recommended wording. |
| Ashbaugh, Benguria (1992) | No issue/pages/DOI; companion paper missing | Ann. of Math. (2) 135(3), 601–628, doi:10.2307/2946578; add Bull. AMS 25(1) (1991) 19–29. |
| Payne, Pólya, Weinberger (1956), "J. Math. Phys. 35" | Wrong abbreviation ("J. Math. Phys." is the AIP *Journal of Mathematical Physics*); no pages | **J. Math. and Phys.** 35 (1956) 289–298, doi:10.1002/sapm1956351289. |
| "Krahn (1926) / Szegő" | Incomplete. The relevant Szegő is not G. Szegő (1954) | Krahn (1926) + Hong (1954) + Pólya (1955). The "Szegő" is Peter Szegő, credited in Pólya (1955). G. Szegő, J. Rational Mech. Anal. 3 (1954) 343–356 is the Neumann (Szegő–Weinberger) result and should not be cited here. |
| "Temple (1928), Kato (1949) … Haynsworth (1968)" | No bibliographic data | Full entries in items 10–12. |
| Liu (2015) | No pages/DOI | Appl. Math. Comput. 267 (2015) 341–355, doi:10.1016/j.amc.2015.03.048. The constant 0.1893 comes from Liu (2015), not from Liu–Kikuchi (2010). |
| (missing) | Cited in text but absent from list | Lehmann–Goerisch, Andrews–Clutterbuck, Allegretto–Piepenbrink, Crouzeix–Raviart, Arb/FLINT, Hadamard-formula reference, software (items 13–20). |
| Endo, Liu, Mariano (2026) | Listed but not cited in the text | Cite it (suggested use under item 5) or drop it. |

## References (paper style)

### 1. Arbon (2022)
- Arbon, R. (2022). Global and local bounds on the fundamental ratio of triangles and quadrilaterals. Preprint, arXiv:2207.05814 [math.SP] (v1, 12 Jul 2022; only version).
  - Used for: prior local result for quadrilaterals (Theorem 1.3, region (t,s) ∈ [0,0.0015]²); critique of its eq. (3.16), which writes the first-order variation of λ2/λ1 under a two-vertex perturbation as the sum of the two one-vertex terms (confirmed in the v1 text, together with (3.18) and the [0,0.0015]² region).
  - Status: VERIFIED https://arxiv.org/abs/2207.05814 (arXiv API; PDF v1). No published version was found in Crossref, OpenAlex (only the arXiv DOI 10.48550/arXiv.2207.05814), Semantic Scholar or MathSciNet. zbMATH lists it as "Preprint, arXiv:2207.05814".

### 2. Arbon, Mannan, Psenka, Ragavan (2022)
- Arbon, R., Mannan, M., Psenka, M., Ragavan, S. (2022). A proof of the triangular Ashbaugh–Benguria–Payne–Pólya–Weinberger inequality. J. Spectr. Theory 12(2), 515–533. doi:10.4171/JST/409. Preprint: arXiv:2009.00927 (v1 2 Sep 2020; v2 31 Mar 2021).
  - Used for: ξ ≤ 7/3 for all triangles, with the equilateral triangle as maximiser. It completes Siudeja's acute case by treating right and obtuse triangles.
  - Status: VERIFIED https://api.crossref.org/works/10.4171/JST/409 (vol. 12, issue 2, pp. 515–533, published 2022-09-21); MathSciNet MR4487482; zbMATH Zbl 1498.35368; https://arxiv.org/abs/2009.00927 (arXiv record carries this DOI). The referee is right: the year is 2022.

### 3. Siudeja (2010)
- Siudeja, B. (2010). Isoperimetric inequalities for eigenvalues of triangles. Indiana Univ. Math. J. 59(3), 1097–1120. doi:10.1512/iumj.2010.59.3744. Preprint: arXiv:0707.3631.
  - Used for: "the triangle case was proved for acute triangles by Siudeja". The result is Theorem 1.1, ineq. (1.9): λ2/λ1(T) ≤ λ2/λ1(equilateral) for acute T. This is the correct Siudeja citation, and no other Siudeja paper is needed. The same paper also states the polygonal conjecture (Conjecture 1.2, ineq. (1.13)); see item 6.
  - Status: VERIFIED https://api.crossref.org/works/10.1512/iumj.2010.59.3744 and MathSciNet MR2779073 (both 1097–1120). The zbMATH review (Zbl 1220.35105, B. Dittmar) confirms the content: "for acute triangles also the equilateral triangle maximizes λ2/λ1". zbMATH prints 1087–1120, which is a zbMATH error. Theorem and conjecture numbers were checked in arXiv v1 (24 Jul 2007) and v2 (17 Jul 2008) only, because the IUMJ site did not resolve.
  - In-text caveat: wherever the draft invokes "the triangle theorem" for a possibly obtuse triangle (§1 definition of simple quadrilateral; thin-quadrilateral lemma), the operative citation is item 2. Siudeja covers acute triangles only.

### 4. Endo, Osting (2026)
- Endo, R., Osting, B. (2026). Maximizing the fundamental Laplace–Neumann eigenvalue on quadrilaterals. Preprint, arXiv:2609.02784 [math.AP] (v1, 2 Sep 2026).
  - Used for: parametrisation of quadrilaterals Q_p, with displacements δ_i(p) identical to the draft's, area |Q_p| = 1−a²−d² (their eq. (2)) and cross products c_i(p) (eq. (3)). Their Theorem 2.1 is a constructive proof that every convex quadrilateral is similar to Q_p for some p ∈ int P. Also used for the certified box-cover architecture (they certify with MATLAB/INTLAB and the solver veigs).
  - Status: VERIFIED https://arxiv.org/abs/2609.02784 (statement of Theorem 2.1, eqs. (2)–(3) and abstract checked in the v1 PDF). No published version found (Crossref title search; no arXiv journal-ref). Their reference square is [−1/2,1/2]² and the draft's is [0,1]²; the two differ only by a translation.

### 5. Endo, Liu, Mariano (2026)
- Endo, R., Liu, X., Mariano, P. (2026). Sharp Dirichlet eigenvalue inequalities on triangles. Preprint, arXiv:2605.04331 [math.SP] (v1, 5 May 2026; 31 pp.).
  - Used for: no use is stated in the draft, which lists this paper but never cites it. Suggested use: a precedent for certified two-sided finite-element eigenvalue bounds on polygons. Their Lemma 3.3 is Liu's CR lower bound with C_h = 0.1893h (Theorem 2.1 of Liu 2015). Their Lemma 3.4 is a finite-element Lehmann–Goerisch bound with an H(div) flux (ŵ, ∇v) = (û, v), the same form as the draft's Lemma 4.1. Cite at §4.2/§5.3, or remove the entry.
  - Status: VERIFIED https://arxiv.org/abs/2605.04331 (arXiv API; PDF v1). No published version found (Crossref).

### 6. Polygonal Ashbaugh–Benguria–PPW conjecture: attribution
- Antunes, P. R. S., Freitas, P. (2008). A numerical study of the spectral gap. J. Phys. A: Math. Theor. 41(5), 055201 (19 pp.). doi:10.1088/1751-8113/41/5/055201.
  - Used for: attribution of the conjecture that the regular k-gon maximises ξ among k-gons.
  - Status (bibliographic): VERIFIED https://api.crossref.org/works/10.1088/1751-8113/41/5/055201 (online 2008-01-23) and MathSciNet MR2433425.
  - Status (content): UNVERIFIED in the primary text. The IOP page sits behind a bot challenge and the paper is closed access. The abstract (via OpenAlex) confirms only that "a similar study is carried out for the quotient of the first two Dirichlet eigenvalues". The claim that it states the polygonal conjecture is LIKELY: Arbon, Mannan, Psenka and Ragavan (2022) cite "Conjecture 13 of [Antunes–Freitas]", and Arbon (2022) says it was "first stated by Antunes and Freitas".
- Siudeja (2010, item 3), Conjecture 1.2, ineq. (1.13), states λ2/λ1(P(n)) ≤ λ2/λ1(R(n)) for n-gons. This is VERIFIED in arXiv:0707.3631v1, dated 24 Jul 2007, six months before Antunes–Freitas appeared online.
- Laugesen, R. S., Siudeja, B. A. (2017). Triangles and other special domains. In: A. Henrot (ed.), Shape Optimization and Spectral Theory, De Gruyter Open, Warsaw, pp. 149–200. doi:10.1515/9783110550887-006 (open access, CC BY-NC-ND 3.0).
  - Used for (suggested addition): the survey statement of the conjecture (Conjecture 6.31).
  - Status: bibliographic data VERIFIED (Crossref; MathSciNet MR3681150). That Conjecture 6.31 is the polygonal ratio conjecture is LIKELY: Arbon–Mannan–Psenka–Ragavan (2022) say so, and the 2020 open-problem list (http://publish.illinois.edu/eigenvalues2020/files/2020/04/Open-problems.pdf, Open Problem 4) says "the general case was raised by Laugesen and Siudeja [Conjecture 6.31]". The chapter text could not be retrieved because the publisher serves a bot challenge.
- Recommended wording: "conjectured by Siudeja [2010, Conj. 1.2] and, on numerical evidence, by Antunes and Freitas [2008]; see also Laugesen–Siudeja [2017, Conj. 6.31]".

### 7. Ashbaugh, Benguria (1991, 1992)
- Ashbaugh, M. S., Benguria, R. D. (1992). A sharp bound for the ratio of the first two eigenvalues of Dirichlet Laplacians and extensions. Ann. of Math. (2) 135(3), 601–628. doi:10.2307/2946578.
- Ashbaugh, M. S., Benguria, R. D. (1991). Proof of the Payne–Pólya–Weinberger conjecture. Bull. Amer. Math. Soc. (N.S.) 25(1), 19–29. doi:10.1090/S0273-0979-1991-16016-7.
  - Used for: proof of the PPW conjecture (the disk maximises ξ; ξ ≤ j₁,₁²/j₀,₁² in the plane).
  - Status: VERIFIED https://api.crossref.org/works/10.2307/2946578 (JSTOR DOI; Crossref gives first page 601), MathSciNet MR1166646 (601–628); https://api.crossref.org/works/10.1090/S0273-0979-1991-16016-7, MathSciNet MR1085824.

### 8. Payne, Pólya, Weinberger (1956)
- Payne, L. E., Pólya, G., Weinberger, H. F. (1956). On the ratio of consecutive eigenvalues. J. Math. and Phys. 35(1–4), 289–298. doi:10.1002/sapm1956351289.
  - Used for: ξ ≤ 3 for planar domains and the PPW conjecture.
  - Status: VERIFIED https://api.crossref.org/works/10.1002/sapm1956351289; MathSciNet MR84696. The journal is the MIT Journal of Mathematics and Physics, now Studies in Applied Mathematics.

### 9. Krahn–Szegő (Hong–Krahn–Szegő) inequality λ2(Ω) ≥ 2πj₀,₁²/|Ω|
- Krahn, E. (1926). Über Minimaleigenschaften der Kugel in drei und mehr Dimensionen. Acta Comm. Univ. Tartu. (Dorpat.) A 9, 1–44. (No DOI.)
- Hong, I. (1954). On an inequality concerning the eigenvalue problem of membrane. Kōdai Math. Sem. Rep. 6(4), 113–114. doi:10.2996/kmj/1138843535.
- Pólya, G. (1955). On the characteristic frequencies of a symmetric membrane. Math. Z. 63, 331–337. doi:10.1007/BF01187944.
- Optional textbook statement: Henrot (2006) (item 19), Ch. 4 "The second eigenvalue of the Laplacian-Dirichlet", pp. 61–71 (chapter title and pages from Crossref doi:10.1007/3-7643-7706-2_4).
  - Used for: the lower bound λ2(Ω) ≥ 2πj₀,₁²/|Ω|; at fixed area, two equal disks minimise λ2.
  - Status: Krahn VERIFIED via JFM 52.0510.03 (zbMATH Open API). The review states the λ2 result explicitly: "Den kleinsten zweiten Eigenwert erhält man, wenn das ganze Gebiet aus zwei inhaltsgleichen Sphären besteht". Hong VERIFIED (MathSciNet MR70015; Crossref). Pólya VERIFIED bibliographically (Crossref; MathSciNet MR73047).
  - The attribution chain is LIKELY: Krahn in the 1920s, Hong in 1954, then Pólya (1955), whose final remark credits Peter Szegő. Source: Brasco, L., De Philippis, G. (2017), "Spectral inequalities in quantitative form", in Henrot (ed.), Shape Optimization and Spectral Theory, pp. 201–281, doi:10.1515/9783110550887-007, footnote 10 (arXiv:1604.05072); Pólya's text was not accessed.
  - Krahn (1926) works in dimension ≥ 3. The nodal-domain argument is dimension-free, but for the planar statement Hong (1954) or Pólya (1955) is the direct citation. The planar Faber–Krahn input is Krahn, E. (1925), Über eine von Rayleigh formulierte Minimaleigenschaft des Kreises, Math. Ann. 94, 97–100, doi:10.1007/BF01208645 (VERIFIED Crossref, JFM 51.0356.05).
  - Not to cite here: G. Szegő (1954), J. Rational Mech. Anal. 3, 343–356 (doi:10.1512/iumj.1954.3.53017). It is the Neumann/free-membrane inequality (Brasco–De Philippis §4).

### 10. Temple (1928)
- Temple, G. (1928). The theory of Rayleigh's principle as applied to continuous systems. Proc. R. Soc. Lond. A 119(782), 276–293. doi:10.1098/rspa.1928.0098.
  - Used for: Temple's inequality λ1 ≥ ρ − ε²/(β − ρ) (§3.1, primitive P2), and the Temple step in §5.3.
  - Status: VERIFIED https://api.crossref.org/works/10.1098/rspa.1928.0098.

### 11. Kato (1949)
- Kato, T. (1949). On the upper and lower bounds of eigenvalues. J. Phys. Soc. Japan 4(4–6), 334–339. doi:10.1143/JPSJ.4.334.
  - Used for: the Kato–Temple bound (§3.1).
  - Status: VERIFIED https://api.crossref.org/works/10.1143/JPSJ.4.334; MathSciNet MR38738.

### 12. Haynsworth (1968)
- Haynsworth, E. V. (1968). Determination of the inertia of a partitioned Hermitian matrix. Linear Algebra Appl. 1(1), 73–81. doi:10.1016/0024-3795(68)90050-5.
  - Used for: inertia additivity In(A) = In(A₁₁) + In(A/A₁₁) (§3, §4, primitive P1).
  - Status: VERIFIED https://api.crossref.org/works/10.1016/0024-3795(68)90050-5; MathSciNet MR223392.

### 13. Lehmann–Goerisch method
Core citations:
- Lehmann, N. J. (1949). Beiträge zur numerischen Lösung linearer Eigenwertprobleme. I. Z. Angew. Math. Mech. 29(11–12), 341–356. doi:10.1002/zamm.19502911005.
- Lehmann, N. J. (1950). Beiträge zur numerischen Lösung linearer Eigenwertprobleme. II. Z. Angew. Math. Mech. 30(1–2), 1–16. doi:10.1002/zamm.19500300101.
- Goerisch, F., Haunhorst, H. (1985). Eigenwertschranken für Eigenwertaufgaben mit partiellen Differentialgleichungen. Z. Angew. Math. Mech. 65(3), 129–135. doi:10.1002/zamm.19850650302.
- Behnke, H., Goerisch, F. (1994). Inclusions for eigenvalues of selfadjoint problems. In: J. Herzberger (ed.), Topics in Validated Computations (Oldenburg, 1993), Stud. Comput. Math. 5, North-Holland/Elsevier, Amsterdam, pp. 277–322. ISBN 0-444-81685-2. (No DOI.)

Optional:
- Lehmann, N. J. (1963). Optimale Eigenwerteinschließungen. Numer. Math. 5, 246–272. doi:10.1007/BF01385896.
- Nakao, M. T., Plum, M., Watanabe, Y. (2019). Numerical Verification Methods and Computer-Assisted Proofs for Partial Differential Equations. Springer Ser. Comput. Math. 53, Springer, Singapore. Ch. 10 "Eigenvalue bounds for self-adjoint eigenvalue problems", pp. 349–411. doi:10.1007/978-981-13-7669-6_10.
- Liu, X. (2024). Guaranteed Computational Methods for Self-Adjoint Differential Eigenvalue Problems. SpringerBriefs in Mathematics, Springer, Singapore. Ch. 5 "Lehmann–Goerisch method for high-precision eigenvalue bounds", pp. 77–92. doi:10.1007/978-981-97-3577-8_5.

Notes:
  - Used for: Lemma 4.1 (Lehmann–Goerisch with n = 1 and a flux w ∈ L²(Ω)², (∇f, w) = (f, u)). The Goerisch extension is Goerisch–Haunhorst (1985). Per the zbMATH review, it relaxes the boundary conditions on the auxiliary function w in the Lehmann–Maehly procedure, which makes the method practical for PDEs. Behnke–Goerisch (1994) is the standard survey. The two optional books state the flux (H(div)) form used in the draft.
  - Status: all VERIFIED via Crossref and MathSciNet (MR34511, MR34512, MR789949, MR1318957, MR155425, MR3971222, MR4807655) and zbMATH (Zbl 0034.37901, 0568.73051, 0838.65060). Record defects to avoid copying:
    - Crossref misspells Lehmann I's author as "Lekmann" and drops "II" from Lehmann II's title; use the forms above.
    - The DOI 10.1002/zamm.19502911005 contains "1950" but is the 1949 vol. 29 paper.
    - MathSciNet gives Behnke–Goerisch the DOI 10.1016/0021-8502(94)90369-7. That DOI is wrong: Crossref resolves it to an unrelated *Journal of Aerosol Science* abstract. Cite the paper without a DOI.

### 14. Liu (2015); source of the CR constant 0.1893
- Liu, X. (2015). A framework of verified eigenvalue bounds for self-adjoint differential operators. Appl. Math. Comput. 267, 341–355. doi:10.1016/j.amc.2015.03.048.
  - Used for: the guaranteed lower bound λ_k ≥ λ_{k,h}/(1 + C_h²λ_{k,h}) with C_h = 0.1893h for Crouzeix–Raviart elements (§5.3).
  - Status: bibliographic data VERIFIED (https://api.crossref.org/works/10.1016/j.amc.2015.03.048; MathSciNet MR3399052). That the constant 0.1893 originates here is LIKELY: four independent works attribute it to Liu (2015), including the author's own later papers.
    - Endo–Liu (arXiv:2209.13415): "In [Liu 2015], it is shown through rigorous numerical computation that C(T) ≤ 0.1893h_T for triangles of arbitrary shapes".
    - Endo–Liu–Mariano (2026), Lemma 3.3: "Theorem 2.1 of [Liu 2015]".
    - Gómez-Serrano–Orriols (arXiv:1911.06758): "[Liu 2015, Theorem 2.1 and Remark 2.2]".
    - Vejchodský (arXiv:1602.08961): "optimal bound κ ≤ 0.1893 from [Liu 2015]".
  - The primary text was not accessed (Elsevier paywall), so the Theorem 2.1 / Remark 2.2 numbering is LIKELY rather than VERIFIED.
- Liu, X., Kikuchi, F. (2010). Analysis and estimation of error constants for P0 and P1 interpolations over triangular finite elements. J. Math. Sci. Univ. Tokyo 17(1), 27–78. (No DOI.)
  - Status: VERIFIED (MathSciNet MR2676659; zbMATH Zbl 1248.65118). This is **not** the source of 0.1893. Per the zbMATH summary, it treats the P0 and P1 interpolation constants. Cite it only if interpolation-constant analysis is discussed.
- Optional, for the sharpness of 0.1893: Endo, R., Liu, X. (2023). Shape optimization for the Laplacian eigenvalue over triangles and its application to interpolation error analysis. J. Differential Equations 376, 750–772. doi:10.1016/j.jde.2023.09.016. It proves that the equilateral triangle maximises the CR constant C(T). VERIFIED Crossref; content checked in arXiv:2209.13415.
- Optional, for priority: Carstensen, C., Gedicke, J. (2014). Guaranteed lower bounds for eigenvalues. Math. Comp. 83(290), 2605–2629. doi:10.1090/S0025-5718-2014-02833-0. This earlier CR lower bound uses κ² = 1/8 + j₁,₁⁻² (κ ≈ 0.4394), not 0.1893. VERIFIED (AMS PDF; MathSciNet MR3246802).

### 15. Andrews, Clutterbuck (2011)
- Andrews, B., Clutterbuck, J. (2011). Proof of the fundamental gap conjecture. J. Amer. Math. Soc. 24(3), 899–916. doi:10.1090/S0894-0347-2011-00699-1.
  - Used for: λ2 − λ1 ≥ 3π²/D² for convex domains, used in the choice of shift β (§4.2).
  - Status: VERIFIED https://api.crossref.org/works/10.1090/S0894-0347-2011-00699-1; MathSciNet MR2784332.

### 16. Allegretto (1974); Piepenbrink (1974)
- Allegretto, W. (1974). On the equivalence of two types of oscillation for elliptic operators. Pacific J. Math. 55(2), 319–328. doi:10.2140/pjm.1974.55.319.
- Piepenbrink, J. (1974). Nonoscillatory elliptic equations. J. Differential Equations 15(3), 541–550. doi:10.1016/0022-0396(74)90072-2.
  - Used for: the Allegretto–Piepenbrink argument (a positive solution gives a spectral lower bound) in Sturm shooting on unbounded intervals (§5.2).
  - Status: VERIFIED via Crossref (both DOIs) and MathSciNet MR374628, MR342829.

### 17. Johansson (2017); FLINT
- Johansson, F. (2017). Arb: efficient arbitrary-precision midpoint-radius interval arithmetic. IEEE Trans. Comput. 66(8), 1281–1292. doi:10.1109/TC.2017.2690633.
- The FLINT team (2026). FLINT: Fast Library for Number Theory. Version 3.6.0. https://flintlib.org.
  - Used for: ball arithmetic underlying all certificates. Arb is now part of FLINT 3.
  - Status: Johansson VERIFIED https://api.crossref.org/works/10.1109/TC.2017.2690633 and MathSciNet MR3681746. The FLINT citation form is VERIFIED at https://flintlib.org/citation.html, which also recommends citing Johansson (2017) for ball arithmetic. FLINT v3.6.0 was released 2026-06-29 (GitHub API).

### 18. Crouzeix, Raviart (1973)
- Crouzeix, M., Raviart, P.-A. (1973). Conforming and nonconforming finite element methods for solving the stationary Stokes equations I. Rev. Française Automat. Informat. Recherche Opérationnelle Sér. Rouge 7(R-3), 33–75. doi:10.1051/m2an/197307R300331.
  - Used for: the Crouzeix–Raviart element (§5.3).
  - Status: VERIFIED https://api.crossref.org/works/10.1051/m2an/197307R300331; MathSciNet MR343661.

### 19. Hadamard derivative of (possibly multiple) Dirichlet eigenvalues
- Henrot, A. (2006). Extremum Problems for Eigenvalues of Elliptic Operators. Frontiers in Mathematics, Birkhäuser, Basel. doi:10.1007/3-7643-7706-2.
- Henrot, A., Pierre, M. (2018). Shape Variation and Optimization: A Geometrical Analysis. EMS Tracts in Mathematics 28, European Mathematical Society, Zürich. doi:10.4171/178.
- Classical primary source (optional): Rousselet, B. (1983). Shape design sensitivity of a membrane. J. Optim. Theory Appl. 40(4), 595–623. doi:10.1007/BF00933973.
  - Used for: Hadamard's first-order formula. For a multiple eigenvalue, the one-sided derivatives are the eigenvalues of the matrix (−∫_∂Ω ∂ₙu_i ∂ₙu_j V·n). This is used in the 2×2 matrices M_k of §1 and Lemma 3.2, and in the critique of Arbon's (3.16).
  - Status: bibliographic data VERIFIED (Crossref; MathSciNet MR2251558, MR3791463, MR717179). Content pointers:
    - Henrot (2006), Ch. 2 "Tools" (pp. 17–44, Crossref doi:10.1007/3-7643-7706-2_2), covers derivatives of eigenvalues with respect to the domain (zbMATH review, Zbl 1109.35081).
    - Henrot–Pierre: Ch. 5 is on shape derivatives (zbMATH review of the 2005 French edition, Zbl 1098.49001). The 2018 English edition updates that chapter (Zbl 1392.49001); its chapter numbering was not checked.
  - UNVERIFIED: the exact theorem numbers for the multiple-eigenvalue statement, and whether Rousselet treats multiple eigenvalues (no text access). Check these before citing a theorem number, including that the stated hypotheses cover convex polygons.

### 20. Software
- python-flint: Johansson, F., Benjamin, O., et al. python-flint: Python bindings for FLINT, version 0.9.0 (2026). https://github.com/flintlib/python-flint (MIT licence).
  - Status: VERIFIED by the GitHub API (repository exists) and PyPI JSON (0.9.0 uploaded 2026-07-03; authors listed as Fredrik Johansson and Oscar Benjamin). The local `analysis-py313` environment has python-flint 0.9.0, bundling FLINT 3.6.0; `bandB/BAND_SPEC.md` records "python-flint 0.9". The project has no DOI or CITATION file; cite it together with FLINT and Johansson (2017).
- Arblib.jl: Dahne, J., Kaluba, M., Timme, S., et al. (2026). Arblib.jl, version 1.8.1. Zenodo. doi:10.5281/zenodo.22073429 (concept DOI for all versions: 10.5281/zenodo.10081694). https://github.com/kalmarek/Arblib.jl (MIT licence).
  - Status: VERIFIED via the Zenodo API (v1.8.1 record dated 2026-08-24; creators as listed, plus further contributors) and the GitHub API. The version used, Arblib 1.8.1 with Julia 1.12.7, is recorded in `global2/README_global2.md` and the replay receipts.
  - Optional Julia citation: Bezanson, J., Edelman, A., Karpinski, S., Shah, V. B. (2017). Julia: a fresh approach to numerical computing. SIAM Rev. 59(1), 65–98. doi:10.1137/141000671 (VERIFIED Crossref).
- scikit-fem: Gustafsson, T., McBain, G. D. (2020). scikit-fem: A Python package for finite element assembly. J. Open Source Softw. 5(52), 2369. doi:10.21105/joss.02369. https://github.com/kinnala/scikit-fem (BSD-3-Clause).
  - Status: VERIFIED. This is the preferred citation in the repository's CITATION.cff, and the DOI was checked at https://api.crossref.org/works/10.21105/joss.02369. Version 12.0.2 is installed in `analysis-py313` (PyPI upload 2026-06-05). No scikit-fem version is recorded in the project logs, so the version used in the final runs should be confirmed and stated.

## Items not fully verified (summary)
1. Antunes–Freitas (2008) stating the polygonal conjecture: LIKELY. The primary text is inaccessible and Siudeja (2007/2010) predates it.
2. Laugesen–Siudeja (2017), Conjecture 6.31: LIKELY. Chapter text not retrieved.
3. Pólya (1955), final remark crediting Peter Szegő: LIKELY. Secondary source only.
4. Liu (2015), Theorem 2.1 / Remark 2.2 numbering for 0.1893: LIKELY. Paywalled; four citing sources agree.
5. Siudeja (2010), theorem and conjecture numbering in the published version: checked on arXiv v1/v2 only.
6. Theorem numbers for the multiple-eigenvalue Hadamard formula in Henrot (2006), Henrot–Pierre (2018) and Rousselet (1983): not checked.
