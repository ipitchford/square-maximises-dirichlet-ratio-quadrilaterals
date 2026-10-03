# Citation audit

Checked on 2 October 2026 against Crossref, MathSciNet (MRef), zbMATH Open (including JFM), arXiv, the publishers' pages and, for
software, GitHub, PyPI and Zenodo. Every DOI resolves to a record whose authors, title, volume and pages match. Full record:
[referee_response/REFERENCES_VERIFIED.md](referee_response/REFERENCES_VERIFIED.md).

## Load-bearing sources and the relationship checked

| Source | Used for | Checked |
|---|---|---|
| Endo and Osting (2026), arXiv:2609.02784 | parametrisation $Q_p$; Theorem 2.1 (every convex quadrilateral is similar to some $Q_p$) | statement of Theorem 2.1 and eqs. (2)–(3) read in the v1 PDF |
| Arbon, Mannan, Psenka and Ragavan (2022), *J. Spectr. Theory* 12(2) 515–533 | $\xi(T)\le7/3$ for every triangle (Lemma A; degenerate quadrilaterals) | bibliographic record; the theorem covers all triangles (completing Siudeja's acute case) |
| Siudeja (2010), *Indiana Univ. Math. J.* 59(3) 1097–1120 | acute triangle case; the polygonal conjecture (Conjecture 1.2, (1.13)) | Theorem 1.1 and Conjecture 1.2 read in arXiv:0707.3631 v1 (24 July 2007) and v2 |
| Arbon (2022), arXiv:2207.05814 | the earlier local result; the two gaps in its argument (§1) | Theorem 1.3, region $[0,0.0015]^2$, eqs. (3.16) and (3.18) read in the v1 text |
| Krahn (1926); Hong (1954); Pólya (1955) | $\lambda_2\ge2\pi j_{0,1}^2/\lvert\Omega\rvert$ (shift levels) | Krahn's λ₂ result confirmed in the JFM review; Hong and Pólya bibliographic records |
| Andrews and Clutterbuck (2011) | the gap $\lambda_2-\lambda_1\ge3\pi^2/D^2$ for convex domains | bibliographic record |
| Hersch (1960) | $\lambda_1\ge\pi^2/(4r^2)$ for convex domains | bibliographic record (Crossref) |
| Liu (2015) | the guaranteed Crouzeix–Raviart lower bound with $C_h=0.1893h$ | bibliographic record; the attribution of 0.1893 to this paper is confirmed by four independent citing works (primary text paywalled) |
| Lehmann (1949, 1950); Goerisch and Haunhorst (1985); Behnke and Goerisch (1994) | the Lehmann–Goerisch method (Lemma 4.1 is proved in the paper for $n=1$) | bibliographic records; a known DOI error in MathSciNet for Behnke–Goerisch avoided |
| Temple (1928); Kato (1949); Haynsworth (1968) | Temple's inequality; inertia additivity | bibliographic records (the results are also proved or stated with proof sketches in the paper) |
| Allegretto (1974); Piepenbrink (1974) | positive solutions bound the spectrum from below (one-dimensional tails) | bibliographic records |
| Johansson (2017); FLINT 3.6.0; python-flint 0.9.0; Arblib.jl 1.8.1 | ball arithmetic | records and versions recorded in the receipts |

## Context-only sources
Payne–Pólya–Weinberger (1956), Ashbaugh–Benguria (1991, 1992), Antunes–Freitas (2008), Laugesen–Siudeja (2017), Endo–Liu–Mariano
(2026), Crouzeix–Raviart (1973), Henrot (2006), Henrot–Pierre (2018), Meurer et al. (2017), Gustafsson–McBain (2020).

## Not fully verified (recorded, not load-bearing)
- That Antunes–Freitas (2008) states the polygonal conjecture: supported by two citing papers; the primary text was not accessible.
  The paper also cites Siudeja's earlier statement, which was read.
- The theorem numbers inside Liu (2015) and Henrot (2006) / Henrot–Pierre (2018); the paper cites these works without theorem
  numbers.
- That Laugesen–Siudeja (2017) number the polygonal conjecture 6.31; the paper cites the survey without a conjecture number.

## Corrections made during the audit
- Arbon, Mannan, Psenka and Ragavan: year 2022 (the draft gave the 2020 preprint year).
- Payne–Pólya–Weinberger: *J. Math. and Phys.* 35 (1956) 289–298. This is not the AIP *J. Math. Phys.*
- The polygonal conjecture is attributed to Siudeja as well as to Antunes and Freitas.
- The Krahn–Szegő inequality is cited to Krahn, Hong and Pólya, not to G. Szegő (1954), which is the Neumann result.
- The constant 0.1893 is cited to Liu (2015), not to Liu–Kikuchi (2010).
