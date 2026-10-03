"""Six wording corrections (and one citation year) in darts/LEMMAS.md and identically in paper/SUPPLEMENT_S.md, found while
typesetting the supplement (2 October 2026).  No statement used by a certificate changes.  Exact single matches, else abort."""
import sys
R = "/Users/admin/claude-sandbox/quad-ppw-2026-09-30/"
E = [  # (LEMMAS old, LEMMAS new, SUPPLEMENT old, SUPPLEMENT new)
 ("Weyl gives λ₂(A_ex) ≥ −(e_A + ‖ΔA‖₂);", "Weyl gives λ₂(A_ex) > −(e_A + ‖ΔA‖₂) (strict, since λ₂(A_fl + ΔA) > 0);",
  "Weyl gives $\\lambda_2(A_{\\mathrm{ex}}) \\ge -(e_A + \\|\\Delta A\\|_2)$;",
  "Weyl gives $\\lambda_2(A_{\\mathrm{ex}}) > -(e_A + \\|\\Delta A\\|_2)$ (strict, since $\\lambda_2(A_{\\mathrm{fl}} + \\Delta A) > 0$);"),
 ("Arbon–Mannan–Psenka–Ragavan 2020", "Arbon–Mannan–Psenka–Ragavan 2022", "Arbon–Mannan–Psenka–Ragavan 2020", "Arbon–Mannan–Psenka–Ragavan 2022"),
 ("1.5 to Q_L (reflected so that the apex is at the origin):", "1.5 to Q_L (whose apex C is at the origin):",
  "1.5 to $Q_L$ (reflected so that the apex is at the origin):", "1.5 to $Q_L$ (whose apex $C$ is at the origin):"),
 ("If x* ≥ d_x apply 1.5 to Q_R with x* replaced by", "If x* ≥ d_x apply 1.5 to Q_R, reflected by x ↦ 1 − x so that its apex A is at the origin, with x* replaced by",
  "If $x^* \\ge d_x$ apply 1.5 to $Q_R$ with $x^*$ replaced by",
  "If $x^* \\ge d_x$ apply 1.5 to $Q_R$, reflected by $x \\mapsto 1 - x$ so that its apex $A$ is at the origin, with $x^*$ replaced by"),
 ("the part with b_X < 39.23° or > 140.77° is open.",
  "outside [40°, 140°] Lemma 9.5 does not apply (its bound cannot reach 5/2 at all when b_X ≤ 39.23° or b_X ≥ 140.77°), and there the face is covered by the BO2 families alone.",
  "the part with $b_X < 39.23^\\circ$ or $> 140.77^\\circ$ is open.",
  "outside $[40^\\circ, 140^\\circ]$ Lemma 9.5 does not apply (its bound cannot reach $5/2$ at all when $b_X \\le 39.23^\\circ$ or $b_X \\ge 140.77^\\circ$), and there the face is covered by the BO2 families alone."),
 ("BQ₁ and DQ₂ are parallel, of positive lengths L and L − k₀u)", "BQ₁ and DQ₂ are parallel, of positive lengths L√(1 + u²) and (L − k₀u)√(1 + u²))",
  "$BQ_1$ and $DQ_2$ are parallel, of positive lengths $L$ and $L - k_0u$)",
  "$BQ_1$ and $DQ_2$ are parallel, of positive lengths $L\\sqrt{1 + u^2}$ and $(L - k_0u)\\sqrt{1 + u^2}$)"),
 ("For the x-derivative: v² is smooth where positive and", "For the x-derivative: v² is smooth, v is smooth where positive, and",
  "For the $x$-derivative: $v^2$ is smooth where positive and", "For the $x$-derivative: $v^2$ is smooth, $v$ is smooth where positive, and"),
 ("the second eigenvalue is j_{ν,2} because", "the second eigenvalue is j²_{ν,2}/x*² because",
  "the second eigenvalue is $j_{\\nu,2}$ because", "the second eigenvalue is $j^2_{\\nu,2}/x^{*2}$ because"),
 ("μ > E is certified by", "μ ≥ E is certified by", "$\\mu > E$ is certified by", "$\\mu \\ge E$ is certified by"),
]
L = open(R + "darts/LEMMAS.md").read(); S = open(R + "paper/SUPPLEMENT_S.md").read()
for lo, ln, so, sn in E:
    if L.count(lo) != 1: sys.exit(f"ABORT: LEMMAS match count {L.count(lo)}: {lo[:60]!r}")
    if S.count(so) != 1: sys.exit(f"ABORT: SUPPLEMENT match count {S.count(so)}: {so[:60]!r}")
    L = L.replace(lo, ln); S = S.replace(so, sn)
open(R + "darts/LEMMAS.md", "w").write(L); open(R + "paper/SUPPLEMENT_S.md", "w").write(S)
print(f"applied {len(E)} edits to both files")
