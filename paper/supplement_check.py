#!/usr/bin/env python3
"""Mechanical fidelity check: darts/LEMMAS.md (reviewed source) vs paper/SUPPLEMENT_S.md (typeset rendering).

Standard library only.  Usage (from any directory):

    python3 supplement_check.py            # prints the report and writes supplement_check.log next to this file

Checks
  (0) front matter: title/author/date and the single italic opening paragraph are exactly as specified;
      the comparison body of SUPPLEMENT_S.md starts at its first '#' header (front matter and opening
      paragraph are excluded from (a)-(e)).
  (a) the multiset of numeric literals (integers and decimals) is IDENTICAL.  Unicode sub/superscript digits
      in the source are normalised to ASCII digits, each script run forming its own literal (so N₁₂² gives
      12 and 2, as N_{12}^2 does); LaTeX command names (\\lambda, \\circ, ...) are skipped, so 174.29^\\circ
      counts as 174.29; signs are not part of a literal.  The ORDERED sequence of literals is compared too.
  (b) the ordered list of section headers, and the ordered list of labels Lemma/Lemmas x.y, Theorem x.y|X,
      Proposition/Prop x.y, Corollary, Remark(s), are identical.
  (c) every backticked code name of LEMMAS.md appears, backticked, in SUPPLEMENT_S.md (multisets compared).
  (d) token skeleton: both files are reduced to a sequence of tokens (letters, numeric literals, named
      mathematical symbols, and the punctuation ( ) [ ] | ‖ < > = + - / * ' " , ; : . § – —), with Unicode
      symbols mapped to their LaTeX names (≤ -> \\le, λ -> \\lambda, ∎ -> \\square, ...) and pure LaTeX
      layout/font markup dropped (\\mathrm, \\text, \\operatorname, \\bigl, \\qquad, \\newline, \\, braces, $, &, \\\\).
      The sequences are compared section by section; every difference is printed.  A small number of
      notation canonicalisations are applied and listed in the log (matrix notation, the range 0..30).
  (e) pdflatex safety: no non-ASCII character inside math; outside math and code only — – § é … occur.
Exit status 0 iff (0)-(e) all pass.
"""

from __future__ import annotations

import collections
import difflib
import hashlib
import re
import sys
import unicodedata
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "darts" / "LEMMAS.md"
SUP = HERE / "SUPPLEMENT_S.md"
LOG = HERE / "supplement_check.log"

EXPECTED_YAML = {
    "title": "Supplement S. Non-convex quadrilaterals (darts): lemmas and proofs",
    "author": "Anonymous",
    "date": "2026-10-02",
}
EXPECTED_OPENING = (
    "*This supplement is a typeset rendering of `darts/LEMMAS.md` in the research archive; "
    "that file is the reviewed source and governs in case of any discrepancy.*"
)

OUT: list[str] = []


def say(line: str = "") -> None:
    OUT.append(line)


# ----------------------------------------------------------------------------------------------
# Source-side (Unicode) normalisation
# ----------------------------------------------------------------------------------------------
SUB = {"₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4", "₅": "5", "₆": "6", "₇": "7", "₈": "8", "₉": "9",
       "₊": "+", "₋": "-"}
SUP_ = {"⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9",
        "⁻": "-", "⁺": "+", "ⁱ": "i", "ᴬ": "A", "ᶜ": "C", "ᵀ": "T"}
# ᶜ (MODIFIER LETTER SMALL C) is the source's superscript for the vertex C (Unicode has no capital C
# superscript); it is rendered as ^C, next to ᴬ -> ^A.

UNI = {
    "≤": r"\le", "≥": r"\ge", "−": "-", "λ": r"\lambda", "ρ": r"\rho", "π": r"\pi", "θ": r"\theta",
    "∈": r"\in", "⊂": r"\subset", "Ω": r"\Omega", "∠": r"\angle", "ℓ": r"\ell", "τ": r"\tau",
    "°": r"\circ", "μ": r"\mu", "‖": r"\|", "∪": r"\cup", "σ": r"\sigma", "→": r"\to", "ξ": r"\xi",
    "β": r"\beta", "∂": r"\partial", "η": r"\eta", "φ": r"\varphi", "Σ": r"\Sigma", "·": r"\cdot",
    "∫": r"\int", "∎": r"\square", "∩": r"\cap", "Π": r"\Pi", "ν": r"\nu", "α": r"\alpha",
    "⟨": r"\langle", "⟩": r"\rangle", "Γ": r"\Gamma", "′": "'", "∇": r"\nabla", "δ": r"\delta",
    "×": r"\times", "…": r"\dots", "𝓕": "F", "ℝ": "R", "ε": r"\varepsilon", "Δ": r"\Delta",
    "√": r"\sqrt", "γ": r"\gamma", "∞": r"\infty", "⊥": r"\perp", "𝓑": "B", "∖": r"\setminus",
    "≠": r"\ne", "⊃": r"\supset", "±": r"\pm", "∓": r"\mp", "⊇": r"\supseteq", "↦": r"\mapsto",
    "⇒": r"\Rightarrow", "⊕": r"\oplus", "∅": r"\emptyset", "∝": r"\propto", "κ": r"\kappa",
    "↔": r"\leftrightarrow",
}
COMBINING = {"\u0304": r"\bar", "\u0302": r"\hat", "\u0301": r"\acute"}

# ----------------------------------------------------------------------------------------------
# Supplement-side (LaTeX) normalisation
# ----------------------------------------------------------------------------------------------
DROP_CMDS = {"mathrm", "operatorname", "text", "mathsf", "mathbb", "mathcal", "bigl", "bigr", "Bigl",
             "Bigr", "big", "Big", "qquad", "quad", "left", "right", "limits", "displaystyle",
             "newline"}  # \newline only moves the date of the first header to its own line
WORD_CMDS = {"sin", "cos", "tan", "cot", "arccos", "arctan", "det", "inf", "sup", "min", "max"}
RENAME_CMDS = {"sum": [r"\Sigma"], "ldots": [r"\dots"], "dots": [r"\dots"], "iint": [r"\int", r"\int"]}

KEEP_PUNCT = set("()[]|<>=+-/*'\",;:.§–—")
SKIP_CHARS = set("{}$&^_`#~\\")

# Notation canonicalisations applied to the supplement before tokenising (listed in the log).
SMALLMATRIX = re.compile(
    r"\\bigl\[\\begin\{smallmatrix\}(.+?)\s*&\s*(.+?)\s*\\\\\s*(.+?)\s*&\s*(.+?)\\end\{smallmatrix\}\\bigr\]")


def sqrt_parens(s: str) -> str:
    """\\sqrt{X} -> \\sqrt(X) when X is compound (contains + or -), as the source writes √(a² + 4b²);
    a simple argument (\\sqrt{\\lambda_k}) keeps no parentheses, as the source writes √λ_k."""
    out, i = [], 0
    while True:
        j = s.find("\\sqrt{", i)
        if j < 0:
            out.append(s[i:])
            return "".join(out)
        k, depth = j + 6, 1
        while depth:
            depth += {"{": 1, "}": -1}.get(s[k], 0)
            k += 1
        arg = s[j + 6:k - 1]
        out.append(s[i:j])
        out.append("\\sqrt(" + arg + ")" if re.search(r"[+-]", arg) else "\\sqrt{" + arg + "}")
        i = k


CANON_RULES = [
    ("2x2 smallmatrix [a b; c d] -> source list notation [[a, b],[c, d]]",
     lambda s: SMALLMATRIX.sub(lambda m: f"[[{m.group(1)}, {m.group(2)}],[{m.group(3)}, {m.group(4)}]]", s)),
    ("range 'k = 0, \\dots, 30' -> source 'k = 0..30' (source '..' is read as an ellipsis)",
     lambda s: s.replace(r"k = 0, \dots, 30", r"k = 0\dots 30")),
    ("\\sqrt{compound argument} -> \\sqrt(argument), as in the source's √(...)",
     sqrt_parens),
]

NUM = r"\d+(?:\.\d+)?"
TOKEN_RE = re.compile(
    r"(?P<dbl>\\\\)|(?P<benv>\\(?:begin|end)\{[A-Za-z*]+\})|(?P<cmd>\\[A-Za-z]+)|(?P<esc>\\[^A-Za-z])"
    rf"|(?P<num>{NUM})|(?P<ws>\s+)|(?P<ch>.)",
    re.S,
)


def nfd(text: str) -> str:
    """NFD (separates ē/é into base + accent) but keeps ≠ as one symbol."""
    return unicodedata.normalize("NFD", text).replace("=\u0338", "≠")


def tokenize_source(text: str) -> list[str]:
    text = nfd(text)
    # script runs -> separated ASCII
    out = []
    i = 0
    while i < len(text):
        c = text[i]
        table = SUB if c in SUB else SUP_ if c in SUP_ else None
        if table is None:
            out.append(c)
            i += 1
            continue
        j = i
        while j < len(text) and text[j] in table:
            j += 1
        out.append(" " + "".join(table[x] for x in text[i:j]) + " ")
        i = j
    text = "".join(out)
    text = re.sub(r"\.{2,}", r" \\dots ", text)      # '0..30' and '...' are ellipses
    text = re.sub(r"(?<= )\\(?= )", r"\\setminus", text)  # 'Q' \ σ' (ASCII set minus, one place)
    toks: list[str] = []
    for m in TOKEN_RE.finditer(text):
        kind = m.lastgroup
        s = m.group()
        if kind == "num":
            toks.append(s)
        elif kind == "cmd":  # the source contains no LaTeX except the \dots/\setminus inserted above
            toks.append(s)
        elif kind in ("ws", "dbl", "benv", "esc"):
            continue
        else:
            c = s
            if c in COMBINING:
                toks.insert(len(toks) - 1 if toks else 0, COMBINING[c])  # accent before its base
            elif c in UNI:
                mapped = UNI[c]
                if mapped.startswith("\\") or mapped in KEEP_PUNCT:
                    toks.append(mapped)
                else:
                    toks.extend(list(mapped))
            elif c in SKIP_CHARS:
                continue
            elif c.isalpha():
                toks.append(c)
            elif c in KEEP_PUNCT:
                toks.append(c)
            else:
                toks.append("<?" + c + ">")
    return toks


def tokenize_latex(text: str) -> list[str]:
    text = nfd(text)
    for _, rule in CANON_RULES:
        text = rule(text)
    toks: list[str] = []
    for m in TOKEN_RE.finditer(text):
        kind = m.lastgroup
        s = m.group()
        if kind in ("ws", "dbl", "benv"):
            continue
        if kind == "num":
            toks.append(s)
        elif kind == "cmd":
            name = s[1:]
            if name in DROP_CMDS:
                continue
            if name in WORD_CMDS:
                toks.extend(list(name))
            elif name in RENAME_CMDS:
                toks.extend(RENAME_CMDS[name])
            else:
                toks.append("\\" + name)
        elif kind == "esc":
            c = s[1]
            if c == "|":
                toks.append(r"\|")
            elif c in ",;:! {}":
                continue
            elif c in KEEP_PUNCT:
                toks.append(c)
            # '\_' and other escapes of skipped characters vanish
        else:
            c = s
            if c in COMBINING:
                toks.insert(len(toks) - 1 if toks else 0, COMBINING[c])
            elif c in UNI:
                mapped = UNI[c]
                toks.append(mapped) if (mapped.startswith("\\") or mapped in KEEP_PUNCT) else toks.extend(list(mapped))
            elif c in SKIP_CHARS:
                continue
            elif c.isalpha():
                toks.append(c)
            elif c in KEEP_PUNCT:
                toks.append(c)
            else:
                toks.append("<?" + c + ">")
    return toks


def is_num(t: str) -> bool:
    return re.fullmatch(NUM, t) is not None


# ----------------------------------------------------------------------------------------------
# Document splitting
# ----------------------------------------------------------------------------------------------
HEADER_RE = re.compile(r"^(#{1,6})\s+(.*\S)\s*$")


def split_sections(text: str) -> list[tuple[str, str]]:
    """[(header_line, section_text_including_header)] -- text before the first header is dropped."""
    sections: list[tuple[str, list[str]]] = []
    for line in text.splitlines():
        if HEADER_RE.match(line):
            sections.append((line, [line]))
        elif sections:
            sections[-1][1].append(line)
    return [(h, "\n".join(lines)) for h, lines in sections]


def supplement_parts(text: str) -> tuple[dict, str, str]:
    m = re.match(r"\A---\n(.*?)\n---\n(.*)\Z", text, re.S)
    if not m:
        return {}, "", text
    yaml_block, rest = m.group(1), m.group(2)
    fields = {}
    for line in yaml_block.splitlines():
        km = re.match(r'^(\w+):\s*"?(.*?)"?\s*$', line)
        if km:
            fields[km.group(1)] = km.group(2)
    idx = re.search(r"^#\s", rest, re.M)
    pre, body = (rest[: idx.start()], rest[idx.start():]) if idx else (rest, "")
    return fields, pre.strip(), body


# ----------------------------------------------------------------------------------------------
# Context helpers
# ----------------------------------------------------------------------------------------------
def lines_with_number(text: str, literal: str, tokenizer) -> list[str]:
    hits = []
    for no, line in enumerate(text.splitlines(), 1):
        if literal in [t for t in tokenizer(line) if is_num(t)]:
            hits.append(f"      l.{no}: {line.strip()[:170]}")
    return hits


def show(tokens: list[str]) -> str:
    return " ".join(tokens)


# ----------------------------------------------------------------------------------------------
def main() -> int:
    src = SRC.read_text(encoding="utf-8")
    sup = SUP.read_text(encoding="utf-8")
    say("supplement_check.py -- mechanical fidelity check")
    say(f"source     : {SRC}  sha256={hashlib.sha256(src.encode()).hexdigest()}")
    say(f"supplement : {SUP}  sha256={hashlib.sha256(sup.encode()).hexdigest()}")
    say("")
    verdict: dict[str, bool] = {}

    # (0) front matter and opening paragraph
    fields, opening, body = supplement_parts(sup)
    ok0 = fields == EXPECTED_YAML and opening == EXPECTED_OPENING
    say("(0) front matter and opening paragraph")
    say(f"    YAML fields: {fields}")
    say(f"    opening paragraph matches the specified text: {opening == EXPECTED_OPENING}")
    say(f"    => {'PASS' if ok0 else 'FAIL'}")
    verdict["0"] = ok0
    say("")

    src_tok = tokenize_source(src)
    sup_tok = tokenize_latex(body)

    # (a) numeric literals
    src_nums = [t for t in src_tok if is_num(t)]
    sup_nums = [t for t in sup_tok if is_num(t)]
    cs, cp = collections.Counter(src_nums), collections.Counter(sup_nums)
    say("(a) numeric literals")
    say(f"    source: {len(src_nums)} literals ({len(cs)} distinct); supplement body: {len(sup_nums)} literals ({len(cp)} distinct)")
    only_src, only_sup = cs - cp, cp - cs
    ok_a = not only_src and not only_sup
    for lit, n in sorted(only_src.items()):
        say(f"    MISSING in supplement: {lit!r} x{n}")
        say("\n".join(lines_with_number(src, lit, tokenize_source)) or "      (no line context)")
    for lit, n in sorted(only_sup.items()):
        say(f"    EXTRA in supplement: {lit!r} x{n}")
        say("\n".join(lines_with_number(body, lit, tokenize_latex)) or "      (no line context)")
    ordered_same = src_nums == sup_nums
    say(f"    multiset identical: {ok_a}; ordered sequence identical: {ordered_same}")
    if ok_a and not ordered_same:
        sm = difflib.SequenceMatcher(None, src_nums, sup_nums, autojunk=False)
        for op, i1, i2, j1, j2 in sm.get_opcodes():
            if op != "equal":
                say(f"    order difference {op}: source {src_nums[max(0,i1-3):i2+3]} vs supplement {sup_nums[max(0,j1-3):j2+3]}")
    ok_a = ok_a and ordered_same
    say(f"    => {'PASS' if ok_a else 'FAIL'}")
    verdict["a"] = ok_a
    say("")

    # (b) headers and labels
    say("(b) section headers and statement labels")
    src_secs = split_sections(src)
    sup_secs = split_sections(body)
    src_heads = [(len(HEADER_RE.match(h).group(1)), tuple(tokenize_source(HEADER_RE.match(h).group(2)))) for h, _ in src_secs]
    sup_heads = [(len(HEADER_RE.match(h).group(1)), tuple(tokenize_latex(HEADER_RE.match(h).group(2)))) for h, _ in sup_secs]
    ok_h = src_heads == sup_heads
    say(f"    headers: source {len(src_heads)}, supplement {len(sup_heads)}; identical (level + token skeleton): {ok_h}")
    for (hs, _), (hp, _) in zip(src_secs, sup_secs):
        mark = "  " if (len(HEADER_RE.match(hs).group(1)), tuple(tokenize_source(HEADER_RE.match(hs).group(2)))) == \
                       (len(HEADER_RE.match(hp).group(1)), tuple(tokenize_latex(HEADER_RE.match(hp).group(2)))) else "!!"
        say(f"    {mark} {hs}")
        if mark == "!!":
            say(f"       vs {hp}")
    label_re = re.compile(r"\b(Lemmas?|Theorem|Proposition|Prop|Corollary|Remarks?)\b(?:\s+(\d+(?:\.\d+)*|[A-Z](?![A-Za-z])))?")
    src_labels = [" ".join(x for x in m.groups() if x) for m in label_re.finditer(src)]
    sup_labels = [" ".join(x for x in m.groups() if x) for m in label_re.finditer(body)]
    ok_l = src_labels == sup_labels
    say(f"    labels: source {len(src_labels)}, supplement {len(sup_labels)}; ordered lists identical: {ok_l}")
    if not ok_l:
        for line in difflib.unified_diff(src_labels, sup_labels, "LEMMAS.md", "SUPPLEMENT_S.md", lineterm="", n=2):
            say("      " + line)
    say("    label sequence: " + "; ".join(src_labels))
    ok_b = ok_h and ok_l
    say(f"    => {'PASS' if ok_b else 'FAIL'}")
    verdict["b"] = ok_b
    say("")

    # (c) code names
    say("(c) backticked code names")
    code_re = re.compile(r"`([^`\n]+)`")
    src_code = collections.Counter(code_re.findall(src))
    sup_code = collections.Counter(code_re.findall(body))
    missing = [c for c in src_code if c not in sup_code]
    ok_c = not missing
    say(f"    source: {sum(src_code.values())} spans ({len(src_code)} distinct); supplement body: {sum(sup_code.values())} spans ({len(sup_code)} distinct)")
    for c in missing:
        say(f"    MISSING: `{c}`")
    say(f"    multisets identical: {src_code == sup_code}")
    if src_code != sup_code:
        say(f"    count differences: source-only {dict(src_code - sup_code)}, supplement-only {dict(sup_code - src_code)}")
    say(f"    => {'PASS' if ok_c else 'FAIL'}")
    verdict["c"] = ok_c
    say("")

    # (d) token skeleton, section by section
    say("(d) token skeleton (letters, numbers, named symbols, punctuation), section by section")
    say("    notation canonicalisations applied to the supplement before tokenising:")
    counts = [len(SMALLMATRIX.findall(body)), body.count(r"k = 0, \dots, 30"),
              sqrt_parens(body).count("\\sqrt(") - body.count("\\sqrt(")]
    for (desc, _), n_hits in zip(CANON_RULES, counts):
        say(f"      - {desc}  [{n_hits} occurrence(s)]")
    ok_d = len(src_secs) == len(sup_secs)
    n_diff = 0
    total_src = total_sup = 0
    for (hs, ts), (hp, tp) in zip(src_secs, sup_secs):
        a, b = tokenize_source(ts), tokenize_latex(tp)
        total_src += len(a)
        total_sup += len(b)
        if a == b:
            continue
        ok_d = False
        sm = difflib.SequenceMatcher(None, a, b, autojunk=False)
        say(f"    section {hs[:70]!r}")
        for op, i1, i2, j1, j2 in sm.get_opcodes():
            if op == "equal":
                continue
            n_diff += 1
            say(f"      {op}: source[{show(a[i1:i2])}] supplement[{show(b[j1:j2])}]")
            say(f"        context source    : ...{show(a[max(0, i1 - 12):i2 + 12])}...")
            say(f"        context supplement: ...{show(b[max(0, j1 - 12):j2 + 12])}...")
    unknown = sorted({t for t in src_tok + sup_tok if t.startswith("<?")})
    if unknown:
        ok_d = False
        say(f"    unmapped characters: {unknown}")
    say(f"    tokens: source {total_src}, supplement {total_sup}; differing blocks: {n_diff}")
    say(f"    => {'PASS' if ok_d else 'FAIL'}")
    verdict["d"] = ok_d
    say("")

    # (e) pdflatex safety
    say("(e) pdflatex safety (no Unicode inside math; only — – § é … outside math/code)")
    allowed_text = set("—–§é…")
    stripped = re.sub(r"`[^`\n]*`", "", body)
    math_spans = re.findall(r"\$\$.+?\$\$|\$[^$\n]+?\$", stripped, re.S)
    bad_math = sorted({c for span in math_spans for c in span if ord(c) > 127})
    text_only = re.sub(r"\$\$.+?\$\$|\$[^$\n]+?\$", " ", stripped, flags=re.S)
    bad_text = sorted({c for c in text_only if ord(c) > 127 and c not in allowed_text})
    ok_e = not bad_math and not bad_text
    say(f"    math spans: {len(math_spans)}; non-ASCII inside math: {bad_math}; disallowed non-ASCII in text: {bad_text}")
    say(f"    => {'PASS' if ok_e else 'FAIL'}")
    verdict["e"] = ok_e
    say("")

    overall = all(verdict.values())
    say("SUMMARY: " + ", ".join(f"({k}) {'PASS' if v else 'FAIL'}" for k, v in verdict.items()))
    say(f"OVERALL: {'PASS' if overall else 'FAIL'}")
    text = "\n".join(OUT) + "\n"
    sys.stdout.write(text)
    LOG.write_text(text, encoding="utf-8")
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
