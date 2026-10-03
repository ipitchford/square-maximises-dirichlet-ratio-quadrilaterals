"""Rigour lint of the verifier sources (adversarial review 2026-10-02, gate 1).  Exit status 1 on any unannotated hazard.

Sources linted: the transitive local-import closure (binding.closure) of every replay / accounting / driver entry point:
  replay2.py replay3.py accounting.py check_charts.py check_cbox_cover.py hull_to_box.py lemma_sector.py lemma_longarm.py
  merge_pool.py replay_ok.py remanifest.py rigor_lint.py replay_all.py bootstrap.py test_containment.py make_pin.py summarize.py
Flagged constructs (AST walk):
  F1  float(...), fin(...), round(...) calls and .mid() method calls      (a nearest conversion that could reach a decision)
  F2  x ** k with a non-literal base                                        (python-flint ball ** k is NaN for balls around 0)
  F3  builtin max(...) / min(...)                                           (ball comparisons can be undecided; NaN is dropped)
  F4  math.nextafter / np.nextafter / math.ulp outside nanguard.py          (endpoint conversions must use fdown / fup)
  F5  NumPy casts and reductions, under ANY alias (import numpy as x; from numpy import f [as g]; from numpy.linalg import ...):
      max/min/amax/amin/sum/prod/mean/var/std/all/any/maximum/minimum/round/rint/floor/ceil/float64/float32/asarray/array/
      linspace/percentile/median/fromiter/dot/matmul/einsum, ufunc methods .reduce/.accumulate/.reduceat/.outer on a NumPy
      attribute, and the methods .astype/.item/.max/.min/.sum/.prod/.mean/.var/.std/.all/.any   (rounding hidden in NumPy)
Every flagged node must be covered by an annotation, else it is a FAILURE:
  * on its own line, or on the first line of its enclosing statement:   # RIGOR: <tag>: <reason>
  * on the def line of an enclosing function:                           # RIGOR-FUNCTION: <tag>: <reason>
  * as a module comment (line starting with it):                        # RIGOR-MODULE: <tag>: <reason>
  Function- and module-level annotations are not blanket exemptions: every site they cover is listed individually in
  rigor_allowlist.json (file, scope, construct, sha256 of the source line); the lint fails unless the covered sites and the
  allowlist agree exactly (a new or changed site in an exempted scope must be audited and added with --update-allowlist).
  A RIGOR-MODULE file must be one of replay_common.EXEMPT_MODULES, which every replay checks were never loaded.
  tags:  heuristic  the value only steers a search or is a free choice (test vectors, shifts, split directions, float
                    pre-checks that can only cause a failure, node placement); it never decides a rigorous inequality
         info       logged / printed only
         exact      the operation is exact (float of a float or of an exactly representable rational; min/max over exact
                    floats or exact arb points; an exactly checked conversion)
         accounted  the rounding is covered by an explicit error term stated in the comment / docstring (cite it)
         helper     inside the audited directed-rounding helpers (nanguard.py)
  A reason (non-empty text after the tag) is required.
--list prints every annotation with its location (for review)."""
import ast, os, sys, re, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from binding import closure
ENTRIES = ['replay2.py', 'replay3.py', 'accounting.py', 'check_charts.py', 'check_cbox_cover.py', 'hull_to_box.py',
           'lemma_sector.py', 'lemma_longarm.py', 'merge_pool.py', 'replay_ok.py', 'remanifest.py', 'rigor_lint.py',
           'replay_all.py', 'bootstrap.py', 'test_containment.py', 'make_pin.py', 'summarize.py']
TAGS = ('heuristic', 'info', 'exact', 'accounted', 'helper')
import io, tokenize
NP_TAINT = {'max', 'min', 'amax', 'amin', 'nanmax', 'nanmin', 'sum', 'prod', 'mean', 'maximum', 'minimum', 'round', 'around', 'rint',
            'floor', 'ceil', 'float64', 'float32', 'float16', 'longdouble', 'asarray', 'array', 'linspace', 'percentile', 'median',
            'fromiter', 'nextafter', 'var', 'std', 'nanvar', 'nanstd', 'all', 'any', 'dot', 'matmul', 'einsum', 'nansum', 'cumsum',
            'cumprod', 'average', 'fsum'}
UFUNC_METHODS = {'reduce', 'accumulate', 'reduceat', 'outer'}
METHOD_TAINT = {'astype', 'item', 'max', 'min', 'sum', 'prod', 'mean', 'var', 'std', 'all', 'any'}
ALLOWLIST = os.path.join(HERE, 'rigor_allowlist.json')
ANN = re.compile(r'#\s*RIGOR(-FUNCTION|-MODULE)?:\s*([a-z]+)\s*:\s*(\S.*)$')

def annotations(src):
    """annotations found in real COMMENT tokens only (tokenize): a string literal that looks like an annotation is ignored"""
    out = {}
    for tok in tokenize.generate_tokens(io.StringIO(src).readline):
        if tok.type == tokenize.COMMENT:
            m = ANN.search(tok.string)
            if m: out[tok.start[0]] = (m.group(1) or '', m.group(2), m.group(3).strip())
    return out

def numpy_aliases(tree):
    """(module aliases of numpy and its submodules, {local name: numpy attribute} of names imported from numpy)"""
    mods = set(); names = {}
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            for a in n.names:
                if a.name == 'numpy' or a.name.startswith('numpy.'): mods.add(a.asname or a.name.split('.')[0])
        elif isinstance(n, ast.ImportFrom) and n.module and (n.module == 'numpy' or n.module.startswith('numpy.')):
            for a in n.names: names[a.asname or a.name] = a.name
    return mods, names

def _root(node):
    while isinstance(node, ast.Attribute): node = node.value
    return node.id if isinstance(node, ast.Name) else None

def flagged(tree):
    mods, npnames = numpy_aliases(tree)
    mods |= {'np', 'numpy'}
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            f = node.func
            if isinstance(f, ast.Name) and f.id in ('float', 'fin', 'round'): yield node, 'F1 ' + f.id + '()'
            elif isinstance(f, ast.Attribute) and f.attr == 'mid' and not node.args: yield node, 'F1 .mid()'
            elif isinstance(f, ast.Name) and f.id in ('max', 'min'): yield node, 'F3 ' + f.id + '()'
            elif isinstance(f, ast.Attribute) and f.attr in ('nextafter', 'ulp') and isinstance(f.value, ast.Name) and f.value.id in ('math', 'np', 'numpy'):
                yield node, 'F4 ' + f.value.id + '.' + f.attr
            elif isinstance(f, ast.Attribute) and _root(f) in mods and f.attr in NP_TAINT:
                yield node, 'F5 np.' + f.attr
            elif isinstance(f, ast.Attribute) and f.attr in UFUNC_METHODS and _root(f) in mods:
                yield node, 'F5 ufunc.' + f.attr
            elif isinstance(f, ast.Name) and f.id in npnames and npnames[f.id] in NP_TAINT:
                yield node, 'F5 ' + f.id + ' (from numpy)'
            elif isinstance(f, ast.Attribute) and f.attr in METHOD_TAINT and not (isinstance(f.value, ast.Name) and f.value.id in (mods | {'math'})):
                yield node, 'F5 .' + f.attr + '()'
        elif isinstance(node, ast.BinOp) and isinstance(node.op, ast.Pow):
            b = node.left
            lit = isinstance(b, ast.Constant) and isinstance(b.value, (int, float))
            if not lit: yield node, 'F2 **'

def lint_file(fn, listing):
    src = open(os.path.join(HERE, fn)).read(); lines = src.split('\n')
    tree = ast.parse(src)
    ann = annotations(src)
    mod = [v for k, v in ann.items() if v[0] == '-MODULE' and lines[k - 1].lstrip().startswith('#')]
    # statement start lines and function def lines for every node
    parents = {}
    for p in ast.walk(tree):
        for c in ast.iter_child_nodes(p): parents[c] = p
    def covering(node):
        lns = [node.lineno]
        p = node
        while p in parents:
            p = parents[p]
            if isinstance(p, ast.stmt): lns.append(p.lineno); break
        q = node; fdefs = []
        while q in parents:
            q = parents[q]
            if isinstance(q, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)) and hasattr(q, 'lineno'):
                fdefs.append((q.lineno, getattr(q, 'name', '<lambda>')))
        return lns, fdefs
    bad = []; used = []; blanket = []
    if mod:
        import replay_common
        if os.path.splitext(fn)[0] not in replay_common.EXEMPT_MODULES:
            bad.append((0, 'RIGOR-MODULE annotation in a module that is not in replay_common.EXEMPT_MODULES', fn))
    for node, kind in flagged(tree):
        lns, fdefs = covering(node)
        hit = None
        for l in lns:
            a = ann.get(l)
            if a and a[0] == '': hit = (l, a); break
        if hit is None:
            for l, fname in fdefs:
                a = ann.get(l)
                if a and a[0] == '-FUNCTION': hit = (l, a); blanket.append(dict(file=fn, scope='function:' + fname, kind=kind.split(' ', 1)[0],
                    line_sha256=hashlib.sha256(lines[node.lineno - 1].strip().encode()).hexdigest())); break
        if hit is None and mod:
            hit = (0, mod[0]); blanket.append(dict(file=fn, scope='module', kind=kind.split(' ', 1)[0],
                                                   line_sha256=hashlib.sha256(lines[node.lineno - 1].strip().encode()).hexdigest()))
        if hit is None: bad.append((node.lineno, kind, lines[node.lineno - 1].strip()[:110])); continue
        tag, reason = hit[1][1], hit[1][2]
        if tag not in TAGS or not reason: bad.append((node.lineno, kind + f' (bad annotation {tag!r})', lines[node.lineno - 1].strip()[:110])); continue
        if kind.startswith('F4') and fn != 'nanguard.py' and tag != 'exact':
            bad.append((node.lineno, kind + ' outside nanguard needs an exact-tag justification', lines[node.lineno - 1].strip()[:110])); continue
        used.append((node.lineno, kind, tag, reason))
    for k, v in ann.items():
        if v[1] not in TAGS or not v[2]: bad.append((k, f'malformed annotation tag {v[1]!r}', lines[k - 1].strip()[:110]))
    if listing:
        for (ln, kind, tag, reason) in sorted(set(used)):
            print(f'  {fn}:{ln}: {kind}: {tag}: {reason[:120]}')
    return bad, len(used), blanket

def main():
    import collections, json
    listing = '--list' in sys.argv
    files = closure(ENTRIES)
    nbad = 0; nann = 0; covered = []
    for fn in files:
        bad, n, blanket = lint_file(fn, listing); nann += n; covered += blanket
        for (ln, kind, txt) in bad:
            print(f'FAIL {fn}:{ln}: {kind}: {txt}'); nbad += 1
    key = lambda d: (d['file'], d['scope'], d['kind'], d['line_sha256'])
    if '--update-allowlist' in sys.argv:
        json.dump(sorted(covered, key=key), open(ALLOWLIST, 'w'), indent=0)
        print(f'allowlist written: {len(covered)} sites in function / module scopes')
    allow = json.load(open(ALLOWLIST)) if os.path.exists(ALLOWLIST) else []
    have = collections.Counter(key(d) for d in covered); want = collections.Counter(key(d) for d in allow)
    for k in sorted((have - want).elements()): print(f'FAIL unlisted site in an exempted scope: {k[0]} {k[1]} {k[2]} {k[3][:12]}'); nbad += 1
    for k in sorted((want - have).elements()): print(f'FAIL stale allowlist entry: {k[0]} {k[1]} {k[2]} {k[3][:12]}'); nbad += 1
    print(f'rigor_lint: {len(files)} files, {nann} annotated sites ({len(covered)} of them in function / module scopes, each listed in '
          f'rigor_allowlist.json), {nbad} failures')
    sys.exit(1 if nbad else 0)

if __name__ == '__main__':
    main()
