"""General multi-piece affine charts with EXACT reference data.

A chart is a list of triangular pieces.  Piece q has exact rational reference vertices R_q = (r0, r1, r2) and a physical
vertex map p -> (P0(p), P1(p), P2(p)) given by polynomials in the parameters (see gcert.py).  The physical mesh of piece q is the
image of a fixed reference mesh of R_q under the affine map x -> P0 + E(p) R^{-1} (x - r0), E = [P1-P0, P2-P0], R = [r1-r0, r2-r0].
Reference meshes: k uniform refinements, optionally graded toward one vertex ('grade' = local vertex index) by s -> s^kappa along
rays from that vertex (kappa integer, so all coordinates stay dyadic rationals).
Gluing: pieces share edges, given by global vertex labels; along a shared edge the node fractions coincide (uniform, or graded
toward a common endpoint), so nodes are identified exactly by (global edge, fraction).
Boundary edges of the union carry a label 'D' (Dirichlet) or 'N' (Neumann) per family ('cr' and 'p2' may differ).
Matrices: for each piece and kl in {11, 12, 22}:  Khat^{kl}_q = int_ref d_k u d_l v (12: symmetrised), Mhat_q = int_ref u v,
assembled into global CR (edge dofs) and P2 (node + edge dofs) numbering, computed exactly and rounded once to doubles."""
import numpy as np, pickle, os, hashlib
from fractions import Fraction as Fr
from skfem import MeshTri
from refexact import grad_bary, p2_universal

CACHE = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'refcache')

# ---------------------------------------------------------------------------------------------- pinned chart caches (proof runs)
# In a proof run the environment variable DARTS_PIN names the driver's private copy of the reviewed pin list.  Then the cache
# directory must contain exactly the pinned refcache/* files with the pinned sha256 (verify_cache_inventory), a cache file is
# unpickled only after its bytes have been read ONCE and their sha256 equals the pin (no unreviewed pickle is ever executed), and
# nothing is written to the cache (an uncached chart is rebuilt in memory from the exact assembly).
def _trusted_caches():
    pin = os.environ.get('DARTS_PIN')
    if not pin: return None
    out = {}
    for line in open(pin, 'rb').read().decode().split('\n'):
        if not line.strip(): continue
        h, f = line.split('  ', 1)
        if f.startswith('refcache/'): out[f[len('refcache/'):]] = h
    return out
TRUSTED = _trusted_caches()
def verify_cache_inventory():
    """problems (empty = ok): the cache directory equals the pinned inventory, byte for byte"""
    if TRUSTED is None: return ['no pinned chart-cache inventory (DARTS_PIN unset)']
    have = sorted(os.listdir(CACHE)) if os.path.isdir(CACHE) else []
    bad = []
    if set(have) != set(TRUSTED): bad.append(f'cache inventory differs from the pin: extra {sorted(set(have) - set(TRUSTED))}, missing {sorted(set(TRUSTED) - set(have))}')
    for f in have:
        if f in TRUSTED and hashlib.sha256(open(os.path.join(CACHE, f), 'rb').read()).hexdigest() != TRUSTED[f]: bad.append(f'cache file {f} differs from the pin')
    return bad
def load_pinned(fn):
    """the unpickled cache state, read once and verified against the pin BEFORE unpickling; None if the file is not pinned"""
    name = os.path.basename(fn)
    if name not in TRUSTED: return None
    raw = open(fn, 'rb').read()
    if hashlib.sha256(raw).hexdigest() != TRUSTED[name]: raise ValueError(f'chart cache {name}: bytes differ from the pin (not unpickled)')
    return pickle.loads(raw)

def ref_unit_mesh(k):
    m = MeshTri.init_refdom().refined(k); n = 2**k
    P = np.rint(m.p*n).astype(np.int64)   # RIGOR: exact: the refined reference nodes are multiples of 1/n (dyadic), rint of exact integers
    return [(Fr(int(i), n), Fr(int(j), n)) for (i, j) in P.T], m.t.copy()

def piece_mesh(R, k, grade=None, kappa=2):
    """exact reference mesh of the triangle R = (r0, r1, r2) (tuples of Fractions).  Unit-triangle barycentrics (l1, l2) map to
    r0 + l1 (r1 - r0) + l2 (r2 - r0).  Grading toward local vertex g: s = 1 - lambda_g -> s^kappa."""
    upts, T = ref_unit_mesh(k)
    pts = []; bary = []
    for (u, v) in upts:
        l = [1 - u - v, u, v]
        if grade is not None:
            g = grade; s = 1 - l[g]
            if s != 0:
                sk = s**kappa; l = [li*sk/s for li in l]; l[g] = 1 - sk   # RIGOR: exact: Fraction ** int
        x = R[0][0] + l[1]*(R[1][0] - R[0][0]) + l[2]*(R[2][0] - R[0][0])
        y = R[0][1] + l[1]*(R[1][1] - R[0][1]) + l[2]*(R[2][1] - R[0][1])
        pts.append((x, y)); bary.append(tuple(l))
    # orient ccw
    T = T.copy()
    for e in range(T.shape[1]):
        a, b, c = (pts[T[i, e]] for i in range(3))
        d = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
        assert d != 0
        if d < 0: T[1, e], T[2, e] = T[2, e], T[1, e]
    return pts, bary, T

def cache_file(pieces, k, bc=None, kappa=2, tag=None):
    """the refcache file of a chart: a hash of its exact definition (pieces, level k, grading kappa, boundary conditions)"""
    bc = bc or {}
    key = tag or hashlib.sha1(repr((sorted((p['name'], tuple(map(str, p['ref'])), p['labels'], p['grade']) for p in pieces),
                                     k, kappa, sorted((tuple(sorted(e)), tuple(sorted(v.items()))) for e, v in bc.items()))).encode()).hexdigest()[:16]
    return os.path.join(CACHE, f'chart_{key}_k{k}.pkl')

class Chart:
    """pieces: list of dict(name, ref=(r0,r1,r2) Fractions, labels=(L0,L1,L2) global vertex labels, grade=local index or None)
    bc: dict mapping a frozenset of two global labels (a boundary edge of the union) -> {'cr': 'D'|'N', 'p2': 'D'|'N'}
    (unlisted boundary edges are Dirichlet for both)."""
    def __init__(self, pieces, k, bc=None, kappa=2, tag=None, fresh=False):
        """fresh=True: build from the exact assembly without reading or writing the cache (check_charts.py compares the two)"""
        self.pieces = pieces; self.pieces_meta = pieces; self.k = k; self.kappa = kappa; self.bc = bc or {}
        fn = cache_file(pieces, k, self.bc, kappa, tag)
        if fresh:
            self._build(); return
        if TRUSTED is not None:                                      # proof run: pinned caches only, no writes
            st = load_pinned(fn)
            if st is None: self._build(); return
            self.__dict__.update(st); self.pieces_meta = pieces; return
        os.makedirs(CACHE, exist_ok=True)
        if os.path.exists(fn):
            with open(fn, 'rb') as f: self.__dict__.update(pickle.load(f))
            self.pieces_meta = pieces; return
        self._build()
        with open(fn, 'wb') as f:
            pickle.dump(self.cached_state(), f)
    def cached_state(self):
        return {kk: v for kk, v in self.__dict__.items() if kk not in ('pieces', 'bc', 'pieces_meta')}

    def _build(self):
        nodes = {}          # key -> global node id ; keys: ('V', label) | ('E', frozenset(l_a,l_b), fraction from min label) | ('I', piece, i)
        node_xy = []        # (piece, reference coords) of a representative, for geometry
        pieces_data = []
        for q, pc in enumerate(self.pieces):
            R = tuple(tuple(Fr(c) for c in v) for v in pc['ref'])
            pts, bary, T = piece_mesh(R, self.k, pc['grade'], self.kappa)
            lab = pc['labels']; gid = []
            for i, l in enumerate(bary):
                nz = [j for j in range(3) if l[j] != 0]
                if len(nz) == 1: key = ('V', lab[nz[0]])
                elif len(nz) == 2:
                    a, b = nz; la, lb = lab[a], lab[b]
                    # fraction measured from the endpoint with the smaller label
                    if la < lb: key = ('E', (la, lb), l[b])
                    else: key = ('E', (lb, la), l[a])
                else: key = ('I', q, i)
                if key not in nodes: nodes[key] = len(nodes); node_xy.append((q, pts[i]))
                gid.append(nodes[key])
            pieces_data.append(dict(R=R, pts=pts, T=T, gid=np.array(gid, dtype=np.int64), bary=bary))   # RIGOR: exact: integer node ids
        self.nnode = len(nodes)
        # edges
        loc = [(1, 2), (2, 0), (0, 1)]
        eid = {}; edge_pieces = {}; edge_glabel = {}
        for q, pdd in enumerate(pieces_data):
            T = pdd['T']; g = pdd['gid']; ne = T.shape[1]
            t2e = np.empty((3, ne), dtype=np.int64)
            for e in range(ne):
                for i, (a, b) in enumerate(loc):
                    u, v = int(g[T[a, e]]), int(g[T[b, e]]); key = (u, v) if u < v else (v, u)
                    if key not in eid: eid[key] = len(eid)
                    t2e[i, e] = eid[key]; edge_pieces.setdefault(eid[key], set()).add(q)
                    # which reference edge of the piece (global label pair) contains this mesh edge (if on the piece boundary)
                    la = pdd['bary'][T[a, e]]; lb = pdd['bary'][T[b, e]]
                    zero = [j for j in range(3) if la[j] == 0 and lb[j] == 0]
                    if zero:
                        j = zero[0]; lab = self.pieces[q]['labels']
                        edge_glabel[eid[key]] = frozenset(lab[m] for m in range(3) if m != j)
            pdd['t2e'] = t2e
        self.nedge = len(eid)
        cnt = np.zeros(self.nedge, dtype=np.int64)
        for q, pdd in enumerate(pieces_data): np.add.at(cnt, pdd['t2e'].ravel(), 1)
        bnd = cnt == 1
        edge_nodes = np.empty((2, self.nedge), dtype=np.int64)
        for key, i in eid.items(): edge_nodes[:, i] = key
        # boundary condition per family
        def bc_of(i, fam):
            if not bnd[i]: return 'I'
            gl = edge_glabel.get(i)
            v = self.bc.get(gl, {}) if gl is not None else {}
            return v.get(fam, 'D')
        self.cr_bc = np.array([bc_of(i, 'cr') for i in range(self.nedge)])   # RIGOR: exact: boundary labels (strings)
        self.p2_edge_bc = np.array([bc_of(i, 'p2') for i in range(self.nedge)])   # RIGOR: exact: boundary labels (strings)
        # P2 node bc: a node is Dirichlet if it lies on a Dirichlet boundary edge
        p2_node_D = np.zeros(self.nnode, bool)
        for i in range(self.nedge):
            if self.p2_edge_bc[i] == 'D': p2_node_D[edge_nodes[:, i]] = True
        # CR node set irrelevant.  free dofs:
        self.crfree = np.where(self.cr_bc != 'D')[0]
        p2D = np.concatenate([p2_node_D, self.p2_edge_bc == 'D'])
        self.p2free = np.where(~p2D)[0]
        # piece-set of each CR dof (for residual grouping)
        self.edge_pieces = [tuple(sorted(edge_pieces[i])) for i in range(self.nedge)]
        # ---- exact matrices ----
        S2, Mu2, _ = p2_universal()
        self.cr = {}; self.p2 = {}; self.evecs = {}
        n2 = self.nnode + self.nedge
        for q, pdd in enumerate(pieces_data):
            T = pdd['T']; g = pdd['gid']; t2e = pdd['t2e']; pts = pdd['pts']; ne = T.shape[1]
            accC = {'11': {}, '12': {}, '22': {}}; mC = {}
            accP = {'11': {}, '12': {}, '22': {}, 'M': {}}
            ev = []
            for e in range(ne):
                verts = [pts[T[i, e]] for i in range(3)]
                gr, ar = grad_bary(*verts)
                for (u, v) in ((verts[0], verts[1]), (verts[1], verts[2]), (verts[2], verts[0])): ev.append((float(v[0]-u[0]), float(v[1]-u[1])))   # RIGOR: accounted: nearest float of an exact edge (|ev - e| <= u|ev|), covered in gcert.hmax_box
                for i in range(3):
                    di = int(t2e[i, e]); mC[di] = mC.get(di, 0) + ar/3
                    for j in range(3):
                        dj = int(t2e[j, e]); key = (di, dj)
                        for nm, val in (('11', 4*ar*gr[i][0]*gr[j][0]), ('22', 4*ar*gr[i][1]*gr[j][1]),
                                        ('12', 4*ar*(gr[i][0]*gr[j][1] + gr[i][1]*gr[j][0]))):
                            d = accC[nm]; d[key] = d.get(key, 0) + val
                dofs = [int(g[T[0, e]]), int(g[T[1, e]]), int(g[T[2, e]])] + [self.nnode + int(t2e[i, e]) for i in range(3)]
                for a in range(6):
                    for b in range(6):
                        v11 = 0; v12 = 0; v22 = 0
                        for m in range(3):
                            for n_ in range(3):
                                s_ = S2[a, b, m, n_]
                                if s_ == 0: continue
                                v11 += s_*gr[m][0]*gr[n_][0]; v22 += s_*gr[m][1]*gr[n_][1]
                                v12 += s_*(gr[m][0]*gr[n_][1] + gr[m][1]*gr[n_][0])
                        key = (dofs[a], dofs[b])
                        for nm, val in (('11', v11*ar), ('12', v12*ar), ('22', v22*ar), ('M', Mu2[a, b]*ar)):
                            d = accP[nm]; d[key] = d.get(key, 0) + val
            for nm, d in accC.items():
                self.cr[(q, nm)] = (np.array([kk[0] for kk in d], dtype=np.int64), np.array([kk[1] for kk in d], dtype=np.int64),   # RIGOR: accounted: integer indices; the values are the once-rounded exact entries (2.0000001u terms of gcert)
                                    np.array([float(x) for x in d.values()]))   # RIGOR: accounted: reference matrices rounded once to doubles; covered by the 2.0000001u / 3u terms of gcert
            md = np.zeros(self.nedge)
            for i, v in mC.items(): md[i] = float(v)   # RIGOR: accounted: lumped mass rounded once; covered by the 3u terms of gcert
            self.cr[(q, 'M')] = md
            for nm, d in accP.items():
                self.p2[(q, nm)] = (np.array([kk[0] for kk in d], dtype=np.int64), np.array([kk[1] for kk in d], dtype=np.int64),   # RIGOR: accounted: integer indices; the values are the once-rounded exact entries (2.0000001u terms of gcert)
                                    np.array([float(x) for x in d.values()]))   # RIGOR: accounted: reference matrices rounded once to doubles; covered by the 2.0000001u / 3u terms of gcert
            self.evecs[q] = np.array(ev)   # RIGOR: accounted: array of the nearest floats of the exact edges (gcert.hmax_box)
        self.Rref = [tuple(tuple(float(c) for c in v) for v in pdd['R']) for pdd in pieces_data]   # RIGOR: info: unused float copy (hmax_box uses Rref_exact)
        self.Rref_exact = [pdd['R'] for pdd in pieces_data]
        self.npiece = len(pieces_data); self.np2 = n2
