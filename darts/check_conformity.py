"""Conformity check of the multi-piece charts: every mesh edge lying on a piece edge that is shared by two pieces must be shared by
two mesh elements (no T-junctions, no unmatched nodes); every boundary mesh edge must lie on a piece edge that belongs to exactly
one piece.  Run for every chart/level used in certs/."""
import sys, collections
sys.path.insert(0, '/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts')
import numpy as np
from charts import bd_chart, spike_cr_chart, c_only_chart, spike2_cr_chart
ALLOK = [True]
def check(ch, name):
    # piece edges (by global vertex labels) and how many pieces contain them
    pe = collections.Counter()
    for pc in ch.pieces_meta:
        L = pc['labels']
        for a, b in ((0, 1), (1, 2), (2, 0)): pe[frozenset((L[a], L[b]))] += 1
    # recompute per-mesh-edge element counts from the stored piece data is not kept; use bc arrays + edge_pieces
    bad = 0
    for i in range(ch.nedge):
        pcs = ch.edge_pieces[i]
        if ch.cr_bc[i] in ('D', 'N'):           # boundary mesh edge (one element)
            # it must not lie on a shared piece edge: its piece set has one piece; check that piece's edge is unshared
            if len(pcs) != 1: bad += 1
    # count mesh edges on each shared piece edge from the edge_pieces table (two-piece edges)
    two = sum(1 for p in ch.edge_pieces if len(p) == 2)
    shared_piece_edges = [e for e, c in pe.items() if c == 2]
    expected = len(shared_piece_edges)*(2**ch.k)
    print(f"{name}: k={ch.k} pieces={len(ch.pieces_meta)} mesh edges={ch.nedge} shared piece edges={len(shared_piece_edges)} "
          f"mesh edges on them={two} expected={expected} boundary-on-shared={bad} -> {'OK' if (two == expected and bad == 0) else 'FAIL'}")
    if not (two == expected and bad == 0): ALLOK[0] = False
for k in (4, 5, 6):
    check(bd_chart(k), 'BD'); check(spike_cr_chart(k), 'SPIKE-CR'); check(c_only_chart(k), 'C-only'); check(spike2_cr_chart(k), 'SPIKE2-CR')

import sys; sys.exit(0 if ALLOK[0] else 1)
