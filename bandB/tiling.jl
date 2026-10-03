# tiling.jl — exact tiling check of a chart domain by leaf paths (used by replay.jl and tests.jl).  Author: Anonymous.
#
# A set of leaf paths tiles the root box exactly iff, in the binary trie of the paths, (i) no leaf is a prefix of
# another leaf (no overlaps) and (ii) every internal node has exactly the two children 'lower half of d' and
# 'upper half of d' for ONE coordinate d (nothing missing).  Then, by induction from the leaves, every internal node's
# box is the union of its two children's boxes, which meet only in a face (measure zero), so the leaves cover the root
# and have disjoint interiors.  Checked on the sorted path list by recursion over contiguous ranges.

lowerof(ch::Char) = 'A' <= ch <= 'D' ? Char(ch - 'A' + 'a') : ch
ispair(c1::Char, c2::Char) = c1 != c2 && lowerof(c1) == lowerof(c2)

function check_range(L::Vector{String}, i::Int, j::Int, n::Int)
    # all L[i..j] share the first n characters; L sorted, no duplicates
    if i == j && length(L[i]) == n
        return true, ""
    end
    length(L[i]) == n && return false, "leaf $(L[i]) overlaps deeper leaves"
    c1 = L[i][n+1]; c2 = L[j][n+1]
    ispair(c1, c2) || return false, "node $(L[i][1:n]) has children '$c1'..'$c2' (missing sibling or mixed directions)"
    m = i
    while m < j && L[m+1][n+1] == c1
        m += 1
    end
    all(k -> L[k][n+1] == c2, m+1:j) || return false, "node $(L[i][1:n]) has more than two child letters"
    ok, msg = check_range(L, i, m, n + 1); ok || return ok, msg
    return check_range(L, m + 1, j, n + 1)
end

"""
    check_tiling(chart, leaves) -> (ok, message)

Exact combinatorial tiling test (see header) plus the redundant exact volume identity Σ vol(leaf) = vol(root).
"""
function check_tiling(chart::Int, leaves::Vector{String})
    isempty(leaves) && return false, "no leaves"
    for p in leaves
        all(c -> ('a' <= c <= 'd') || ('A' <= c <= 'D'), p) || return false, "bad path $p"
    end
    L = sort(leaves)
    for k = 2:length(L)
        L[k] == L[k-1] && return false, "duplicate leaf $(L[k])"
    end
    ok, msg = check_range(L, 1, length(L), 0)
    ok || return false, msg
    vol = sum(volume(box_from_path(chart, p)) for p in L)
    vol == volume(root_box(chart)) || return false, "volume mismatch $(vol) vs $(volume(root_box(chart)))"
    return true, "ok: $(length(L)) leaves tile the chart-$chart root exactly"
end
