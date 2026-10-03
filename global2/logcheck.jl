# logcheck.jl — fail-closed loading of cover logs and hash binding of verification results.
#
# load_log_strict(logf, schema) parses EVERY line with the strict JSON parser, applies the per-status schema
# (status whitelist; exact key sets; field types and ranges), detects duplicate paths, and returns
# (recs, rawsha, problems): recs[path] = record (numbers Float64, arrays Vector{Float64}, null -> NaN),
# rawsha[path] = SHA-256 of the record's exact line text (the canonical record), problems = list of strings.
# Unknown statuses (including FAIL, ERROR, pilot_pass, pending) are problems.  Callers must refuse to certify
# anything when problems is non-empty.
using SHA
include("strictjson.jl")

sha_hex(s::AbstractString) = bytes2hex(sha256(String(s)))
file_sha(path) = open(io -> bytes2hex(sha256(io)), path)

"""
    verifier_identity(files; extra) -> (combined, perfile)

SHA-256 of every file, and a combined SHA-256 over the sorted lines "name<TAB>sha" plus the `extra` strings
(Julia and package versions).  Results are valid only for the exact verifier identity under which they were made.
"""
function verifier_identity(files::Vector{String}; extra::Vector{String} = String[])
    perfile = Dict(basename(f) => file_sha(f) for f in files)
    lines = sort([k * "\t" * v for (k, v) in perfile])
    combined = sha_hex(join(vcat(lines, extra), "\n"))
    return combined, perfile
end

_isnum(x) = x isa Float64 && isfinite(x)
_isint(x, vals) = _isnum(x) && isinteger(x) && Int(x) in vals
_isbox(lo, hi) = lo isa Vector && hi isa Vector && length(lo) == 4 && length(hi) == 4 &&
                 all(_isnum, lo) && all(_isnum, hi) && all(lo[i] < hi[i] for i = 1:4)
"compound tree grammar: leaf = E1|E2|E3|E3T|E4, node = [abcd](tree,tree)"
function _tree_ok(t::AbstractString)
    t in ("E1", "E2", "E3", "E3T", "E4") && return true
    (length(t) >= 6 && t[1] in "abcd" && t[2] == '(' && t[end] == ')') || return false
    inner = t[3:end-1]; depth = 0; cut = 0
    for (i, ch) in enumerate(inner)
        ch == '(' && (depth += 1); ch == ')' && (depth -= 1)
        depth < 0 && return false
        if ch == ',' && depth == 0; cut = i; break; end
    end
    cut == 0 && return false
    return _tree_ok(inner[1:cut-1]) && _tree_ok(inner[cut+1:end])
end
_tonative(v) = v === nothing ? NaN : (v isa Vector && all(x -> x isa Float64, v)) ? Float64.(v) :
               v isa Dict ? Dict{String,Any}(k => _tonative(w) for (k, w) in v) : v

const KEYS_EXCL = Set(["hi", "lo", "path", "status", "t"])
const KEYS_SPLIT = Set(["dim", "path", "status"])
const KEYS_CERT_T = Set(["L0", "a11up", "beta", "bound", "cert", "dlo", "fo", "g1lo", "g2lo", "hi", "kappa_cert", "kappa_f64",
                         "lo", "mode", "n", "path", "shift", "status", "t", "tF", "ups", "which"])
const KEYS_CERT_S = Set(["L0", "a11up", "beta", "g1lo", "g2lo", "hi", "kappa_cert", "kpt", "lo", "mode", "n", "path", "status",
                         "t", "tF", "ups", "which", "which_f64"])
const RX_BOUND = r"^(const|boot|KS|strip[1-4]|rect[1-7]|AC\+[A-Za-z0-9]+|bootTM2?:[A-Za-z0-9+]+)$"
_numornull(x) = x === nothing || _isnum(x)

"schema of the transport log (cover_t.jl); returns an error string or nothing"
function schema_transport(r::Dict{String,Any})
    st = r["status"]; ks = Set(keys(r))
    occursin(r"^T([1-9]|1[0-9]|2[0-4])(\.[abcd][01])*$", r["path"]) || return "bad path"
    if st == "split"
        ks == KEYS_SPLIT || return "split: keys $(sort(collect(ks)))"
        _isint(r["dim"], 1:4) || return "split: dim"
    elseif st in ("E1", "E2", "E3", "E3T", "E4", "SHELL")
        ks == KEYS_EXCL || return "$st: keys $(sort(collect(ks)))"
        _isbox(r["lo"], r["hi"]) || return "$st: box"
    elseif st == "COMPOUND"
        ks == union(KEYS_EXCL, Set(["tree"])) || return "COMPOUND: keys"
        _isbox(r["lo"], r["hi"]) || return "COMPOUND: box"
        (r["tree"] isa String && _tree_ok(r["tree"])) || return "COMPOUND: tree"
    elseif st == "certified"
        ks == KEYS_CERT_T || return "certified: keys $(sort(collect(ks)))"
        _isbox(r["lo"], r["hi"]) || return "certified: box"
        r["cert"] == "transport1" || return "certified: cert"
        r["fo"] isa Bool || return "certified: fo"
        _isint(r["n"], (8, 12)) || return "certified: n"
        _isint(r["shift"], (0, 1)) || return "certified: shift"
        r["mode"] in ("LG", "MONO") || return "certified: mode"
        (r["bound"] isa String && occursin(RX_BOUND, r["bound"])) || return "certified: bound"
        (_isnum(r["kappa_cert"]) && 0 < r["kappa_cert"] <= 2.5) || return "certified: kappa_cert"
        (_numornull(r["beta"]) && _numornull(r["L0"]) && _numornull(r["g1lo"]) && _numornull(r["g2lo"]) &&
         _isnum(r["a11up"]) && _isnum(r["kappa_f64"])) || return "certified: numbers"
        r["which"] in ("tr", "det", "zsz", "s11", "s22") || return "certified: which"
        r["ups"] isa Dict || return "certified: ups"
    else
        return "status not allowed: " * st
    end
    return nothing
end
"schema of the shell log (cover.jl)"
function schema_shell(r::Dict{String,Any})
    st = r["status"]; ks = Set(keys(r))
    occursin(r"^R[0-9]+(\.[abcd][01])*$", r["path"]) || return "bad path"
    if st == "split"
        ks == KEYS_SPLIT || return "split: keys"
        _isint(r["dim"], 1:4) || return "split: dim"
    elseif st in ("E1", "E2", "E3", "E4")
        ks == KEYS_EXCL || return "$st: keys"
        _isbox(r["lo"], r["hi"]) || return "$st: box"
    elseif st == "certified"
        ks == KEYS_CERT_S || return "certified: keys $(sort(collect(ks)))"
        _isbox(r["lo"], r["hi"]) || return "certified: box"
        _isint(r["n"], 4:16) || return "certified: n"
        r["mode"] in ("LG", "MONO") || return "certified: mode"
        (_isnum(r["beta"]) && _isnum(r["L0"]) && _isnum(r["kappa_cert"]) && 0 < r["kappa_cert"] <= 2.5) || return "certified: numbers"
        r["which"] in ("tr", "det", "zsz", "s11", "s22") || return "certified: which"
    else
        return "status not allowed: " * st
    end
    return nothing
end

"""
    load_log_strict(logf, schema) -> (recs, rawsha, problems)
"""
function load_log_strict(logf::AbstractString, schema::Function)
    recs = Dict{String,Dict{String,Any}}(); rawsha = Dict{String,String}(); problems = String[]
    for (k, ln) in enumerate(eachline(logf))
        d = try
            strict_json_object(ln)
        catch err
            push!(problems, "line $k: unparsable ($(sprint(showerror, err)))"); continue
        end
        if !(get(d, "path", nothing) isa String && get(d, "status", nothing) isa String)
            push!(problems, "line $k: no path/status"); continue
        end
        msg = schema(d)
        msg === nothing || (push!(problems, "line $k ($(d["path"])): $msg"); continue)
        p = d["path"]
        if haskey(recs, p)
            push!(problems, "line $k: duplicate path $p"); continue
        end
        recs[p] = Dict{String,Any}(kk => _tonative(v) for (kk, v) in d)
        rawsha[p] = sha_hex(ln)
    end
    return recs, rawsha, problems
end
