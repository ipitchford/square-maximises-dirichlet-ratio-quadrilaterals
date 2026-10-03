# strictjson.jl — a strict parser for the JSON subset written by cover.jl / cover_t.jl (fail closed).
#
# Accepts exactly one JSON object per line: objects, arrays, strings WITHOUT escape sequences, numbers, true, false,
# null.  Anything else (trailing garbage, duplicate keys, escapes, missing commas, NaN literals, ...) throws.
struct StrictJSONError <: Exception
    msg::String
end
Base.showerror(io::IO, e::StrictJSONError) = print(io, "StrictJSONError: ", e.msg)

mutable struct _SJ
    s::String
    i::Int
end
_peek(p::_SJ) = p.i <= ncodeunits(p.s) ? p.s[p.i] : '\0'
function _ws!(p::_SJ)
    while p.i <= ncodeunits(p.s) && p.s[p.i] in (' ', '\t')
        p.i += 1
    end
end
function _expect!(p::_SJ, c::Char)
    _ws!(p)
    _peek(p) == c || throw(StrictJSONError("expected '$c' at byte $(p.i)"))
    p.i += 1
end
function _value!(p::_SJ)
    _ws!(p)
    c = _peek(p)
    if c == '{'
        p.i += 1; d = Dict{String,Any}()
        _ws!(p)
        if _peek(p) == '}'
            p.i += 1; return d
        end
        while true
            _ws!(p)
            k = _string!(p)
            haskey(d, k) && throw(StrictJSONError("duplicate key \"$k\""))
            _expect!(p, ':')
            d[k] = _value!(p)
            _ws!(p)
            if _peek(p) == ','
                p.i += 1
            elseif _peek(p) == '}'
                p.i += 1; return d
            else
                throw(StrictJSONError("expected ',' or '}' at byte $(p.i)"))
            end
        end
    elseif c == '['
        p.i += 1; a = Any[]
        _ws!(p)
        if _peek(p) == ']'
            p.i += 1; return a
        end
        while true
            push!(a, _value!(p))
            _ws!(p)
            if _peek(p) == ','
                p.i += 1
            elseif _peek(p) == ']'
                p.i += 1; return a
            else
                throw(StrictJSONError("expected ',' or ']' at byte $(p.i)"))
            end
        end
    elseif c == '"'
        return _string!(p)
    elseif startswith(SubString(p.s, p.i), "true")
        p.i += 4; return true
    elseif startswith(SubString(p.s, p.i), "false")
        p.i += 5; return false
    elseif startswith(SubString(p.s, p.i), "null")
        p.i += 4; return nothing
    else
        m = match(r"^-?(0|[1-9][0-9]*)(\.[0-9]+)?([eE][-+]?[0-9]+)?", SubString(p.s, p.i))
        m === nothing && throw(StrictJSONError("unexpected character at byte $(p.i)"))
        p.i += ncodeunits(m.match)
        x = parse(Float64, m.match)
        isfinite(x) || throw(StrictJSONError("non-finite number"))
        return x
    end
end
function _string!(p::_SJ)
    _peek(p) == '"' || throw(StrictJSONError("expected string at byte $(p.i)"))
    j = findnext('"', p.s, p.i + 1)
    j === nothing && throw(StrictJSONError("unterminated string"))
    t = p.s[p.i+1:j-1]
    occursin('\\', t) && throw(StrictJSONError("escape sequences are not allowed"))
    p.i = j + 1
    return String(t)
end
"parse one line holding exactly one JSON object (throws StrictJSONError otherwise)"
function strict_json_object(line::AbstractString)
    p = _SJ(String(line), 1)
    _ws!(p)
    _peek(p) == '{' || throw(StrictJSONError("line does not start with '{'"))
    v = _value!(p)
    _ws!(p)
    p.i == ncodeunits(p.s) + 1 || throw(StrictJSONError("trailing characters after the object"))
    return v
end
