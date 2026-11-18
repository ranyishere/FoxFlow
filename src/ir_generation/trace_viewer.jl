"""
TraceViewer: renders the rule-fire trace of a FoxFlow simulation as a self-contained HTML report.

A simulation records its trace when the simulations section asks for it:

    simulations Microtubule {
        ...
        moveMicrotubules : Simulation := RunSimulation(...)
        trace(moveMicrotubules, "crossover", "zippering_case1")   # "*" traces every rule
    }

The generated model then writes `trace_<simulation>.jsonl` (one JSON record per rule fire:
matched nodes before/after, their 1-hop neighbourhood, edge changes, propensity and the
where-clause locals), and `ir.jl run` renders `trace_<simulation>.html` once the run ends.
`ir.jl trace` renders an existing trace again, with filters.

The compiler writes `trace_manifest.json` next to the generated sources; it lists, per traced
simulation, the rules, the .fflow file defining them and the files to read and write.

Dependency free (Printf and OrderedCollections only).
"""
module TraceViewer

using Printf
using OrderedCollections: OrderedDict

export write_manifest, render_manifest, render_trace_report, load_trace, trace_main,
       parse_json, to_json

const MANIFEST_NAME = "trace_manifest.json"

const PALETTE = ["#4c78a8", "#f58518", "#54a24b", "#e45756", "#72b7b2", "#b279a2",
                 "#ff9da6", "#9d755d", "#bab0ac", "#eeca3b", "#79706e", "#59a14f"]

# ---------------------------------------------------------------------------
# Minimal JSON reader / writer
# ---------------------------------------------------------------------------

function _skip_ws!(b::Vector{UInt8}, p::Ref{Int})
    while p[] <= length(b) && (b[p[]] == 0x20 || b[p[]] == 0x09 || b[p[]] == 0x0a || b[p[]] == 0x0d)
        p[] += 1
    end
end

function _json_string!(b::Vector{UInt8}, p::Ref{Int})
    p[] += 1  # opening quote
    out = UInt8[]
    while p[] <= length(b)
        c = b[p[]]
        if c == UInt8('"')
            p[] += 1
            return String(out)
        elseif c == UInt8('\\')
            p[] += 1
            e = Char(b[p[]])
            if e == 'n'; push!(out, UInt8('\n'))
            elseif e == 't'; push!(out, UInt8('\t'))
            elseif e == 'r'; push!(out, UInt8('\r'))
            elseif e == 'b'; push!(out, UInt8('\b'))
            elseif e == 'f'; push!(out, UInt8('\f'))
            elseif e == 'u'
                code = parse(UInt32, String(b[p[]+1:p[]+4]); base=16)
                append!(out, codeunits(string(Char(code))))
                p[] += 4
            else
                push!(out, UInt8(e))
            end
            p[] += 1
        else
            push!(out, c)
            p[] += 1
        end
    end
    error("JSON: unterminated string")
end

function _json_value!(b::Vector{UInt8}, p::Ref{Int})
    _skip_ws!(b, p)
    p[] <= length(b) || error("JSON: unexpected end of input")
    c = Char(b[p[]])
    if c == '{'
        obj = OrderedDict{String,Any}()
        p[] += 1
        _skip_ws!(b, p)
        if Char(b[p[]]) == '}'
            p[] += 1
            return obj
        end
        while true
            _skip_ws!(b, p)
            key = _json_string!(b, p)
            _skip_ws!(b, p)
            Char(b[p[]]) == ':' || error("JSON: expected ':' at byte $(p[])")
            p[] += 1
            obj[key] = _json_value!(b, p)
            _skip_ws!(b, p)
            d = Char(b[p[]])
            p[] += 1
            d == ',' && continue
            d == '}' && return obj
            error("JSON: expected ',' or '}' at byte $(p[] - 1)")
        end
    elseif c == '['
        arr = Any[]
        p[] += 1
        _skip_ws!(b, p)
        if Char(b[p[]]) == ']'
            p[] += 1
            return arr
        end
        while true
            push!(arr, _json_value!(b, p))
            _skip_ws!(b, p)
            d = Char(b[p[]])
            p[] += 1
            d == ',' && continue
            d == ']' && return arr
            error("JSON: expected ',' or ']' at byte $(p[] - 1)")
        end
    elseif c == '"'
        return _json_string!(b, p)
    elseif c == 't'
        p[] += 4
        return true
    elseif c == 'f'
        p[] += 5
        return false
    elseif c == 'n'
        p[] += 4
        return nothing
    else
        start = p[]
        while p[] <= length(b) && Char(b[p[]]) in "+-0123456789.eE"
            p[] += 1
        end
        text = String(b[start:p[]-1])
        isempty(text) && error("JSON: unexpected character '$c' at byte $start")
        x = tryparse(Int, text)
        return x === nothing ? parse(Float64, text) : x
    end
end

"""Parse a JSON document. Objects become `OrderedDict{String,Any}` (field order kept)."""
function parse_json(src::AbstractString)
    b = Vector{UInt8}(codeunits(src))
    return _json_value!(b, Ref(1))
end

function _json_escape(s::AbstractString)
    io = IOBuffer()
    for c in s
        if c == '"'; print(io, "\\\"")
        elseif c == '\\'; print(io, "\\\\")
        elseif c == '\n'; print(io, "\\n")
        elseif c == '\t'; print(io, "\\t")
        elseif c == '\r'; print(io, "\\r")
        elseif c < ' '; print(io, @sprintf("\\u%04x", Int(c)))
        else print(io, c)
        end
    end
    return String(take!(io))
end

_write_json(io::IO, ::Nothing) = print(io, "null")
_write_json(io::IO, x::Bool) = print(io, x ? "true" : "false")
_write_json(io::IO, x::Integer) = print(io, x)
function _write_json(io::IO, x::AbstractFloat)
    if isnan(x)
        print(io, "\"nan\"")
    elseif isinf(x)
        print(io, x > 0 ? "\"inf\"" : "\"-inf\"")
    else
        print(io, x)
    end
end
_write_json(io::IO, x::AbstractString) = print(io, '"', _json_escape(x), '"')
function _write_json(io::IO, x::AbstractVector)
    print(io, '[')
    for (i, v) in enumerate(x)
        i > 1 && print(io, ',')
        _write_json(io, v)
    end
    print(io, ']')
end
function _write_json(io::IO, x::AbstractDict)
    print(io, '{')
    first = true
    for (k, v) in x
        first || print(io, ',')
        first = false
        _write_json(io, string(k))
        print(io, ':')
        _write_json(io, v)
    end
    print(io, '}')
end

"""Serialize `x` as compact JSON."""
function to_json(x)
    io = IOBuffer()
    _write_json(io, x)
    return String(take!(io))
end

# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

"""Read a `*.jsonl` rule-fire trace (one record per non-empty line)."""
function load_trace(path::AbstractString)
    isfile(path) || error("trace file not found: $path")
    recs = OrderedDict{String,Any}[]
    for line in eachline(path)
        isempty(strip(line)) && continue
        push!(recs, parse_json(line))
    end
    return recs
end

"""Node types in variant order (== the `type` index stored in the trace records)."""
function type_names(types_header::AbstractString)
    isfile(types_header) || return String[]
    m = match(r"Spatial\w*<([^<>]*)>"s, read(types_header, String))
    m === nothing && return String[]
    return [String(strip(t)) for t in split(m.captures[1], ',') if !isempty(strip(t))]
end

type_name(names, t::Integer) = (1 <= t + 1 <= length(names)) ? names[t+1] : string(t)

_pos(n) = Float64[x for x in n["pos"]]
_esc(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;", "'" => "&#x27;")
_f1(x) = @sprintf("%.1f", x)

function _moved(a, b; eps=1e-9)
    return maximum(abs.(_pos(a) .- _pos(b))) > eps
end

function trace_moved(r)
    before = Dict{Int,Any}(n["id"] => n for n in r["before"])
    out = Int[]
    for n in r["after"]
        b = get(before, n["id"], nothing)
        b !== nothing && _moved(n, b) && push!(out, n["id"])
    end
    return out
end

# ---------------------------------------------------------------------------
# Rule definition (read back from the .fflow source)
# ---------------------------------------------------------------------------

struct RuleInfo
    rule::String
    file::String
    line::Int
    comment::String
    lhs_nodes::OrderedDict{String,String}
    lhs_edges::Vector{Tuple{String,String}}
    rhs_nodes::OrderedDict{String,String}
    rhs_edges::Vector{Tuple{String,String}}
    lhs_bind::Vector{NTuple{3,String}}
    rhs_bind::Vector{NTuple{3,String}}
    propensity::String
    where_block::String
    params::Vector{Tuple{String,String}}
end

function _parse_graph(text::AbstractString)
    nodes = OrderedDict{String,String}()
    edges = Set{Tuple{String,String}}()
    for line in split(text, '\n')
        prev = nothing
        pending = false
        for m in eachmatch(r"\(\s*(\w+)\s*:\s*(\w+)\s*\)|--", line)
            if m.match == "--"
                pending = prev !== nothing
                continue
            end
            al, ty = String(m.captures[1]), String(m.captures[2])
            haskey(nodes, al) || (nodes[al] = ty)
            if pending && prev != al
                push!(edges, prev < al ? (prev, al) : (al, prev))
            end
            prev = al
            pending = false
        end
    end
    return nodes, sort!(collect(edges))
end

function _bindings(text::AbstractString)
    out = NTuple{3,String}[]
    for blk in eachmatch(r"<<(.*?)>>"s, text)
        for m in eachmatch(r"(\w+)\s*:\s*(\w+)::(\w+)", blk.captures[1])
            push!(out, (String(m.captures[1]), String(m.captures[2]), String(m.captures[3])))
        end
    end
    return out
end

function _dedent(lines)
    nonblank = [l for l in lines if !isempty(strip(l))]
    ind = isempty(nonblank) ? 0 : minimum(length(l) - length(lstrip(l)) for l in nonblank)
    out = [rstrip(chop(l; head=min(ind, length(l)), tail=0)) for l in lines]
    while !isempty(out) && isempty(out[1]); popfirst!(out); end
    while !isempty(out) && isempty(out[end]); pop!(out); end
    return join(out, "\n")
end

_indent(l) = length(l) - length(lstrip(l))

"""
    load_rule(rules_file, rule) -> Union{RuleInfo,Nothing}

Read rule `rule` back from `rules_file` (graph patterns, bindings, propensity, where block and
the parameters it uses). Returns `nothing` when the rule is not defined there.
"""
function load_rule(rules_file::AbstractString, rule::AbstractString)
    isfile(rules_file) || return nothing
    lines = split(read(rules_file, String), '\n')
    pat = Regex("^(\\s*)\\Q" * rule * "\\E\\s*:=")
    start = findfirst(l -> occursin(pat, l), lines)
    start === nothing && return nothing
    ind = length(match(pat, lines[start]).captures[1])

    stop = start + 1
    while stop <= length(lines) && (isempty(strip(lines[stop])) || _indent(lines[stop]) > ind)
        stop += 1
    end
    raw = lines[start+1:stop-1]

    head = String[]
    j = start - 1
    while j >= 1 && startswith(strip(lines[j]), "#")
        pushfirst!(head, strip(lstrip(strip(lines[j]), '#')))
        j -= 1
    end

    code = replace(join(raw, "\n"), r"#[^\n]*" => "")
    arrow = findfirst("->", code)
    arrow === nothing && return nothing
    lhs_t = code[1:prevind(code, first(arrow))]
    rest = code[nextind(code, last(arrow)):end]

    prop = ""
    rhs_t = rest
    w = match(r"\bwith\s*\(", rest)
    if w !== nothing
        rhs_t = rest[1:prevind(rest, w.offset)]
        open_ix = w.offset + ncodeunits(w.match) - 1   # index of the "("
        depth = 0
        close_ix = lastindex(rest)
        for k in open_ix:lastindex(rest)
            isvalid(rest, k) || continue
            c = rest[k]
            depth += (c == '(') - (c == ')')
            if depth == 0
                close_ix = k
                break
            end
        end
        prop = _dedent(split(rest[nextind(rest, open_ix):prevind(rest, close_ix)], '\n'))
    end

    rawtxt = join(raw, "\n")
    where_block = ""
    wm = match(r"\bwhere\s*\{", rawtxt)
    if wm !== nothing
        body_start = wm.offset + ncodeunits(wm.match)
        last_brace = findlast('}', rawtxt)
        if last_brace !== nothing && body_start <= last_brace
            where_block = _dedent(split(rawtxt[body_start:prevind(rawtxt, last_brace)], '\n'))
        end
    end

    ln, le = _parse_graph(replace(lhs_t, r"<<.*?>>"s => ""))
    rn, re_ = _parse_graph(replace(rhs_t, r"<<.*?>>"s => ""))

    # Parameters: UPPER_CASE names defined in the sibling .fflow files.
    defs = Dict{String,String}()
    for f in sort(filter(f -> endswith(f, ".fflow"), readdir(dirname(abspath(rules_file)); join=true)))
        for m in eachmatch(r"^\s*([A-Za-z_]\w*)\s*:\s*\w+\s*:=\s*(.+)$"m, read(f, String))
            get!(defs, String(m.captures[1]), String(strip(m.captures[2])))
        end
    end
    scan = prop * replace(where_block, r"#[^\n]*" => "")
    used = sort!(unique([m.match for m in eachmatch(r"\b[A-Z][A-Z0-9_]{2,}\b", scan)]))
    params = [(String(u), defs[u]) for u in used if haskey(defs, u)]

    return RuleInfo(String(rule), String(rules_file), start, join(head, " "),
                    ln, le, rn, re_, _bindings(lhs_t), _bindings(rhs_t), prop, where_block, params)
end

"""system node id -> rule alias, using the engine's ordinal scheme
(LHS aliases 1..n in order of appearance, RHS-only aliases after them)."""
function rule_alias_for_record(r, info::RuleInfo)
    new = [x for x in keys(info.rhs_nodes) if !haskey(info.lhs_nodes, x)]
    order = vcat(collect(keys(info.lhs_nodes)), new)
    al = Dict{Int,String}()
    for (k, sid) in r["match"]
        ki = parse(Int, k)
        1 <= ki <= length(order) && (al[sid] = order[ki])
    end
    for (sid, name) in zip(sort(Int.(r["created"])), new)
        al[sid] = name
    end
    return al
end

# ---------------------------------------------------------------------------
# SVG rendering
# ---------------------------------------------------------------------------

function trace_svg(r, names; size=440, alias=Dict{Int,String}())
    before = Dict{Int,Any}(n["id"] => n for n in r["before"])
    after = Dict{Int,Any}(n["id"] => n for n in r["after"])
    ctx = Dict{Int,Any}(n["id"] => n for n in r["context"])
    created = Set{Int}(Int.(r["created"]))
    destroyed = Set{Int}(Int.(r["destroyed"]))

    allpts = [_pos(n) for d in (before, after, ctx) for n in values(d)]
    xs = [p[1] for p in allpts]
    ys = [p[2] for p in allpts]
    x0, x1, y0, y1 = minimum(xs), maximum(xs), minimum(ys), maximum(ys)
    span = max(x1 - x0, y1 - y0, 1e-9)
    pad = 0.12 * span
    x0 -= pad
    y0 -= pad
    span += 2 * pad
    sc = size / span

    P(p) = ((p[1] - x0) * sc, size - (p[2] - y0) * sc)

    function panel(nodes, edges, other, which)
        s = String["<svg width=\"$size\" height=\"$size\" viewBox=\"0 0 $size $size\" class=\"p\"><text x=\"6\" y=\"16\">$which</text>"]
        pos = merge(ctx, nodes)
        for (a_, b_) in edges
            if haskey(pos, a_) && haskey(pos, b_)
                (xa, ya), (xb, yb) = P(_pos(pos[a_])), P(_pos(pos[b_]))
                push!(s, "<line x1=\"$(_f1(xa))\" y1=\"$(_f1(ya))\" x2=\"$(_f1(xb))\" y2=\"$(_f1(yb))\" stroke=\"#999\" stroke-width=\"1.5\"/>")
            end
        end
        for n in sort!(collect(values(ctx)); by=n -> n["id"])
            x, y = P(_pos(n))
            push!(s, "<circle cx=\"$(_f1(x))\" cy=\"$(_f1(y))\" r=\"4\" fill=\"#ccc\"/><text x=\"$(_f1(x + 5))\" y=\"$(_f1(y - 5))\" class=\"c\">$(n["id"])</text>")
        end
        for n in sort!(collect(values(nodes)); by=n -> n["id"])
            x, y = P(_pos(n))
            col = PALETTE[mod(n["type"], length(PALETTE))+1]
            o = get(other, n["id"], nothing)
            ring = ""
            if which == "after" && n["id"] in created
                ring = " stroke=\"#00a000\" stroke-width=\"3\""
            elseif which == "before" && n["id"] in destroyed
                ring = " stroke=\"#d00000\" stroke-width=\"3\""
            end
            if o !== nothing && which == "after" && _moved(n, o)
                gx, gy = P(_pos(o))
                push!(s, "<line x1=\"$(_f1(gx))\" y1=\"$(_f1(gy))\" x2=\"$(_f1(x))\" y2=\"$(_f1(y))\" stroke=\"#888\" stroke-dasharray=\"3,3\"/>" *
                         "<circle cx=\"$(_f1(gx))\" cy=\"$(_f1(gy))\" r=\"5\" fill=\"none\" stroke=\"#888\" stroke-dasharray=\"2,2\"/>")
            end
            d = get(n["data"], "Direction", nothing)
            if d isa AbstractVector && length(d) >= 2
                push!(s, "<line x1=\"$(_f1(x))\" y1=\"$(_f1(y))\" x2=\"$(_f1(x + 24 * d[1]))\" y2=\"$(_f1(y - 24 * d[2]))\" stroke=\"$col\" stroke-width=\"2\"/>")
            end
            push!(s, "<circle cx=\"$(_f1(x))\" cy=\"$(_f1(y))\" r=\"6\" fill=\"$col\"$ring><title>$(_esc(to_json(n["data"])))</title></circle>")
            tn = type_name(names, n["type"])
            if haskey(alias, n["id"])
                push!(s, "<text x=\"$(_f1(x + 8))\" y=\"$(_f1(y + 12))\" class=\"n\"><tspan font-weight=\"bold\">$(alias[n["id"]])</tspan> $(n["id"]) $tn</text>")
            else
                push!(s, "<text x=\"$(_f1(x + 8))\" y=\"$(_f1(y + 12))\" class=\"n\">$(n["id"]) $tn</text>")
            end
        end
        push!(s, "</svg>")
        return join(s)
    end

    eb = [(Int(e[1]), Int(e[2])) for e in r["edges_before"]]
    ea = [(Int(e[1]), Int(e[2])) for e in r["edges_after"]]
    return panel(before, eb, after, "before") * panel(after, ea, before, "after")
end

function trace_attr_changes(r, names)
    before = Dict{Int,Any}(n["id"] => n for n in r["before"])
    rows = Tuple{Any,String,String}[]
    for n in r["after"]
        b = get(before, n["id"], nothing)
        tn = type_name(names, n["type"])
        if b === nothing
            push!(rows, (n["id"], "created as " * tn, ""))
            continue
        end
        ch = String[]
        if b["type"] != n["type"]
            push!(ch, "type $(type_name(names, b["type"])) -> $tn")
        end
        if _moved(n, b)
            bp, np = _pos(b), _pos(n)
            push!(ch, @sprintf("pos (%.4f, %.4f, %.4f) -> (%.4f, %.4f, %.4f)",
                               bp[1], bp[2], bp[3], np[1], np[2], np[3]))
        end
        for (k, v) in n["data"]
            bv = get(b["data"], k, nothing)
            bv != v && push!(ch, "$k $(to_json(bv)) -> $(to_json(v))")
        end
        isempty(ch) || push!(rows, (n["id"], "changed", join(ch, "; ")))
    end
    for i in r["destroyed"]
        push!(rows, (i, "destroyed", ""))
    end
    return rows
end

# --- rule graph (LHS/RHS) ---

_type_color(ty, names) = (ix = findfirst(==(ty), names); ix === nothing ? "#999" : PALETTE[mod(ix - 1, length(PALETTE))+1])

function _layout(nodes::Vector{String}, edges, w, h; margin=46)
    n = length(nodes)
    n == 0 && return Dict{String,Tuple{Float64,Float64}}()
    n == 1 && return Dict(nodes[1] => (w / 2, h / 2))
    p = Dict(a => [cos(2pi * (i - 1) / n) + 0.01 * (i - 1), sin(2pi * (i - 1) / n)] for (i, a) in enumerate(nodes))
    for it in 0:299
        t = 0.12 * (1 - it / 300)
        disp = Dict(a => [0.0, 0.0] for a in nodes)
        for (i, a) in enumerate(nodes)
            for b in nodes[i+1:end]
                dx, dy = p[a][1] - p[b][1], p[a][2] - p[b][2]
                d = hypot(dx, dy)
                d == 0 && (d = 1e-3)
                f = 1 / (d * d)
                disp[a][1] += dx / d * f
                disp[a][2] += dy / d * f
                disp[b][1] -= dx / d * f
                disp[b][2] -= dy / d * f
            end
        end
        for (a, b) in edges
            (haskey(p, a) && haskey(p, b)) || continue
            dx, dy = p[a][1] - p[b][1], p[a][2] - p[b][2]
            d = hypot(dx, dy)
            d == 0 && (d = 1e-3)
            f = d * d
            disp[a][1] -= dx / d * f
            disp[a][2] -= dy / d * f
            disp[b][1] += dx / d * f
            disp[b][2] += dy / d * f
        end
        for a in nodes
            m = hypot(disp[a][1], disp[a][2])
            m == 0 && (m = 1e-9)
            s = min(m, t) / m
            p[a][1] += disp[a][1] * s
            p[a][2] += disp[a][2] * s
        end
    end
    xs = [v[1] for v in values(p)]
    ys = [v[2] for v in values(p)]
    x0, x1, y0, y1 = minimum(xs), maximum(xs), minimum(ys), maximum(ys)
    sx = (w - 2margin) / max(x1 - x0, 1e-9)
    sy = (h - 2margin) / max(y1 - y0, 1e-9)
    sc = min(sx, sy)
    ox = (w - sc * (x1 - x0)) / 2
    oy = (h - sc * (y1 - y0)) / 2
    return Dict(a => (ox + (p[a][1] - x0) * sc, oy + (p[a][2] - y0) * sc) for a in nodes)
end

function rule_graph_svg(info::RuleInfo, names; w=400, h=290)
    union_nodes = merge(info.lhs_nodes, info.rhs_nodes)
    all_edges = sort!(collect(union(Set(info.lhs_edges), Set(info.rhs_edges))))
    pos = _layout(collect(keys(union_nodes)), all_edges, w, h)

    function panel(title, nodes, edges_, other_nodes, other_edges, side)
        s = String["<svg width=\"$w\" height=\"$h\" class=\"p\"><text x=\"6\" y=\"16\" font-weight=\"bold\">$title</text>"]
        for (a_, b_) in edges_
            col, dash = "#999", ""
            if !((a_, b_) in other_edges)
                col = side == "lhs" ? "#d00000" : "#00a000"
                dash = " stroke-dasharray=\"5,3\""
            end
            push!(s, "<line x1=\"$(_f1(pos[a_][1]))\" y1=\"$(_f1(pos[a_][2]))\" x2=\"$(_f1(pos[b_][1]))\" y2=\"$(_f1(pos[b_][2]))\" stroke=\"$col\" stroke-width=\"2\"$dash/>")
        end
        for (al, ty) in nodes
            x, y = pos[al]
            stroke, sw, dash = "#333", 1, ""
            if !haskey(other_nodes, al)
                stroke, sw = side == "lhs" ? ("#d00000", 3) : ("#00a000", 3)
            elseif other_nodes[al] != ty
                stroke, sw, dash = "#ff8c00", 3, " stroke-dasharray=\"3,2\""
            end
            push!(s, "<circle cx=\"$(_f1(x))\" cy=\"$(_f1(y))\" r=\"11\" fill=\"$(_type_color(ty, names))\" stroke=\"$stroke\" stroke-width=\"$sw\"$dash/>")
            push!(s, "<text x=\"$(_f1(x))\" y=\"$(_f1(y - 16))\" class=\"n\" text-anchor=\"middle\"><tspan font-weight=\"bold\">$al</tspan></text>" *
                     "<text x=\"$(_f1(x))\" y=\"$(_f1(y + 24))\" class=\"c\" text-anchor=\"middle\">$ty</text>")
        end
        push!(s, "</svg>")
        return join(s)
    end

    le, re_ = Set(info.lhs_edges), Set(info.rhs_edges)
    return panel("LHS (matched)", info.lhs_nodes, info.lhs_edges, info.rhs_nodes, re_, "lhs") *
           panel("RHS (result)", info.rhs_nodes, info.rhs_edges, info.lhs_nodes, le, "rhs")
end

function rule_section_html(info::RuleInfo, names; base_dir="")
    order = vcat(collect(keys(info.lhs_nodes)), [x for x in keys(info.rhs_nodes) if !haskey(info.lhs_nodes, x)])
    ids = join(["$i=$a" for (i, a) in enumerate(order)], ", ")
    shown_file = isempty(base_dir) ? info.file : relpath(info.file, base_dir)
    o = String["<section class=\"rule\"><h2>Rule <code>$(_esc(info.rule))</code></h2>",
               "<div class=\"src\">$(_esc(shown_file)):$(info.line)</div>"]
    isempty(info.comment) || push!(o, "<p>$(_esc(info.comment))</p>")
    push!(o, "<div>" * rule_graph_svg(info, names) * "</div>")
    push!(o, "<p class=\"c2\">Pattern ids: $(_esc(ids)). Green = created/new edge, red = destroyed/removed edge, orange dashed = retyped in place.</p>")
    rows = max(length(info.lhs_bind), length(info.rhs_bind))
    if rows > 0
        push!(o, "<table><tr><th colspan=2>LHS bindings<th colspan=2>RHS bindings</tr>")
        cell(b) = b === nothing ? "<td><td>" : "<td>$(_esc(b[1]))<td>$(_esc(b[2]))::$(_esc(b[3]))"
        for i in 1:rows
            l = i <= length(info.lhs_bind) ? info.lhs_bind[i] : nothing
            r = i <= length(info.rhs_bind) ? info.rhs_bind[i] : nothing
            push!(o, "<tr>$(cell(l))$(cell(r))</tr>")
        end
        push!(o, "</table>")
    end
    if !isempty(info.params)
        push!(o, "<table><tr><th>parameter<th>definition</tr>" *
                 join(["<tr><td>$(_esc(k))<td>$(_esc(v))" for (k, v) in info.params]) * "</table>")
    end
    isempty(info.propensity) || push!(o, "<h4>propensity (with)</h4><pre>$(_esc(info.propensity))</pre>")
    isempty(info.where_block) || push!(o, "<details><summary>where block</summary><pre>$(_esc(info.where_block))</pre></details>")
    push!(o, "</section>")
    return join(o)
end

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

"""Indices (1-based) of the records matching the filters."""
function select_records(recs; rule=nothing, sim=nothing, index=nothing, node=nothing)
    idx = collect(1:length(recs))
    rule === nothing || filter!(i -> recs[i]["rule"] == rule, idx)
    sim === nothing || filter!(i -> recs[i]["sim_step"] == sim, idx)
    index === nothing || filter!(i -> i - 1 == index, idx)
    node === nothing || filter!(i -> any(n["id"] == node for n in vcat(recs[i]["before"], recs[i]["after"])), idx)
    return idx
end

const REPORT_CSS = """
body{font-family:sans-serif;margin:16px}.p{border:1px solid #ccc;margin:2px;background:#fff}
.n{font-size:10px}.c{font-size:9px;fill:#888}table{border-collapse:collapse;font-size:12px}
td,th{border:1px solid #ddd;padding:2px 6px;text-align:left}h3{margin-top:28px}
pre{background:#f6f6f6;padding:8px;overflow:auto;font-size:12px}.src{color:#888;font-size:12px}
.c2{font-size:12px;color:#666}.rule{border:1px solid #bbb;padding:8px 14px;margin-bottom:16px;
background:#fcfcfc}.lg span{display:inline-block;margin-right:12px;font-size:12px}
.lg i{display:inline-block;width:10px;height:10px;border-radius:5px;margin-right:4px}"""

"""
    render_trace_report(recs, names, rules_file, out; title, max, rule, sim, index, node, base_dir)

Write an HTML report of the selected records: one section per traced rule (its graph pattern,
bindings, propensity and where block, read from `rules_file`) followed by before/after panels
for each fire. Returns the number of records shown.
"""
function render_trace_report(recs, names, rules_file::Union{AbstractString,Nothing}, out::AbstractString;
                             title::AbstractString="FoxFlow rule trace", max::Integer=40,
                             rule=nothing, sim=nothing, index=nothing, node=nothing, base_dir="")
    idx = select_records(recs; rule=rule, sim=sim, index=index, node=node)
    shown = idx[1:min(length(idx), max)]

    parts = String["<!doctype html><meta charset=utf-8><title>$(_esc(title))</title><style>$REPORT_CSS</style>",
                   "<h1>$(_esc(title))</h1>",
                   "<p class=\"c2\">$(length(shown)) of $(length(idx)) selected fires shown ($(length(recs)) recorded)" *
                   (length(idx) > length(shown) ? "; raise the limit with --max" : "") * "</p>"]
    if !isempty(names)
        push!(parts, "<div class=\"lg\">" * join(["<span><i style=\"background:$(PALETTE[mod(k - 1, length(PALETTE))+1])\"></i>$(_esc(n))</span>"
                                                 for (k, n) in enumerate(names)]) * "</div>")
    end

    infos = Dict{String,Union{RuleInfo,Nothing}}()
    for rname in unique(recs[i]["rule"] for i in shown)
        info = rules_file === nothing ? nothing : load_rule(rules_file, rname)
        infos[rname] = info
        if info === nothing
            println("[trace] note: rule '$rname' not found in $(rules_file === nothing ? "any .fflow file" : rules_file)")
        else
            push!(parts, rule_section_html(info, names; base_dir=base_dir))
        end
    end

    for i in shown
        r = recs[i]
        info = infos[r["rule"]]
        alias = info === nothing ? Dict{Int,String}() : rule_alias_for_record(r, info)
        push!(parts, "<h3>#$(i - 1) $(_esc(r["rule"])) &mdash; sim step $(r["sim_step"]), ssa step $(r["ssa_step"]), propensity $(_esc(r["propensity"]))</h3>")
        push!(parts, "<div>" * trace_svg(r, names; alias=alias) * "</div>")
        push!(parts, "<p>match (pattern key &rarr; node id): " *
                     join(["$k&rarr;$v" for (k, v) in r["match"]], ", ") * "</p>")
        rows = trace_attr_changes(r, names)
        if !isempty(rows)
            push!(parts, "<table><tr><th>node<th>what<th>detail</tr>" *
                         join(["<tr><td>$a<td>$b<td>$(_esc(c))" for (a, b, c) in rows]) * "</table>")
        end
        if !isempty(r["locals"])
            push!(parts, "<table><tr><th>where-local<th>value</tr>" *
                         join(["<tr><td>$k<td>$(_esc(to_json(v)))" for (k, v) in r["locals"]]) * "</table>")
        end
    end

    mkpath(dirname(abspath(out)))
    write(out, join(parts))
    return length(shown)
end

function print_records(recs, idx; max=40)
    @printf("%5s %-22s %6s %6s %10s %7s %9s %5s\n", "idx", "rule", "sim", "ssa", "propensity", "created", "destroyed", "moved")
    for i in idx[1:min(length(idx), max)]
        r = recs[i]
        p = r["propensity"] isa Real ? Float64(r["propensity"]) : NaN
        @printf("%5d %-22s %6d %6d %10.4g %7d %9d %5d\n", i - 1, r["rule"], r["sim_step"], r["ssa_step"], p,
                length(r["created"]), length(r["destroyed"]), length(trace_moved(r)))
    end
    println("$(min(length(idx), max)) of $(length(idx)) records", length(idx) <= max ? "" : " (use --max)")
end

# ---------------------------------------------------------------------------
# Manifest (written by the compiler, read back after the run)
# ---------------------------------------------------------------------------

"""
    write_manifest(output_dir; input_dir, types_header, traces)

Record the `trace(...)` statements of the simulations section. Each entry of `traces` is a
`Dict` with `simulation`, `rules`, `rules_file` and `trace_file`. Removes a stale manifest when
there are no traces.
"""
function write_manifest(output_dir::AbstractString; input_dir::AbstractString,
                        types_header::AbstractString, traces::AbstractVector)
    path = joinpath(output_dir, MANIFEST_NAME)
    if isempty(traces)
        isfile(path) && rm(path)
        return nothing
    end
    doc = OrderedDict{String,Any}(
        "input_dir" => abspath(input_dir),
        "types_header" => abspath(types_header),
        "traces" => traces,
    )
    write(path, to_json(doc) * "\n")
    return path
end

"""
    render_manifest(output_dir, run_dir; kwargs...) -> Vector{String}

Render the HTML report of every traced simulation listed in `<output_dir>/trace_manifest.json`
from the `trace_<simulation>.jsonl` the model wrote into `run_dir`. Keywords `simulation`,
`rule`, `sim`, `index`, `node`, `max` filter the report. Returns the HTML paths written.
"""
function render_manifest(output_dir::AbstractString, run_dir::AbstractString;
                         simulation=nothing, rule=nothing, sim=nothing, index=nothing,
                         node=nothing, max::Integer=40, output=nothing)
    mpath = joinpath(output_dir, MANIFEST_NAME)
    isfile(mpath) || return String[]
    manifest = parse_json(read(mpath, String))
    names = type_names(manifest["types_header"])
    written = String[]
    for t in manifest["traces"]
        simulation === nothing || t["simulation"] == simulation || continue
        jsonl = joinpath(run_dir, t["trace_file"])
        if !isfile(jsonl)
            println("[trace] $(t["simulation"]): no trace file at $jsonl (did the run reach this simulation?)")
            continue
        end
        recs = load_trace(jsonl)
        if isempty(recs)
            println("[trace] $(t["simulation"]): no fires of $(join(t["rules"], ", ")) were recorded")
            continue
        end
        html_path = output !== nothing ? String(output) :
                    joinpath(run_dir, "trace_$(t["simulation"]).html")
        title = "Rule trace: $(t["simulation"]) ($(join(t["rules"], ", ")))"
        n = render_trace_report(recs, names, t["rules_file"], html_path;
                                title=title, max=max, rule=rule, sim=sim, index=index,
                                node=node, base_dir=manifest["input_dir"])
        println("[trace] wrote $html_path ($n of $(length(recs)) records)")
        push!(written, html_path)
    end
    return written
end

# ---------------------------------------------------------------------------
# CLI:  ir.jl trace <input_dir> [output_dir] [options]
# ---------------------------------------------------------------------------

function print_trace_usage()
    println("""
    Render rule-fire traces recorded by a model that ran with trace(sim, "rule", ...).

    Usage:
      julia ir.jl trace <input_dir> [output_dir] [options]

    Options:
      --build <dir>     Directory the model ran in (default: <output_dir>/build).
      --sim-name <n>    Only this simulation (the name used in trace(<n>, ...)).
      --rule <name>     Only fires of this rule.
      --step <n>        Only fires at this simulation step.
      --index <n>       Only this record index.
      --node <id>       Only fires involving this node id.
      --max <n>         Show at most n fires (default 40).
      --list            Print the selected records instead of writing HTML.
      -o <file.html>    Output path (default: <build>/trace_<simulation>.html).
    """)
end

function trace_main(args)
    if any(a -> a in ("-h", "--help", "help"), args)
        print_trace_usage()
        return 0
    end
    build = nothing; simulation = nothing; rule = nothing; step = nothing; index = nothing
    node = nothing; max = 40; list = false; output = nothing
    positionals = String[]
    i = 1
    intopt(name) = (i + 1 <= length(args) && tryparse(Int, args[i+1]) !== nothing) ?
                   parse(Int, args[i+1]) : error("$name expects an integer")
    while i <= length(args)
        a = args[i]
        if a == "--build"; build = args[i+1]; i += 2
        elseif a == "--sim-name"; simulation = args[i+1]; i += 2
        elseif a == "--rule"; rule = args[i+1]; i += 2
        elseif a == "--step"; step = intopt(a); i += 2
        elseif a == "--index"; index = intopt(a); i += 2
        elseif a == "--node"; node = intopt(a); i += 2
        elseif a == "--max"; max = intopt(a); i += 2
        elseif a == "--list"; list = true; i += 1
        elseif a == "-o" || a == "--output"; output = args[i+1]; i += 2
        elseif startswith(a, "-")
            println("Unknown option: $a")
            print_trace_usage()
            return 1
        else
            push!(positionals, a); i += 1
        end
    end
    isempty(positionals) && (println("Error: missing <input_dir>."); print_trace_usage(); return 1)
    input_dir = positionals[1]
    output_dir = length(positionals) >= 2 ? positionals[2] : joinpath(input_dir, "generated")
    build === nothing && (build = joinpath(output_dir, "build"))

    if !isfile(joinpath(output_dir, MANIFEST_NAME))
        println("Error: no $MANIFEST_NAME in $output_dir.")
        println("Add trace(<simulation>, \"<rule>\") to simulation.fflow, then run `gen`/`run` first.")
        return 1
    end

    if list
        manifest = parse_json(read(joinpath(output_dir, MANIFEST_NAME), String))
        for t in manifest["traces"]
            simulation === nothing || t["simulation"] == simulation || continue
            jsonl = joinpath(build, t["trace_file"])
            isfile(jsonl) || (println("no trace file at $jsonl"); continue)
            recs = load_trace(jsonl)
            println("== $(t["simulation"]) ($(jsonl))")
            print_records(recs, select_records(recs; rule=rule, sim=step, index=index, node=node); max=max)
        end
        return 0
    end

    written = render_manifest(output_dir, build; simulation=simulation, rule=rule, sim=step,
                              index=index, node=node, max=max, output=output)
    return isempty(written) ? 1 : 0
end

end # module TraceViewer
