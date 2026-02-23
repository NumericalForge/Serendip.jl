# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

const _typeset_symbols = Dict(
    "alpha" => "α", "beta" => "β", "gamma" => "γ", "delta" => "δ", "epsilon" => "ε", "varepsilon" => "ε",
    "zeta" => "ζ", "eta" => "η", "theta" => "θ", "iota" => "ι", "kappa" => "κ", "lambda" => "λ",
    "mu" => "μ", "nu" => "ν", "xi" => "ξ", "omicron" => "ο", "pi" => "π", "rho" => "ρ",
    "sigma" => "σ", "tau" => "τ", "upsilon" => "υ", "phi" => "φ", "chi" => "χ", "psi" => "ψ", "omega" => "ω",
    "Gamma" => "Γ", "Delta" => "Δ", "Theta" => "Θ", "Lambda" => "Λ", "Xi" => "Ξ", "Pi" => "Π",
    "Sigma" => "Σ", "Upsilon" => "Υ", "Phi" => "Φ", "Psi" => "Ψ", "Omega" => "Ω",
    "times" => "×", "nabla" => "∇", "partial" => "∂", "cdot" => "⋅",
)

const _typeset_functions = Set([
    "sin", "cos", "tan", "asin", "acos", "atan",
    "sinh", "cosh", "tanh",
    "log", "ln", "exp", "sqrt",
    "min", "max",
])

const _typeset_nonitalic_symbols = Dict(
    "times" => true,
    "nabla" => true,
    "partial" => true,
    "cdot" => true,
)
const _typeset_binary_ops = Set(["+", "−", "×", "⋅", "="])

abstract type TypesetNode end

struct TSAtom <: TypesetNode
    text::String
    italic::Bool
    bold::Bool
end

struct TSGroup <: TypesetNode
    items::Vector{TypesetNode}
end

struct TSScripts <: TypesetNode
    base::TypesetNode
    sub::Union{Nothing, TypesetNode}
    sup::Union{Nothing, TypesetNode}
end

struct TSParenGroup <: TypesetNode
    items::Vector{TypesetNode}
end

struct TSSqrt <: TypesetNode
    body::TypesetNode
end

struct TSFraction <: TypesetNode
    num::TypesetNode
    den::TypesetNode
end

struct GlyphElem
    text::String
    x::Float64
    y::Float64
    width::Float64
    height::Float64
    fontsize::Float64
    italic::Bool
    bold::Bool
    xscale::Float64
end

struct RuleElem
    x::Float64
    y::Float64
    width::Float64
    thickness::Float64
end

struct TypesetLayout
    glyphs::Vector{GlyphElem}
    rules::Vector{RuleElem}
    width::Float64
    height::Float64
    min_y::Float64
    max_y::Float64
end

const _paren_xscale_regular = 0.85
const _paren_xscale_min = 0.6
const _paren_outer_pad_factor = 0.06

@inline function _paren_xscale_for_height(body_h::Float64, fontsize::Float64)
    ratio = body_h / max(fontsize, 1e-9)
    scale = _paren_xscale_regular - 0.35 * (ratio - 1.0)
    return clamp(scale, _paren_xscale_min, _paren_xscale_regular)
end

@inline _next(s::AbstractString, i::Int) = nextind(s, i)
@inline _isletter(c::Char) = isletter(c)
@inline _isdigit(c::Char) = isdigit(c)
@inline function _skip_spaces(s::AbstractString, i::Int)
    while i <= lastindex(s) && isspace(s[i])
        i = _next(s, i)
    end
    return i
end

function _atom_from_char(c::Char; in_math::Bool=false)
    if in_math && _isletter(c)
        return TSAtom(string(c), true, false)
    end
    return TSAtom(string(c), false, false)
end

function _read_identifier(s::AbstractString, i::Int)
    j = i
    while j <= lastindex(s) && (_isletter(s[j]) || _isdigit(s[j]))
        j = _next(s, j)
    end
    return s[i:prevind(s, j)], j
end

function _identifier_nodes(name::AbstractString)
    if haskey(_typeset_symbols, name)
        italic = !haskey(_typeset_nonitalic_symbols, name)
        return TypesetNode[TSAtom(_typeset_symbols[name], italic, false)]
    elseif name in _typeset_functions
        return TypesetNode[TSAtom(name, false, false)]
    else
        return [TSAtom(string(c), true, false) for c in name]
    end
end

function _parse_script_arg(s::AbstractString, i::Int)
    i > lastindex(s) && return TSAtom("", false, false), i

    c = s[i]
    if c == '{'
        nodes, j = _parse_math_seq(s, _next(s, i); stopchars=Set(['}']))
        return TSGroup(nodes), j <= lastindex(s) && s[j] == '}' ? _next(s, j) : j
    elseif c == '('
        nodes, j = _parse_math_seq(s, _next(s, i); stopchars=Set([')']))
        return TSGroup(nodes), j <= lastindex(s) && s[j] == ')' ? _next(s, j) : j
    elseif _isletter(c)
        name, j = _read_identifier(s, i)
        if name == "sqrt"
            sqrt_node, k = _parse_sqrt_call(s, j)
            sqrt_node !== nothing && return sqrt_node, k
        end
        nodes = _identifier_nodes(name)
        return length(nodes) == 1 ? nodes[1] : TSGroup(nodes), j
    else
        return _atom_from_char(c, in_math=true), _next(s, i)
    end
end

@inline function _nodes_to_node(nodes::AbstractVector{<:TypesetNode})
    if isempty(nodes)
        return TSAtom("", false, false)
    end
    return length(nodes) == 1 ? nodes[1] : TSGroup(nodes)
end

function _parse_fraction_call(s::AbstractString, i::Int)
    i = _skip_spaces(s, i)
    if i > lastindex(s) || s[i] != '('
        return nothing, i
    end

    start = _next(s, i)
    pos = start
    depth = 1
    comma_pos = 0
    while pos <= lastindex(s)
        c = s[pos]
        if c == '('
            depth += 1
        elseif c == ')'
            depth -= 1
            depth == 0 && break
        elseif c == ',' && depth == 1 && comma_pos == 0
            comma_pos = pos
        end
        pos = _next(s, pos)
    end

    if depth != 0 || comma_pos == 0
        return nothing, i
    end

    num_str = strip(s[start:prevind(s, comma_pos)])
    den_str = strip(s[_next(s, comma_pos):prevind(s, pos)])

    num_nodes, _ = isempty(num_str) ? (TypesetNode[], 1) : _parse_math_seq(num_str, firstindex(num_str))
    den_nodes, _ = isempty(den_str) ? (TypesetNode[], 1) : _parse_math_seq(den_str, firstindex(den_str))

    frac = TSFraction(_nodes_to_node(num_nodes), _nodes_to_node(den_nodes))
    return frac, _next(s, pos)
end

function _parse_sqrt_call(s::AbstractString, i::Int)
    i = _skip_spaces(s, i)
    if i > lastindex(s) || s[i] != '('
        return nothing, i
    end

    start = _next(s, i)
    pos = start
    depth = 1
    while pos <= lastindex(s)
        c = s[pos]
        if c == '('
            depth += 1
        elseif c == ')'
            depth -= 1
            depth == 0 && break
        end
        pos = _next(s, pos)
    end

    depth != 0 && return nothing, i

    body_str = strip(s[start:prevind(s, pos)])
    body_nodes, _ = isempty(body_str) ? (TypesetNode[], 1) : _parse_math_seq(body_str, firstindex(body_str))
    return TSSqrt(_nodes_to_node(body_nodes)), _next(s, pos)
end

function _set_bold(node::TypesetNode)::TypesetNode
    if node isa TSAtom
        a = node::TSAtom
        return TSAtom(a.text, a.italic, true)
    elseif node isa TSGroup
        g = node::TSGroup
        return TSGroup([_set_bold(it) for it in g.items])
    elseif node isa TSParenGroup
        g = node::TSParenGroup
        return TSParenGroup([_set_bold(it) for it in g.items])
    elseif node isa TSScripts
        s = node::TSScripts
        sub = s.sub === nothing ? nothing : _set_bold(s.sub)
        sup = s.sup === nothing ? nothing : _set_bold(s.sup)
        return TSScripts(_set_bold(s.base), sub, sup)
    elseif node isa TSSqrt
        sq = node::TSSqrt
        return TSSqrt(_set_bold(sq.body))
    else
        f = node::TSFraction
        return TSFraction(_set_bold(f.num), _set_bold(f.den))
    end
end

function _parse_bold_call(s::AbstractString, i::Int)
    i = _skip_spaces(s, i)
    if i > lastindex(s) || s[i] != '('
        return nothing, i
    end

    start = _next(s, i)
    pos = start
    depth = 1
    while pos <= lastindex(s)
        c = s[pos]
        if c == '('
            depth += 1
        elseif c == ')'
            depth -= 1
            depth == 0 && break
        end
        pos = _next(s, pos)
    end

    depth != 0 && return nothing, i

    body_str = strip(s[start:prevind(s, pos)])
    body_nodes, _ = isempty(body_str) ? (TypesetNode[], 1) : _parse_math_seq(body_str, firstindex(body_str))
    bold_nodes = [_set_bold(n) for n in body_nodes]
    return _nodes_to_node(bold_nodes), _next(s, pos)
end

function _attach_script(base::TypesetNode, op::Char, arg::TypesetNode)
    if base isa TSScripts
        s = base::TSScripts
        if op == '_' && s.sub === nothing
            return TSScripts(s.base, arg, s.sup)
        elseif op == '^' && s.sup === nothing
            return TSScripts(s.base, s.sub, arg)
        else
            return TSScripts(base, op == '_' ? arg : nothing, op == '^' ? arg : nothing)
        end
    end
    return TSScripts(base, op == '_' ? arg : nothing, op == '^' ? arg : nothing)
end

function _bind_script!(nodes::Vector{TypesetNode}, op::Char, arg::TypesetNode)
    if isempty(nodes)
        push!(nodes, TSAtom(string(op), false, false))
        return
    end
    push!(nodes, _attach_script(pop!(nodes), op, arg))
end

function _parse_math_seq(s::AbstractString, i::Int; stopchars::Set{Char}=Set{Char}())
    nodes = TypesetNode[]

    while i <= lastindex(s)
        c = s[i]
        c in stopchars && return nodes, i

        if isspace(c)
            i = _next(s, i)
        elseif _isletter(c)
            name, i = _read_identifier(s, i)
            if name == "frac"
                frac, j = _parse_fraction_call(s, i)
                if frac === nothing
                    append!(nodes, _identifier_nodes(name))
                else
                    push!(nodes, frac)
                    i = j
                end
            elseif name == "sqrt"
                sqrt_node, j = _parse_sqrt_call(s, i)
                if sqrt_node === nothing
                    append!(nodes, _identifier_nodes(name))
                else
                    push!(nodes, sqrt_node)
                    i = j
                end
            elseif name == "bold"
                bold_node, j = _parse_bold_call(s, i)
                if bold_node === nothing
                    append!(nodes, _identifier_nodes(name))
                else
                    push!(nodes, bold_node)
                    i = j
                end
            else
                append!(nodes, _identifier_nodes(name))
            end
        elseif c == '{'
            inner, j = _parse_math_seq(s, _next(s, i); stopchars=Set(['}']))
            push!(nodes, TSGroup(inner))
            i = j <= lastindex(s) && s[j] == '}' ? _next(s, j) : j
        elseif c == '('
            inner, j = _parse_math_seq(s, _next(s, i); stopchars=Set([')']))
            if j <= lastindex(s) && s[j] == ')'
                push!(nodes, TSParenGroup(inner))
                i = _next(s, j)
            else
                push!(nodes, TSAtom("(", false, false))
                i = _next(s, i)
            end
        elseif c == '/'
            if isempty(nodes)
                push!(nodes, TSAtom("/", false, false))
                i = _next(s, i)
            else
                num = pop!(nodes)
                den, j = _parse_script_arg(s, _skip_spaces(s, _next(s, i)))
                while j <= lastindex(s) && (s[j] == '_' || s[j] == '^')
                    op = s[j]
                    arg, j = _parse_script_arg(s, _next(s, j))
                    den = _attach_script(den, op, arg)
                end
                push!(nodes, TSFraction(num, den))
                i = j
            end
        elseif c == '_'
            arg, j = _parse_script_arg(s, _next(s, i))
            _bind_script!(nodes, '_', arg)
            i = j
        elseif c == '^'
            arg, j = _parse_script_arg(s, _next(s, i))
            _bind_script!(nodes, '^', arg)
            i = j
        elseif c == '-'
            push!(nodes, TSAtom("−", false, false))
            i = _next(s, i)
        else
            push!(nodes, _atom_from_char(c, in_math=true))
            i = _next(s, i)
        end
    end

    return nodes, i
end

function _parse_text_chunk(s::AbstractString)
    nodes = TypesetNode[]
    i = firstindex(s)
    while i <= lastindex(s)
        push!(nodes, TSAtom(string(s[i]), false, false))
        i = _next(s, i)
    end
    return nodes
end

function parse_typeset(s::AbstractString)
    nodes = TypesetNode[]
    isempty(s) && return nodes

    i = firstindex(s)
    start = i

    while i <= lastindex(s)
        if s[i] == '$'
            if start < i
                append!(nodes, _parse_text_chunk(s[start:prevind(s, i)]))
            end

            j = _next(s, i)
            math_nodes, k = _parse_math_seq(s, j; stopchars=Set(['$']))
            if k > lastindex(s) || s[k] != '$'
                append!(nodes, _parse_text_chunk(s[i:lastindex(s)]))
                return nodes
            end

            append!(nodes, math_nodes)
            i = _next(s, k)
            start = i
            continue
        end
        i = _next(s, i)
    end

    if start <= lastindex(s)
        append!(nodes, _parse_text_chunk(s[start:lastindex(s)]))
    end

    return nodes
end

function _measure_text(cc::CairoContext, font::String, text::String, fontsize::Float64, italic::Bool, bold::Bool=false)
    if isempty(text)
        return 0.0, 0.0, 0.0, 0.0
    end

    Cairo.save(cc)
    weight = bold ? Cairo.FONT_WEIGHT_BOLD : Cairo.FONT_WEIGHT_NORMAL
    select_font_face(cc, font, italic ? Cairo.FONT_SLANT_ITALIC : Cairo.FONT_SLANT_NORMAL, weight)
    set_font_size(cc, fontsize)

    te = text_extents(cc, text)
    Cairo.restore(cc)

    width = max(te[5], te[3])
    # Use glyph bearings instead of fixed font fractions so vertical spacing is responsive.
    ascent = max(-te[2], 0.55 * fontsize)
    descent = max(te[4] + te[2], 0.12 * fontsize)
    height = max(te[4], ascent + descent)
    return width, height, ascent, descent
end

function _atom_ink_right(cc::CairoContext, font::String, atom::TSAtom, fontsize::Float64)
    isempty(atom.text) && return 0.0

    Cairo.save(cc)
    weight = atom.bold ? Cairo.FONT_WEIGHT_BOLD : Cairo.FONT_WEIGHT_NORMAL
    select_font_face(cc, font, atom.italic ? Cairo.FONT_SLANT_ITALIC : Cairo.FONT_SLANT_NORMAL, weight)
    set_font_size(cc, fontsize)
    te = text_extents(cc, atom.text)
    Cairo.restore(cc)

    # Right edge of painted glyphs relative to origin.
    return max(te[1] + te[3], 0.0)
end

function _node_width(cc::CairoContext, node::TypesetNode, font::String, fontsize::Float64)
    if node isa TSAtom
        atom = node::TSAtom
        w, _, _, _ = _measure_text(cc, font, atom.text, fontsize, atom.italic, atom.bold)
        return w
    elseif node isa TSGroup
        items = (node::TSGroup).items
        width = 0.0
        for (i, item) in enumerate(items)
            prev_node = i > 1 ? items[i-1] : nothing
            next_node = i < length(items) ? items[i+1] : nothing
            pre, post = _operator_spacing(item, prev_node, next_node, fontsize)
            width += pre + _node_width(cc, item, font, fontsize) + post
        end
        return width
    elseif node isa TSParenGroup
        group = node::TSParenGroup
        body = TSGroup(group.items)
        body_w = _node_width(cc, body, font, fontsize)
        body_min, body_max = _node_bounds(cc, body, font, fontsize)
        body_h = max(body_max - body_min, fontsize)
        paren_xscale = _paren_xscale_for_height(body_h, fontsize)
        paren_size = max(fontsize, 0.9 * body_h)
        left_w_raw, _, _, _ = _measure_text(cc, font, "(", paren_size, false)
        right_w_raw, _, _, _ = _measure_text(cc, font, ")", paren_size, false)
        left_w = left_w_raw * paren_xscale
        right_w = right_w_raw * paren_xscale
        outer_pad = _paren_outer_pad_factor * fontsize
        return left_w + body_w + right_w + 2 * outer_pad
    elseif node isa TSScripts
        s = node::TSScripts
        base_w = _node_width(cc, s.base, font, fontsize)
        trailing_pad = s.base isa TSParenGroup ? _paren_outer_pad_factor * fontsize : 0.0
        base_anchor = if s.base isa TSAtom
            _atom_ink_right(cc, font, s.base::TSAtom, fontsize)
        else
            base_w - trailing_pad
        end
        child_size = 0.7 * fontsize
        sub_w = s.sub === nothing ? 0.0 : _node_width(cc, s.sub, font, child_size)
        sup_w = s.sup === nothing ? 0.0 : _node_width(cc, s.sup, font, child_size)
        return max(base_w, base_anchor + 0.07 * fontsize + max(sub_w, sup_w))
    elseif node isa TSSqrt
        sqrt_node = node::TSSqrt
        body_w = _node_width(cc, sqrt_node.body, font, fontsize)
        body_min, body_max = _node_bounds(cc, sqrt_node.body, font, fontsize)
        body_h = max(body_max - body_min, fontsize)
        root_size = max(fontsize, 0.9 * body_h)
        root_w, _, _, _ = _measure_text(cc, font, "√", root_size, false)
        return root_w + 0.08 * fontsize + body_w
    else
        f = node::TSFraction
        child_size = 0.85 * fontsize
        num_w = _node_width(cc, f.num, font, child_size)
        den_w = _node_width(cc, f.den, font, child_size)
        return max(num_w, den_w) + 0.24 * fontsize
    end
end

function _node_bounds(cc::CairoContext, node::TypesetNode, font::String, fontsize::Float64)
    glyphs = GlyphElem[]
    rules = RuleElem[]
    min_y = Ref(Inf)
    max_y = Ref(-Inf)
    _layout_node!(cc, node, font, glyphs, rules, 0.0, 0.0, fontsize, min_y, max_y)
    if !isfinite(min_y[]) || !isfinite(max_y[])
        return 0.0, 0.0
    end
    return min_y[], max_y[]
end

@inline _is_operator_atom(node::TypesetNode) = node isa TSAtom && ((node::TSAtom).text in _typeset_binary_ops)
@inline _is_opening_delim_atom(node::TypesetNode) = node isa TSAtom && ((node::TSAtom).text in ("(", "[", "{"))

function _operator_spacing(node::TypesetNode, prev_node::Union{Nothing,TypesetNode}, next_node::Union{Nothing,TypesetNode}, fontsize::Float64)
    _is_operator_atom(node) || return 0.0, 0.0
    op = (node::TSAtom).text

    unary_context = prev_node === nothing || _is_operator_atom(prev_node) || _is_opening_delim_atom(prev_node)

    # Unary signs at start of math/group (or after operator/opening delimiter): no extra spacing.
    if unary_context && (op == "+" || op == "−")
        return 0.0, 0.0
    end

    # Only space operators that sit between two operands.
    if prev_node === nothing || next_node === nothing || _is_operator_atom(prev_node) || _is_operator_atom(next_node)
        return 0.0, 0.0
    end

    s = 0.16 * fontsize
    return s, s
end

function _layout_node!(cc::CairoContext, node::TypesetNode, font::String, glyphs::Vector{GlyphElem}, rules::Vector{RuleElem},
    x::Float64, y::Float64, fontsize::Float64, min_y::Base.RefValue{Float64}, max_y::Base.RefValue{Float64})

    if node isa TSAtom
        atom = node::TSAtom
        w, h, ascent, descent = _measure_text(cc, font, atom.text, fontsize, atom.italic, atom.bold)
        push!(glyphs, GlyphElem(atom.text, x, y, w, h, fontsize, atom.italic, atom.bold, 1.0))

        min_y[] = min(min_y[], y - ascent)
        max_y[] = max(max_y[], y + descent)
        return w
    end

    if node isa TSGroup
        group = node::TSGroup
        return _layout_nodes!(cc, group.items, font, glyphs, rules, x, y, fontsize, min_y, max_y)
    end

    if node isa TSParenGroup
        group = node::TSParenGroup
        body = TSGroup(group.items)
        body_w = _node_width(cc, body, font, fontsize)
        body_min, body_max = _node_bounds(cc, body, font, fontsize)
        body_h = max(body_max - body_min, fontsize)
        body_center = 0.5 * (body_min + body_max)
        paren_xscale = _paren_xscale_for_height(body_h, fontsize)

        paren_size = max(fontsize, 0.9 * body_h)
        l_w_raw, l_h, l_a, l_d = _measure_text(cc, font, "(", paren_size, false)
        r_w_raw, r_h, r_a, r_d = _measure_text(cc, font, ")", paren_size, false)
        l_w = l_w_raw * paren_xscale
        r_w = r_w_raw * paren_xscale
        outer_pad = _paren_outer_pad_factor * fontsize
        x0 = x + outer_pad

        l_base = body_center - 0.5 * (l_d - l_a)
        r_base = body_center - 0.5 * (r_d - r_a)

        push!(glyphs, GlyphElem("(", x0, l_base, l_w, l_h, paren_size, false, false, paren_xscale))
        min_y[] = min(min_y[], l_base - l_a)
        max_y[] = max(max_y[], l_base + l_d)

        _layout_nodes!(cc, group.items, font, glyphs, rules, x0 + l_w, y, fontsize, min_y, max_y)

        right_x = x0 + l_w + body_w
        push!(glyphs, GlyphElem(")", right_x, r_base, r_w, r_h, paren_size, false, false, paren_xscale))
        min_y[] = min(min_y[], r_base - r_a)
        max_y[] = max(max_y[], r_base + r_d)

        return l_w + body_w + r_w + 2 * outer_pad
    end

    if node isa TSScripts
        script = node::TSScripts
        base_w = _layout_node!(cc, script.base, font, glyphs, rules, x, y, fontsize, min_y, max_y)
        base_min, base_max = _node_bounds(cc, script.base, font, fontsize)
        base_ascent = max(-base_min, 0.0)
        base_descent = max(base_max, 0.0)

        # If base is a parenthesized group, avoid extra trailing outer pad
        # before scripts to keep ')' tightly attached to sub/sup indices.
        trailing_pad = script.base isa TSParenGroup ? _paren_outer_pad_factor * fontsize : 0.0
        base_anchor = if script.base isa TSAtom
            _atom_ink_right(cc, font, script.base::TSAtom, fontsize)
        else
            base_w - trailing_pad
        end
        script_x = x + base_anchor + 0.07 * fontsize
        sub_w = 0.0
        sup_w = 0.0
        child_size = 0.6 * fontsize

        sup_extra = max(0.0, base_ascent - 0.6 * fontsize)
        sub_extra = max(0.0, base_descent - 0.25 * fontsize)
        group_boost = script.base isa TSParenGroup ? 0.14 * fontsize : script.base isa TSFraction ? 0.08 * fontsize : 0.0
        sup_y = y - (0.42 * fontsize + 0.75 * sup_extra + group_boost)
        sub_y = y + (0.23 * fontsize + 0.25 * sub_extra)

        if script.sub !== nothing
            sub_w = _layout_node!(cc, script.sub, font, glyphs, rules, script_x, sub_y, child_size, min_y, max_y)
        end

        if script.sup !== nothing
            sup_w = _layout_node!(cc, script.sup, font, glyphs, rules, script_x, sup_y, child_size, min_y, max_y)
        end

        return max(base_w, (script_x - x) + max(sub_w, sup_w))
    end

    if node isa TSSqrt
        sqrt_node = node::TSSqrt
        body_w = _node_width(cc, sqrt_node.body, font, fontsize)
        body_min, body_max = _node_bounds(cc, sqrt_node.body, font, fontsize)
        body_h = max(body_max - body_min, fontsize)
        body_center = 0.5 * (body_min + body_max)
        body_ascent = max(-body_min, 0.0)

        root_size = max(fontsize, 0.9 * body_h)
        root_w, root_h, root_a, root_d = _measure_text(cc, font, "√", root_size, false)
        root_base = body_center - 0.5 * (root_d - root_a)
        push!(glyphs, GlyphElem("√", x, root_base, root_w, root_h, root_size, false, false, 1.0))
        min_y[] = min(min_y[], root_base - root_a)
        max_y[] = max(max_y[], root_base + root_d)

        body_x = x + root_w + 0.08 * fontsize
        _layout_node!(cc, sqrt_node.body, font, glyphs, rules, body_x, y, fontsize, min_y, max_y)

        bar_y = y - body_ascent - 0.08 * fontsize
        bar_thickness = max(0.045 * fontsize, 0.2)
        push!(rules, RuleElem(body_x, bar_y, body_w, bar_thickness))
        min_y[] = min(min_y[], bar_y - 0.5 * bar_thickness)

        return root_w + 0.08 * fontsize + body_w
    end

    frac = node::TSFraction
    child_size = 0.85 * fontsize
    num_w = _node_width(cc, frac.num, font, child_size)
    den_w = _node_width(cc, frac.den, font, child_size)
    _, num_max = _node_bounds(cc, frac.num, font, child_size)
    den_min, den_max = _node_bounds(cc, frac.den, font, child_size)
    bar_w = max(num_w, den_w) + 0.24 * fontsize

    num_x = x + 0.5 * (bar_w - num_w)
    den_x = x + 0.5 * (bar_w - den_w)
    # Math axis is usually above baseline; align the fraction bar near '+'/'-' center.
    line_y = y - 0.25 * fontsize
    base_gap = max(0.1 * fontsize, 0.15 * child_size)
    num_gap = base_gap

    # Descender-aware denominator spacing:
    # glyphs with shallow descenders (e.g. "3") need more optical clearance.
    den_descent = max(den_max, 0.0)
    den_ref_descent = max(0.24 * child_size, 1e-9)
    shallow_desc = clamp(1.0 - den_descent / den_ref_descent, 0.0, 1.0)
    den_gap = base_gap + 0.15 * fontsize * shallow_desc

    num_baseline = line_y - num_gap - num_max
    den_baseline = line_y + den_gap - den_min

    _layout_node!(cc, frac.num, font, glyphs, rules, num_x, num_baseline, child_size, min_y, max_y)
    _layout_node!(cc, frac.den, font, glyphs, rules, den_x, den_baseline, child_size, min_y, max_y)

    line_thickness = max(0.05 * fontsize, 0.2)
    push!(rules, RuleElem(x, line_y, bar_w, line_thickness))
    min_y[] = min(min_y[], line_y - 0.5 * line_thickness)
    max_y[] = max(max_y[], line_y + 0.5 * line_thickness)

    return bar_w
end

function _layout_nodes!(cc::CairoContext, nodes::Vector{TypesetNode}, font::String, glyphs::Vector{GlyphElem}, rules::Vector{RuleElem},
    x::Float64, y::Float64, fontsize::Float64, min_y::Base.RefValue{Float64}, max_y::Base.RefValue{Float64})

    current_x = x
    for (i, node) in enumerate(nodes)
        prev_node = i > 1 ? nodes[i-1] : nothing
        next_node = i < length(nodes) ? nodes[i+1] : nothing
        pre, post = _operator_spacing(node, prev_node, next_node, fontsize)
        current_x += pre
        current_x += _layout_node!(cc, node, font, glyphs, rules, current_x, y, fontsize, min_y, max_y)
        current_x += post
    end
    return current_x - x
end

function layout_typeset(cc::CairoContext, text::AbstractString, fontsize::Float64; font::String="NewComputerModern")
    nodes = parse_typeset(text)
    glyphs = GlyphElem[]
    rules = RuleElem[]
    min_y = Ref(Inf)
    max_y = Ref(-Inf)

    width = _layout_nodes!(cc, nodes, font, glyphs, rules, 0.0, 0.0, fontsize, min_y, max_y)
    if !isfinite(min_y[]) || !isfinite(max_y[])
        min_y[] = 0.0
        max_y[] = 0.0
    end
    height = max(max_y[] - min_y[], fontsize)

    return TypesetLayout(glyphs, rules, width, height, min_y[], max_y[])
end

function draw_typeset_layout!(cc::CairoContext, layout::TypesetLayout; font::String="NewComputerModern")
    for glyph in layout.glyphs
        weight = glyph.bold ? Cairo.FONT_WEIGHT_BOLD : Cairo.FONT_WEIGHT_NORMAL
        select_font_face(cc, font, glyph.italic ? Cairo.FONT_SLANT_ITALIC : Cairo.FONT_SLANT_NORMAL, weight)
        set_font_size(cc, glyph.fontsize)
        Cairo.save(cc)
        draw_x = glyph.x
        if glyph.xscale != 1.0
            # Scale around glyph center and compensate translation so visual bounds stay in-place.
            raw_w = glyph.width / glyph.xscale
            draw_x -= 0.5 * raw_w * (1.0 - glyph.xscale)
            Cairo.translate(cc, draw_x, glyph.y)
            Cairo.translate(cc, 0.5 * raw_w, 0.0)
            Cairo.scale(cc, glyph.xscale, 1.0)
            Cairo.translate(cc, -0.5 * raw_w, 0.0)
        else
            Cairo.translate(cc, draw_x, glyph.y)
        end
        move_to(cc, 0.0, 0.0)
        show_text(cc, glyph.text)
        Cairo.restore(cc)
    end
    for rule in layout.rules
        set_line_width(cc, rule.thickness)
        move_to(cc, rule.x, rule.y)
        rel_line_to(cc, rule.width, 0.0)
        stroke(cc)
    end
    return nothing
end
