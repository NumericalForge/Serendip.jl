
const _text_to_unicode = Dict(
    "alpha"    => 'α',
    "beta"     => 'β',
    "gamma"    => 'γ',
    "delta"    => 'δ',
    "epsilon"  => 'ε',
    "varepsilon" => 'ε',
    "zeta"     => 'ζ',
    "eta"      => 'η',
    "theta"    => 'θ',
    "iota"     => 'ι',
    "kappa"    => 'κ',
    "lambda"   => 'λ',
    "mu"       => 'μ',
    "nu"       => 'ν',
    "xi"       => 'ξ',
    "omicron"  => 'ο',
    "pi"       => 'π',
    "rho"      => 'ρ',
    "sigma"    => 'σ',
    "tau"      => 'τ',
    "upsilon"  => 'υ',
    "phi"      => 'φ',
    "chi"      => 'χ',
    "psi"      => 'ψ',
    "omega"    => 'ω',

    "Gamma"    => 'Γ',
    "Delta"    => 'Δ',
    "Theta"    => 'Θ',
    "Lambda"   => 'Λ',
    "Xi"       => 'Ξ',
    "Pi"       => 'Π',
    "Sigma"    => 'Σ',
    "Upsilon"  => 'Υ',
    "Phi"      => 'Φ',
    "Psi"      => 'Ψ',
    "Omega"    => 'Ω',
)



struct Token
    kind::Symbol # :number, :symbol, :text, :literal, :space, :underscore, :caret, :lbracket, :rbracket, :dollar
    text::String # empty for structural tokens

    Token(kind::Symbol, text::String) = new(kind, text)
    Token(kind::Symbol) = new(kind, "")
end


abstract type INode end

struct SeqNode <: INode
    items::Vector{INode}
end

struct BlockNode <: INode
    items::Vector{INode}
    bracket::String
end

struct TextNode <: INode
    s::String      # “σ” or “MPa” or “x”
    italic::Bool   # optional
end

struct MathNode <: INode
    items::Vector{INode}
end

struct OperatorNode <: INode
    s::String
end

struct BracketNode <: INode
    s::String
end

struct SpaceNode <: INode
end

struct ScriptsNode <: INode
    base::INode
    sub::Union{Nothing,INode}
    sup::Union{Nothing,INode}
end

struct GlyphElem
    text::String
    x::Float64
    y::Float64
    width::Float64
    height::Float64
    italic::Bool
end

# Get basic tokens from a string
function tokenize(content::String)
    tokens   = Token[]
    pos      = 1
    next_pos = nextind(content, pos)
    last_pos = lastindex(content)

    while pos <= last_pos
        c = content[pos]
        next_pos = nextind(content, pos, 1) # ~ pos+1
        if next_pos <= last_pos
            next_char = content[next_pos]
        else
            next_char = '\0' # end of string
        end

        if c == '/'
            push!(tokens, string(c))
            pos += nextind(content, pos)
        elseif c == '\"' # string literal
            m = match(r"\"(.*?)\"", content, pos)
            m === nothing && error("Unterminated string literal starting at position $pos")
            text = m.captures[1] # the text inside the quotes
            push!(tokens, Token(:literal, text))
            pos += ncodeunits(text) + 2 # +2 for the quotes
        elseif c in "[{(" # brackets
            push!(tokens, Token(:lbracket, string(c)))
            pos += 1
        elseif c in ")}]" # brackets
            push!(tokens, Token(:rbracket, string(c)))
            pos += 1
        elseif isspace(c)
            m = match(r"\s+", content, pos)
            text = m.match
            push!(tokens, Token(:space, " "))
            pos += ncodeunits(text)
        elseif c == '^'
            push!(tokens, Token(:caret, string(c)))
            pos += 1
        elseif c == '_'
            push!(tokens, Token(:underscore, string(c)))
            pos += 1
        elseif c == '$'
            push!(tokens, Token(:dollar))
            pos += 1
        elseif isletter(c)
            m = match(r"(?:(?![_])\w)+", content, pos) # match word characters but not underscore
            text = m.match
            push!(tokens, Token(:text, text))
            pos += ncodeunits(text)
        elseif isnumeric(c)
            m = match(r"[0-9.]+", content, pos)
            number = m.match
            push!(tokens, Token(:number, number))
            pos += ncodeunits(number)
        else
            # General symbols
            push!(tokens, Token(:symbol, string(c)))
            pos += 1
        end
    end

    tokens
end


function get_token(tokens::Vector{<:Token}, pos::Int)
    if 1 <= pos <= length(tokens) 
        return tokens[pos]
    else
        return Token(:symbol, "\0") # End of input or invalid position
    end
end


# refine tokens according to mode
function refine_tokens(tokens::Vector{<:Token}, pos_ini::Int, pos_end::Int; mode=:TEXT)
    # Parse the tokens into a structured format
    refined = Token[]
    pos = pos_ini
    while pos <= pos_end
        token      = tokens[pos]
        next_token = get_token(tokens, pos+1)

        if mode == :TEXT
            if token.kind == :dollar # Math mode
                cpos = findnext(==(Token(:dollar)), tokens, pos+1)
                cpos === nothing && error("Unmatched dollar sign in input")
                eq_tokens = refine_tokens(tokens, pos+1, cpos-1, mode=:MATH)
                push!(refined, Token(:dollar))
                append!(refined, eq_tokens)
                push!(refined, Token(:dollar))
                pos = cpos + 1
            
            elseif token.kind == :literal
                push!(refined, Token(:text, "\""*token.text*"\""))
                pos += 1
            elseif token.kind == :space
                push!(refined, token)
                pos += 1
            else
                push!(refined, Token(:text, token.text))
                pos += 1
            end
        elseif mode == :MATH
            # refine operators and brackets
            if token.kind == :symbol
                if next_token.text == "="
                    if token.text == ">"
                        token = Token(:operator, ">=")
                    elseif token.text == "<"
                        token = Token(:operator, "<=")
                    elseif token.text == "!"
                        token = Token(:operator, "!=")
                    end
                    pos += 1
                else
                    token = Token(:operator, token.text)
                end
                push!(refined, token)
                pos += 1
            elseif token.kind in (:underscore, :caret)
                push!(refined, Token(:operator, token.text))
                pos += 1
            elseif token.kind == :text
                if length(token.text)==1
                    token = Token(:symbol, token.text)
                elseif token.text in keys(_text_to_unicode)
                    token = Token(:symbol, string(_text_to_unicode[token.text]))
                else
                    token = Token(:literal, token.text)
                end
                push!(refined, token)
                pos += 1
            elseif token.kind == :space # skip spaces
                pos += 1
            else
                push!(refined, token)
                pos += 1
            end
        end

    end

    return refined
end


# Parsing ====================

function get_closure(tokens::Vector{Token}, start::Int)
    opening = get_token(tokens, start)
    left    = opening.text
    right   = left == "(" ? ")" : left == "[" ? "]" : "}"
    
    pos = start
    depth = 0
    while true
        pos > length(tokens) && return nothing

        token = tokens[pos]

        # @show pos, token

        if token.kind == :lbracket && token.text == left
            depth += 1
        elseif token.kind == :rbracket && token.text == right
            depth -= 1
            if depth == 0
                return pos - 1
            end
        end

        pos += 1
    end
end


function bind_scripts(nodes::Vector{INode})
    binded = INode[]
    pos = 1
    while pos <= length(nodes)
        node = nodes[pos]
        if node isa OperatorNode && node.s in ("_", "^")

            # missing base or script
            if isempty(binded) || pos==length(nodes)
                push!(binded, node)
                pos += 1
                continue
            end

            base = pop!(binded)
            arg  = nodes[pos+1]

            if node.s == "_"
                if base isa ScriptsNode && base.sub===nothing
                    node = ScriptsNode(base.base, arg, base.sup)
                else
                    node = ScriptsNode(base, arg, nothing)
                end
            else
                if base isa ScriptsNode && base.sup===nothing
                    node = ScriptsNode(base.base, base.sub, arg)
                else
                    node = ScriptsNode(base, nothing, arg)
                end
            end

            push!(binded, node)
            pos += 2
        else
            push!(binded, node)
            pos += 1
        end
        
    end
    return binded
end


function Base.parse(tokens::Vector{<:Token}, pos_ini::Int, pos_end::Int; mode=:TEXT)
    # Parse the tokens into a structured format
    nodes = INode[]
    pos = pos_ini
    while pos <= pos_end
        token = tokens[pos]

        if mode == :TEXT
            if token==Token(:dollar) # Math mode
                cpos = findnext(==(Token(:dollar)), tokens, pos+1)
                cpos === nothing && error("Unmatched dollar sign")
                seq = parse(tokens, pos+1, cpos-1, mode=:MATH)
                node = MathNode(seq)
                push!(nodes, node)
                pos = cpos + 1 # Skip the dollar signs
            else
                node = TextNode(token.text, false)
                push!(nodes, node)
                pos += 1
            end
        elseif mode == :MATH
            if token.kind == :lbracket
                # find the matching right bracket
                pos_right = get_closure(tokens, pos)
                pos_right === nothing && error("Unmatched left bracket '$(token.text)' while parsing")
                seq = parse(tokens, pos+1, pos_right, mode=:MATH)

                prev_token = get_token(tokens, pos-1)
                bracket = token.text
                if token.text == "("
                    @show prev_token
                    if prev_token==Token(:operator, "_") || prev_token.text==Token(:operator, "^")
                        bracket = ""
                    end
                end

                
                push!(nodes, BlockNode(seq, bracket))
                pos = pos_right + 2 # Skip brackets
            # elseif token.kind in (:underscore, :caret)
            #     prev_token = get_token(tokens, pos-1)

            #     push!(nodes, OperatorNode(token.text))
            #     pos += 1
            elseif token.kind == :operator
                push!(nodes, OperatorNode(token.text))
                pos += 1
            elseif token.kind in (:symbol,)
                push!(nodes, TextNode(token.text, true))
                pos += 1
            else # :number, :literal
                push!(nodes, TextNode(token.text, false))
                pos += 1
            end

        end

    end

    if mode==:MATH
        nodes = bind_scripts(nodes)
    end

    return nodes
end


# TODO: check this function
function layout_nodes(nodes::Vector{INode}, fontsize::Float64, x=0.0, y=0.0)
    glyphs = GlyphElem[]
    current_x = x
    
    for node in nodes
        if node isa TextNode
            # Calculate width based on font metrics (using FreeType or Cairo)
            w, h = get_text_extents(node.s, fontsize, node.italic)
            push!(glyphs, GlyphElem(node.s, current_x, y, w, h, node.italic))
            current_x += w
            
        elseif node isa MathNode
            # Math nodes usually imply italicizing symbols
            math_glyphs = layout_nodes(node.items, fontsize, current_x, y)
            append!(glyphs, math_glyphs)
            current_x = isempty(math_glyphs) ? current_x : math_glyphs[end].x + math_glyphs[end].width
            
        elseif node isa ScriptsNode
            # 1. Layout the base
            base_glyphs = layout_nodes([node.base], fontsize, current_x, y)
            append!(glyphs, base_glyphs)
            base_w = isempty(base_glyphs) ? 0.0 : base_glyphs[end].width
            
            # 2. Layout Subscript (smaller font, shifted down)
            if node.sub !== nothing
                sub_fs = fontsize * 0.7
                sub_y = y + (fontsize * 0.3) 
                sub_glyphs = layout_nodes([node.sub], sub_fs, current_x + base_w, sub_y)
                append!(glyphs, sub_glyphs)
            end
            
            # 3. Layout Superscript (smaller font, shifted up)
            if node.sup !== nothing
                sup_fs = fontsize * 0.7
                sup_y = y - (fontsize * 0.4)
                sup_glyphs = layout_nodes([node.sup], sup_fs, current_x + base_w, sup_y)
                append!(glyphs, sup_glyphs)
            end
            
            # Update X based on the widest script
            # (Simplified: just using base width here for brevity)
            current_x += base_w + (fontsize * 0.5) 
        end
    end
    return glyphs
end


# input = raw"$σ_(x x)^M$ [MPa]"
input = raw"$sigma_y^M$ [MPa]"
tokens = tokenize(input)
display(tokens)
tokens = refine_tokens(tokens, 1, length(tokens))


nodes = parse(tokens, 1, length(tokens))