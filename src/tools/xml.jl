# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export XmlDocument, XmlElement

abstract type XmlNode end


function _xml_attributes(attributes)
    result = OrderedDict{String,String}()
    for (key, value) in attributes
        result[string(key)] = string(value)
    end
    return result
end


"""
    XmlElement(name; attributes=Dict(), children=XmlNode[], content="")

Create an element for Serendip's lightweight XML representation. An element
may contain either child nodes or text. Binary `Vector{UInt8}` content is
supported when the element has `encoding="raw"`.
"""
mutable struct XmlElement <: XmlNode
    name::String
    attributes::OrderedDict{String,String}
    children::Vector{XmlNode}
    content::Union{AbstractString,Vector{UInt8}}

    function XmlElement(
        name::AbstractString;
        attributes::Union{AbstractDict,Tuple}=Dict(),
        children::AbstractVector{<:XmlNode}=XmlNode[],
        content::Union{AbstractString,Vector{UInt8}}="",
    )
        return new(string(name), _xml_attributes(attributes), XmlNode[children...], content)
    end
end


addchild!(node::XmlElement, child::XmlNode) = push!(node.children, child)


mutable struct XmlComment <: XmlNode
    content::String
    XmlComment(content::AbstractString) = new(string(content))
end


haschildren(node::XmlElement) = !isempty(node.children)


mutable struct XmlDocument
    attributes::OrderedDict{String,String}
    children::Vector{XmlNode}
    root::Union{XmlElement,Nothing}

    function XmlDocument(attributes::Union{AbstractDict,Tuple})
        return new(_xml_attributes(attributes), XmlNode[], nothing)
    end

    function XmlDocument(attributes::Union{AbstractDict,Tuple}, root::XmlElement)
        return new(_xml_attributes(attributes), XmlNode[root], root)
    end
end


function addchild!(doc::XmlDocument, child::XmlNode)
    if child isa XmlElement
        doc.root === nothing || error("XmlDocument: a document can only have one root element")
        doc.root = child
    end
    push!(doc.children, child)
    return child
end


# Get a list of all nodes with a given attribute.
function Base.getindex(doc::XmlDocument, p::Pair{String,String})
    doc.root === nothing && return XmlElement[]
    return getindex(doc.root, p)
end


function getallchildren(node::XmlElement)
    collected = XmlElement[]
    function _collect!(parent::XmlElement)
        for child in parent.children
            if child isa XmlElement
                push!(collected, child)
                _collect!(child)
            end
        end
    end
    _collect!(node)
    return collected
end


# Get a node from a nested sequence of names. If siblings share a name, use the last one.
function (node::XmlElement)(names::String...)
    current = node
    for name in names
        found = nothing
        for child in Iterators.reverse(current.children)
            if child isa XmlElement && child.name == name
                found = child
                break
            end
        end
        found === nothing && return nothing
        current = found
    end
    return current
end


# Get all immediate child elements with a given name.
function Base.getindex(node::XmlElement, name::String)
    return XmlElement[
        child for child in node.children if child isa XmlElement && child.name == name
    ]
end


# Get a child according to index.
function Base.getindex(node::XmlElement, index::Int)
    checkbounds(Bool, node.children, index) || return nothing
    return node.children[index]
end


function getchild(node::XmlElement, name::String)
    for child in node.children
        child isa XmlElement && child.name == name && return child
    end
    return nothing
end


# Get all descendant elements with a given attribute and value.
function Base.getindex(node::XmlElement, p::Pair{String,String})
    nodes = XmlElement[]
    attribute, value = p

    get(node.attributes, attribute, nothing) == value && push!(nodes, node)
    for child in node.children
        child isa XmlElement || continue
        append!(nodes, getindex(child, p))
    end
    return nodes
end


mutable struct _XmlParser
    data::Vector{UInt8}
    pos::Int
end


_is_xml_whitespace(byte::UInt8) = byte in (0x09, 0x0a, 0x0d, 0x20)


function _is_name_start(char::Char)
    value = Int(char)
    return char in (':', '_') ||
           'A' <= char <= 'Z' ||
           'a' <= char <= 'z' ||
           0x00c0 <= value <= 0x00d6 ||
           0x00d8 <= value <= 0x00f6 ||
           0x00f8 <= value <= 0x02ff ||
           0x0370 <= value <= 0x037d ||
           0x037f <= value <= 0x1fff ||
           0x200c <= value <= 0x200d ||
           0x2070 <= value <= 0x218f ||
           0x2c00 <= value <= 0x2fef ||
           0x3001 <= value <= 0xd7ff ||
           0xf900 <= value <= 0xfdcf ||
           0xfdf0 <= value <= 0xfffd ||
           0x10000 <= value <= 0xeffff
end


function _is_name_char(char::Char)
    value = Int(char)
    return _is_name_start(char) ||
           char in ('-', '.') ||
           '0' <= char <= '9' ||
           value == 0x00b7 ||
           0x0300 <= value <= 0x036f ||
           0x203f <= value <= 0x2040
end


_is_name_delimiter(byte::UInt8) = _is_xml_whitespace(byte) ||
                                  byte in (UInt8('='), UInt8('/'), UInt8('>'), UInt8('?'), UInt8('<'))


function _xml_error(parser::_XmlParser, message::AbstractString)
    line = 1
    column = 1
    for index in 1:min(parser.pos - 1, length(parser.data))
        if parser.data[index] == UInt8('\n')
            line += 1
            column = 1
        else
            column += 1
        end
    end
    error("XmlDocument: $message at line $line, column $column")
end


function _starts_with(parser::_XmlParser, token::AbstractString)
    bytes = codeunits(token)
    parser.pos + length(bytes) - 1 <= length(parser.data) || return false
    for (offset, byte) in enumerate(bytes)
        parser.data[parser.pos + offset - 1] == byte || return false
    end
    return true
end


function _consume!(parser::_XmlParser, token::AbstractString)
    _starts_with(parser, token) || _xml_error(parser, "expected $(repr(token))")
    parser.pos += ncodeunits(token)
    return nothing
end


function _skip_whitespace!(parser::_XmlParser)
    while parser.pos <= length(parser.data) && _is_xml_whitespace(parser.data[parser.pos])
        parser.pos += 1
    end
    return nothing
end


function _string(parser::_XmlParser, first::Int, last::Int)
    first > last && return ""
    value = String(copy(@view parser.data[first:last]))
    isvalid(value) || _xml_error(parser, "invalid UTF-8 text")
    return value
end


function _read_name!(parser::_XmlParser)
    parser.pos <= length(parser.data) || _xml_error(parser, "expected an XML name")

    first = parser.pos
    while parser.pos <= length(parser.data) && !_is_name_delimiter(parser.data[parser.pos])
        parser.pos += 1
    end
    name = _string(parser, first, parser.pos - 1)
    _is_valid_xml_name(name) || _xml_error(parser, "invalid XML name $(repr(name))")
    return name
end


function _read_attribute!(parser::_XmlParser, attributes::OrderedDict{String,String})
    name = _read_name!(parser)
    haskey(attributes, name) && _xml_error(parser, "duplicate attribute $(repr(name))")

    _skip_whitespace!(parser)
    parser.pos <= length(parser.data) && parser.data[parser.pos] == UInt8('=') ||
        _xml_error(parser, "expected '=' after attribute $(repr(name))")
    parser.pos += 1
    _skip_whitespace!(parser)

    parser.pos <= length(parser.data) || _xml_error(parser, "missing value for attribute $(repr(name))")
    quote_byte = parser.data[parser.pos]
    quote_byte in (UInt8('\''), UInt8('"')) ||
        _xml_error(parser, "attribute $(repr(name)) must use quotes")
    parser.pos += 1

    first = parser.pos
    while parser.pos <= length(parser.data) && parser.data[parser.pos] != quote_byte
        parser.data[parser.pos] == UInt8('<') &&
            _xml_error(parser, "'<' is not allowed in an attribute value")
        parser.pos += 1
    end
    parser.pos <= length(parser.data) || _xml_error(parser, "unclosed attribute $(repr(name))")

    attributes[name] = _xml_unescape(_string(parser, first, parser.pos - 1))
    parser.pos += 1
    return nothing
end


function _parse_declaration!(parser::_XmlParser)
    attributes = OrderedDict{String,String}()
    _consume!(parser, "<?xml")
    parser.pos <= length(parser.data) && _is_xml_whitespace(parser.data[parser.pos]) ||
        _xml_error(parser, "expected whitespace after '<?xml'")

    while true
        _skip_whitespace!(parser)
        if _starts_with(parser, "?>")
            _consume!(parser, "?>")
            return attributes
        end
        _read_attribute!(parser, attributes)
    end
end


function _parse_comment!(parser::_XmlParser)
    _consume!(parser, "<!--")
    closing = collect(codeunits("-->"))
    range = findnext(closing, parser.data, parser.pos)
    range === nothing && _xml_error(parser, "unclosed comment")

    content = _string(parser, parser.pos, first(range) - 1)
    _validate_xml_characters(content, "XML comment")
    (occursin("--", content) || endswith(content, "-")) &&
        _xml_error(parser, "'--' is not allowed inside an XML comment")
    parser.pos = last(range) + 1
    return XmlComment(content)
end


function _has_non_whitespace(bytes::Vector{UInt8})
    return any(byte -> !_is_xml_whitespace(byte), bytes)
end


function _parse_raw_content!(parser::_XmlParser, name::String)
    closing = collect(codeunits("</$name>"))
    range = findnext(closing, parser.data, parser.pos)
    range === nothing && _xml_error(parser, "missing closing tag </$name>")

    content_start = parser.pos
    while content_start < first(range) && _is_xml_whitespace(parser.data[content_start])
        content_start += 1
    end
    content_start < first(range) && parser.data[content_start] == UInt8('_') ||
        _xml_error(parser, "raw XML content must start with '_'")

    content = copy(parser.data[content_start + 1:first(range) - 1])
    parser.pos = last(range) + 1
    return content
end


function _parse_element!(parser::_XmlParser)
    _consume!(parser, "<")
    name = _read_name!(parser)
    attributes = OrderedDict{String,String}()

    self_closing = false
    while true
        _skip_whitespace!(parser)
        if _starts_with(parser, "/>")
            _consume!(parser, "/>")
            self_closing = true
            break
        elseif _starts_with(parser, ">")
            _consume!(parser, ">")
            break
        end
        _read_attribute!(parser, attributes)
    end

    self_closing && return XmlElement(name, attributes=attributes)

    if get(attributes, "encoding", "") == "raw"
        content = _parse_raw_content!(parser, name)
        return XmlElement(name, attributes=attributes, content=content)
    end

    children = XmlNode[]
    text = UInt8[]
    while true
        parser.pos <= length(parser.data) || _xml_error(parser, "missing closing tag </$name>")

        if _starts_with(parser, "</")
            _consume!(parser, "</")
            closing_name = _read_name!(parser)
            _skip_whitespace!(parser)
            _consume!(parser, ">")
            closing_name == name ||
                _xml_error(parser, "expected </$name>, got </$closing_name>")
            break
        elseif _starts_with(parser, "<!--")
            _has_non_whitespace(text) &&
                _xml_error(parser, "mixed text and child nodes are not supported")
            empty!(text)
            push!(children, _parse_comment!(parser))
        elseif _starts_with(parser, "<![CDATA[")
            _xml_error(parser, "CDATA sections are not supported")
        elseif _starts_with(parser, "<!DOCTYPE")
            _xml_error(parser, "DOCTYPE declarations are not supported")
        elseif _starts_with(parser, "<?")
            _xml_error(parser, "processing instructions are not supported")
        elseif parser.data[parser.pos] == UInt8('<')
            _has_non_whitespace(text) &&
                _xml_error(parser, "mixed text and child elements are not supported")
            empty!(text)
            push!(children, _parse_element!(parser))
        else
            next_tag = findnext(==(UInt8('<')), parser.data, parser.pos)
            last = next_tag === nothing ? length(parser.data) : next_tag - 1
            append!(text, @view parser.data[parser.pos:last])
            parser.pos = last + 1
        end
    end

    if isempty(children)
        value = String(copy(text))
        isvalid(value) || _xml_error(parser, "invalid UTF-8 text")
        content = _xml_unescape(value)
        return XmlElement(name, attributes=attributes, content=content)
    end

    _has_non_whitespace(text) &&
        _xml_error(parser, "mixed child nodes and text are not supported")
    return XmlElement(name, attributes=attributes, children=children)
end


function _parse_document(data::Vector{UInt8})
    parser = _XmlParser(data, 1)
    if length(data) >= 3 && data[1:3] == UInt8[0xef, 0xbb, 0xbf]
        parser.pos = 4
    end

    _skip_whitespace!(parser)
    attributes = _starts_with(parser, "<?xml") ? _parse_declaration!(parser) : OrderedDict{String,String}()
    _validate_xml_declaration(attributes)

    children = XmlNode[]
    root = nothing
    while true
        _skip_whitespace!(parser)
        parser.pos > length(parser.data) && break

        if _starts_with(parser, "<!--")
            push!(children, _parse_comment!(parser))
        elseif _starts_with(parser, "<!DOCTYPE")
            _xml_error(parser, "DOCTYPE declarations are not supported")
        elseif _starts_with(parser, "<?")
            _xml_error(parser, "processing instructions are not supported")
        elseif parser.data[parser.pos] == UInt8('<')
            root === nothing || _xml_error(parser, "a document can only have one root element")
            root = _parse_element!(parser)
            push!(children, root)
        else
            _xml_error(parser, "text is not allowed outside the root element")
        end
    end

    root === nothing && _xml_error(parser, "no root element found")
    document = XmlDocument(attributes)
    document.children = children
    document.root = root
    return document
end


function _validate_xml_declaration(attributes::AbstractDict)
    isempty(attributes) && return nothing
    haskey(attributes, "version") || error("XmlDocument: an XML declaration requires a version")
    attributes["version"] == "1.0" ||
        error("XmlDocument: only XML version 1.0 is supported")

    encoding = lowercase(get(attributes, "encoding", "utf-8"))
    encoding in ("utf-8", "utf8") ||
        error("XmlDocument: only UTF-8 encoding is supported")
    return nothing
end


"""
    XmlDocument(input::String)

Read a UTF-8 XML string or file using Serendip's lightweight XML subset. The
parser supports declarations, comments, XML 1.0 element and attribute names,
text-only elements, nested elements, self-closing tags, and VTK raw appended
data. Mixed content, CDATA, DOCTYPE declarations, and additional processing
instructions are not supported.
"""
function XmlDocument(input::String)
    if isfile(input)
        data = read(input)
    else
        data = collect(codeunits(input))
        check = _XmlParser(data, 1)
        length(data) >= 3 && data[1:3] == UInt8[0xef, 0xbb, 0xbf] && (check.pos = 4)
        _skip_whitespace!(check)
        if check.pos > length(data) || data[check.pos] != UInt8('<')
            data = read(input)
        end
    end
    return _parse_document(data)
end


function _is_valid_xml_name(name::AbstractString)
    isempty(name) && return false
    indices = eachindex(name)
    first_index = first(indices)
    _is_name_start(name[first_index]) || return false
    return all(index -> _is_name_char(name[index]), Iterators.drop(indices, 1))
end


function writenode(io::IO, node::XmlComment, level::Int)
    _validate_xml_characters(node.content, "XML comment")
    (occursin("--", node.content) || endswith(node.content, "-")) &&
        error("XmlComment: '--' and a trailing '-' are not allowed in comments")
    println(io, "   "^level, "<!--", node.content, "-->")
    return nothing
end


function writenode(io::IO, node::XmlElement, level::Int)
    _is_valid_xml_name(node.name) || error("XmlElement: invalid name $(repr(node.name))")
    !isempty(node.children) && !isempty(node.content) &&
        error("XmlElement: $(repr(node.name)) cannot contain both text and child nodes")

    indent = "   "^level
    print(io, indent, "<", node.name)
    for (attribute, value) in node.attributes
        _is_valid_xml_name(attribute) ||
            error("XmlElement: invalid attribute name $(repr(attribute))")
        print(io, " $attribute=\"$(_xml_escape(value))\"")
    end

    if isempty(node.children) && isempty(node.content)
        println(io, "/>")
        return nothing
    end

    if get(node.attributes, "encoding", "") == "raw"
        isempty(node.children) || error("XmlElement: raw content cannot have child nodes")
        print(io, ">_")
        write(io, node.content)
        println(io, "</", node.name, ">")
        return nothing
    end

    if !isempty(node.children)
        println(io, ">")
        for child in node.children
            writenode(io, child, level + 1)
        end
        println(io, indent, "</", node.name, ">")
        return nothing
    end

    node.content isa AbstractString ||
        error("XmlElement: binary content requires encoding=\"raw\"")
    escaped = _xml_escape_text(node.content)
    if node.name == "DataArray"
        println(io, ">")
        for line in split(escaped, '\n', keepempty=true)
            println(io, "   "^(level + 1), line)
        end
        println(io, indent, "</", node.name, ">")
    else
        print(io, ">", escaped)
        println(io, "</", node.name, ">")
    end
    return nothing
end


function save(doc::XmlDocument, filename::String)
    doc.root === nothing && error("XmlDocument: cannot save a document without a root element")
    _validate_xml_declaration(doc.attributes)
    open(filename, "w") do io
        if !isempty(doc.attributes)
            print(io, "<?xml version=\"$(_xml_escape(doc.attributes["version"]))\"")
            for (attribute, value) in doc.attributes
                attribute == "version" && continue
                _is_valid_xml_name(attribute) ||
                    error("XmlDocument: invalid declaration attribute $(repr(attribute))")
                print(io, " $attribute=\"$(_xml_escape(value))\"")
            end
            println(io, "?>")
        end

        for child in doc.children
            writenode(io, child, 0)
        end
    end
    return nothing
end


function _is_xml_codepoint(value::Integer)
    return value in (0x09, 0x0a, 0x0d) ||
           0x20 <= value <= 0xd7ff ||
           0xe000 <= value <= 0xfffd ||
           0x10000 <= value <= 0x10ffff
end


function _validate_xml_characters(value::AbstractString, context::AbstractString)
    all(char -> _is_xml_codepoint(Int(char)), value) ||
        error("$context contains an invalid XML character")
    return nothing
end


function _xml_escape_text(value::AbstractString)
    _validate_xml_characters(value, "XML text")
    io = IOBuffer()
    for char in value
        if char == '&'
            print(io, "&amp;")
        elseif char == '<'
            print(io, "&lt;")
        elseif char == '>'
            print(io, "&gt;")
        else
            print(io, char)
        end
    end
    return String(take!(io))
end


function _xml_escape(value::AbstractString)
    escaped = _xml_escape_text(value)
    escaped = replace(escaped, "\"" => "&quot;")
    return replace(escaped, "'" => "&apos;")
end


function _xml_entity(entity::AbstractString)
    entity == "amp" && return "&"
    entity == "lt" && return "<"
    entity == "gt" && return ">"
    entity == "quot" && return "\""
    entity == "apos" && return "'"

    value = if startswith(entity, "#x")
        tryparse(Int, entity[3:end], base=16)
    elseif startswith(entity, "#")
        tryparse(Int, entity[2:end])
    else
        nothing
    end
    value === nothing && error("XmlDocument: unknown XML entity &$(entity);")
    _is_xml_codepoint(value) || error("XmlDocument: invalid numeric XML entity &$(entity);")
    return string(Char(value))
end


function _xml_unescape(value::AbstractString)
    isempty(value) && return ""
    io = IOBuffer()
    pos = firstindex(value)
    while pos <= lastindex(value)
        amp = findnext('&', value, pos)
        if amp === nothing
            print(io, SubString(value, pos))
            break
        end
        amp > pos && print(io, SubString(value, pos, prevind(value, amp)))

        entity_start = nextind(value, amp)
        semicolon = findnext(';', value, entity_start)
        semicolon === nothing && error("XmlDocument: unclosed XML entity")
        entity = SubString(value, entity_start, prevind(value, semicolon))
        print(io, _xml_entity(entity))
        pos = nextind(value, semicolon)
    end
    result = String(take!(io))
    _validate_xml_characters(result, "XML text")
    return result
end


export to_xml_node
const Xsingletype = Union{Number,Bool,AbstractString,Char,Symbol}


function to_xml_node(
    value::Xsingletype,
    name::String="Value",
    attributes::AbstractDict=Dict{String,String}(),
)
    node_attributes = _xml_attributes(attributes)
    node_attributes["type"] = string(typeof(value))
    return XmlElement(name, attributes=node_attributes, content=string(value))
end


function to_xml_node(
    array::AbstractArray,
    name::String="Array",
    attributes::AbstractDict=Dict{String,String}(),
)
    node_attributes = _xml_attributes(attributes)
    if eltype(array) <: Xsingletype
        node_attributes["format"] = "ascii"
        node_attributes["type"] = string(eltype(array))
        node_attributes["components"] = string(size(array, 2))
        if isempty(array)
            content = ""
        else
            representation = repr("text/plain", array)
            newline = findfirst('\n', representation)
            content = newline === nothing ? representation : representation[nextind(representation, newline):end]
            content = replace(content, r"^ "m => "")
        end
        return XmlElement(name, attributes=node_attributes, content=content)
    end

    node_attributes["type"] = string(eltype(array))
    children = XmlNode[to_xml_node(item) for item in array]
    return XmlElement(name, attributes=node_attributes, children=children)
end


function to_xml_node(
    dict::AbstractDict,
    name::String="Dict",
    attributes::AbstractDict=Dict{String,String}(),
)
    node_attributes = _xml_attributes(attributes)
    node_attributes["type"] = split(string(typeof(dict)), ".")[end]
    children = XmlNode[
        to_xml_node(collect(keys(dict)), "keys"),
        to_xml_node(collect(values(dict)), "values"),
    ]
    return XmlElement(name, attributes=node_attributes, children=children)
end


function to_xml_node(obj::Any, name::String=""; exclude::Vector{Symbol}=Symbol[])
    isempty(name) && (name = string(nameof(typeof(obj))))
    attributes = OrderedDict{String,String}()
    children = XmlNode[]

    for field in fieldnames(typeof(obj))
        field in exclude && continue
        field_string = string(field)
        startswith(field_string, "_") && continue

        value = getfield(obj, field)
        if value isa Xsingletype
            attributes[field_string] = string(value)
        else
            push!(children, to_xml_node(value, field_string))
        end
    end
    return XmlElement(name, attributes=attributes, children=children)
end


function XmlElement(obj::Any, name::String=""; exclude::Vector{Symbol}=Symbol[])
    return to_xml_node(obj, name, exclude=exclude)
end
