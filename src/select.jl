# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export any_of, none_of

abstract type SelectorCombiner end

struct AnyOf <: SelectorCombiner
    selectors::Tuple
end

struct NoneOf <: SelectorCombiner
    selectors::Tuple
end

any_of(selectors...) = AnyOf(selectors)
none_of(selectors...) = NoneOf(selectors)


_selector_contains_redirect(selector::Symbol) = selector in (:node, :ip)
_selector_contains_redirect(selector::Tuple) = any(_selector_contains_redirect, selector)
_selector_contains_redirect(selector::SelectorCombiner) = any(_selector_contains_redirect, selector.selectors)
_selector_contains_redirect(::Any) = false


function _validate_selector_combiner(selector::SelectorCombiner)
    _selector_contains_redirect(selector) || return nothing
    throw(ArgumentError("select: any_of(...) and none_of(...) cannot contain :node or :ip selectors"))
end


_selector_branch_selectors(branch) = _flatten_selectors((branch,))
_selector_kind(::Type{Node}) = "node"
_selector_kind(::Type{Ip}) = "ip"


function _filter_selected_by_mask(selected::Vector{Int}, mask::AbstractVector{Bool})
    return Int[i for i in selected if mask[i]]
end


function _exclude_selected_by_mask(selected::Vector{Int}, mask::AbstractVector{Bool})
    return Int[i for i in selected if !mask[i]]
end


function _invert_selected(nitems::Int, selected::Vector{Int})
    mask = falses(nitems)
    for idx in selected
        mask[idx] = true
    end
    return _exclude_selected_by_mask(collect(1:nitems), mask)
end


function _tag_selected!(items, selected::Vector{Int}, tag::Union{String,Nothing})
    tag === nothing && return nothing
    for i in selected
        items[i].tag = tag
    end
    return nothing
end


function _select_combiner_indices(eval_branch, nitems::Int, selector::AnyOf, selected::Vector{Int})
    _validate_selector_combiner(selector)
    mask = falses(nitems)
    for branch in selector.selectors
        branch_selected = eval_branch(_selector_branch_selectors(branch), copy(selected))
        for idx in branch_selected
            mask[idx] = true
        end
    end
    return _filter_selected_by_mask(selected, mask)
end


function _select_combiner_indices(eval_branch, nitems::Int, selector::NoneOf, selected::Vector{Int})
    excluded = _select_combiner_indices(eval_branch, nitems, AnyOf(selector.selectors), selected)
    mask = falses(nitems)
    for idx in excluded
        mask[idx] = true
    end
    return _exclude_selected_by_mask(selected, mask)
end


function _select_point_indices_seq(
    items::Vector{T},
    selectors,
    selected::Vector{Int};
    nearest::Bool=false,
    quiet::Bool=false,
    prefix::AbstractString="select",
) where {T<:Union{Node,Ip}}
    for selector in selectors
        selected = _select_point_indices(items, selector, selected; nearest=nearest, quiet=quiet, prefix=prefix)
    end
    return selected
end


function _select_point_indices(
    items::Vector{T},
    selector::SelectorCombiner,
    selected::Vector{Int};
    nearest::Bool=false,
    quiet::Bool=false,
    prefix::AbstractString="select",
) where {T<:Union{Node,Ip}}
    return _select_combiner_indices(length(items), selector, selected) do branch_selectors, branch_selected
        _select_point_indices_seq(items, branch_selectors, branch_selected; nearest=nearest, quiet=quiet, prefix=prefix)
    end
end


function _select_point_indices(
    ::Vector{T},
    selector::Symbol,
    selected::Vector{Int};
    nearest::Bool=false,
    quiet::Bool=false,
    prefix::AbstractString="select",
) where {T<:Union{Node,Ip}}
    if selector == :all
        return selected
    elseif selector == :none
        return Int[]
    else
        error("select: unknown symbol selector $(repr(selector))")
    end
end


function _select_point_indices(
    items::Vector{T},
    selector::Union{Expr,Symbolic},
    selected::Vector{Int};
    nearest::Bool=false,
    quiet::Bool=false,
    prefix::AbstractString="select",
) where {T<:Union{Node,Ip}}
    filtered = items[selected]
    mask = Bool[]
    for item in filtered
        x, y, z = item.coord.x, item.coord.y, item.coord.z
        push!(mask, evaluate(selector, x=x, y=y, z=z))
    end
    return selected[mask]
end


function _select_point_indices(
    items::Vector{T},
    selector::String,
    selected::Vector{Int};
    nearest::Bool=false,
    quiet::Bool=false,
    prefix::AbstractString="select",
) where {T<:Union{Node,Ip}}
    return Int[i for i in selected if items[i].tag == selector]
end


function _select_point_indices(
    items::Vector{T},
    selector::AbstractVector{<:Real},
    selected::Vector{Int};
    nearest::Bool=false,
    quiet::Bool=false,
    prefix::AbstractString="select",
) where {T<:Union{Node,Ip}}
    kind = _selector_kind(T)
    X = Vec3(selector)
    exact = Int[i for i in selected if norm(items[i].coord - X) < 1e-8]

    if !isempty(exact)
        return exact
    elseif nearest && !isempty(selected)
        candidates = items[selected]
        nearest_item = Serendip.nearest(candidates, X)
        i = findfirst(isequal(nearest_item), candidates)
        resolved = Int[selected[i]]
        if !quiet
            selector_str = replace(string(selector), r"(?<!\,)\s+" => "")
            Xn = round.(items[resolved[1]].coord, sigdigits=4)
            notify("$prefix: No $kind found at $selector_str. Using nearest at $Xn")
        end
        return resolved
    else
        if !quiet
            selector_str = replace(string(selector), r"(?<!\,)\s+" => "")
            notify("$prefix: No $kind found at $selector_str")
        end
        return Int[]
    end
end


function _select_point_indices(
    ::Vector{T},
    selector,
    ::Vector{Int};
    nearest::Bool=false,
    quiet::Bool=false,
    prefix::AbstractString="select",
) where {T<:Union{Node,Ip}}
    error("select: unknown selector type $(typeof(selector))")
end


function _select_points(
    items::Vector{T},
    selectors...;
    invert=false,
    nearest=false,
    quiet=false,
    prefix::AbstractString="select",
    tag::Union{String,Nothing}=nothing,
) where {T<:Union{Node,Ip}}
    selectors = _flatten_selectors(selectors)
    selected = _select_point_indices_seq(items, selectors, collect(eachindex(items)); nearest=nearest, quiet=quiet, prefix=prefix)
    invert && (selected = _invert_selected(length(items), selected))
    _tag_selected!(items, selected, tag)
    return items[selected]
end


"""
    select(nodes::Vector{Node}, selectors...; invert=false, nearest=false, quiet=false, prefix="select", tag=nothing)

Select nodes from a collection using one or more filters. Filters are applied
sequentially (logical AND), starting from all nodes. Tuple selectors are
flattened and applied in order.

Supported selectors:
- `:all` keeps the current selection unchanged.
- `:none` clears the selection.
- `String` matches `node.tag`.
- `Expr` or `Symbolic` evaluates a coordinate condition using `x`, `y`, `z`.
- `AbstractVector{<:Real}` selects by point coordinates.
- `any_of(...)` unions multiple node-selector branches.
- `none_of(...)` excludes the union of multiple node-selector branches.
"""
function select(
    nodes::Vector{Node},
    selectors...;
    invert=false,
    nearest=false,
    quiet=false,
    prefix::AbstractString="select",
    tag::Union{String,Nothing}=nothing,
)
    return _select_points(nodes, selectors...; invert=invert, nearest=nearest, quiet=quiet, prefix=prefix, tag=tag)
end


"""
    select(ips::Vector{Ip}, selectors...; invert=false, nearest=false, quiet=false, prefix="select", tag=nothing)

Select integration points from a collection using one or more filters. Filters
are applied sequentially (logical AND), starting from all integration points.
Tuple selectors are flattened and applied in order.

Supported selectors:
- `:all` keeps the current selection unchanged.
- `:none` clears the selection.
- `String` matches `ip.tag`.
- `Expr` or `Symbolic` evaluates a coordinate condition using `x`, `y`, `z`.
- `AbstractVector{<:Real}` selects by point coordinates.
- `any_of(...)` unions multiple IP-selector branches.
- `none_of(...)` excludes the union of multiple IP-selector branches.
"""
function select(
    ips::Vector{Ip},
    selectors...;
    invert=false,
    nearest=false,
    quiet=false,
    prefix::AbstractString="select",
    tag::Union{String,Nothing}=nothing,
)
    return _select_points(ips, selectors...; invert=invert, nearest=nearest, quiet=quiet, prefix=prefix, tag=tag)
end


function _select_elem_indices_seq(
    elems::Vector{<:AbstractCell},
    selectors,
    selected::Vector{Int};
    nearest::Bool=false,
    quiet::Bool=false,
    prefix::AbstractString="select",
)
    for selector in selectors
        selected = _select_elem_indices(elems, selector, selected; nearest=nearest, quiet=quiet, prefix=prefix)
    end
    return selected
end


function _select_elem_indices(
    elems::Vector{<:AbstractCell},
    selector::SelectorCombiner,
    selected::Vector{Int};
    nearest::Bool=false,
    quiet::Bool=false,
    prefix::AbstractString="select",
)
    return _select_combiner_indices(length(elems), selector, selected) do branch_selectors, branch_selected
        _select_elem_indices_seq(elems, branch_selectors, branch_selected; nearest=nearest, quiet=quiet, prefix=prefix)
    end
end


function _select_elem_indices(
    elems::Vector{<:AbstractCell},
    selector::Symbol,
    selected::Vector{Int};
    nearest::Bool=false,
    quiet::Bool=false,
    prefix::AbstractString="select",
)
    if selector == :all
        return selected
    elseif selector == :none
        return Int[]
    elseif selector in (:bulk, :solid, :line, :surface, :cohesive, :contact, :line_interface, :tip)
        role = selector == :bulk ? :solid : selector
        return Int[i for i in selected if elems[i].role == role]
    elseif selector == :active
        return Int[i for i in selected if elems[i].active]
    elseif selector == :embedded
        return Int[i for i in selected if elems[i].role == :line && length(elems[i].couplings) > 0]
    elseif selector in (:node, :ip)
        throw(ArgumentError("select: :node and :ip are only supported as top-level element selectors, not inside any_of(...) or none_of(...)"))
    else
        error("select: cannot filter array of Cell with symbol $(repr(selector))")
    end
end


function _select_elem_indices(
    elems::Vector{<:AbstractCell},
    selector::String,
    selected::Vector{Int};
    nearest::Bool=false,
    quiet::Bool=false,
    prefix::AbstractString="select",
)
    return Int[i for i in selected if elems[i].tag == selector]
end


function _select_elem_indices(
    elems::Vector{<:AbstractCell},
    selector::Union{Expr,Symbolic},
    selected::Vector{Int};
    nearest::Bool=false,
    quiet::Bool=false,
    prefix::AbstractString="select",
)
    isempty(selected) && return Int[]

    nodes = [node for i in selected for node in elems[i].nodes]
    max_id = maximum(node -> node.id, nodes)
    pointmap = zeros(Int, max_id)
    truth = Bool[]

    for (i, node) in enumerate(nodes)
        pointmap[node.id] = i
        x, y, z = node.coord.x, node.coord.y, node.coord.z
        push!(truth, evaluate(selector, x=x, y=y, z=z))
    end

    return Int[i for i in selected if all(truth[pointmap[node.id]] for node in elems[i].nodes)]
end


function _select_elem_indices(
    elems::Vector{<:AbstractCell},
    selector::Vector{Int},
    selected::Vector{Int};
    nearest::Bool=false,
    quiet::Bool=false,
    prefix::AbstractString="select",
)
    return intersect(selected, selector)
end


function _select_elem_indices(
    ::Vector{<:AbstractCell},
    selector,
    ::Vector{Int};
    nearest::Bool=false,
    quiet::Bool=false,
    prefix::AbstractString="select",
)
    error("select: unknown selector type $(typeof(selector))")
end


"""
    select(elems::Vector{<:AbstractCell}, selectors...; invert=false, nearest=false, quiet=false, prefix="select", tag=nothing)

Select entities from an element list using one or more filters. By default this
function filters cells, but `:node` and `:ip` can redirect selection to nodes or
integration points. Filters are applied sequentially (logical AND), starting
from all indices. Tuple selectors are flattened and applied as independent
filters.

Supported selectors:
- `Symbol`:
  - `:all` keeps the current selection unchanged.
  - `:none` clears the selection.
  - `:solid`, `:bulk` (alias of `:solid`), `:line`, `:surface`, `:contact`,
    `:cohesive`, `:line_interface`, `:tip` filter by `cell.role`.
  - `:node` returns `select(get_nodes(elems), ...)` using the remaining selectors.
  - `:ip` returns `select(get_ips(elems), ...)` using the remaining selectors.
  - `:active` keeps only active cells.
  - `:embedded` keeps embedded line cells (`role == :line` with couplings).
- `String`: keep cells whose `tag` matches exactly.
- `Expr` or `Symbolic`: keep cells whose nodes all satisfy the coordinate
  expression (variables `x`, `y`, `z`).
- `Vector{Int}`: intersect the current selection with explicit element indices.
- `any_of(...)`: union multiple element-selector branches.
- `none_of(...)`: exclude the union of multiple element-selector branches.
"""
function select(
    elems::Vector{<:AbstractCell},
    selectors...;
    invert=false,
    nearest=false,
    quiet=false,
    prefix::AbstractString="select",
    tag::Union{String,Nothing}=nothing,
)
    selectors = _flatten_selectors(selectors)
    selected = collect(eachindex(elems))

    for (i, selector) in enumerate(selectors)
        if selector === :ip
            return select(get_ips(elems[selected]), selectors[i + 1:end]...; invert=invert, nearest=nearest, quiet=quiet, prefix=prefix, tag=tag)
        elseif selector === :node
            return select(get_nodes(elems[selected]), selectors[i + 1:end]...; invert=invert, nearest=nearest, quiet=quiet, prefix=prefix, tag=tag)
        else
            selected = _select_elem_indices(elems, selector, selected; nearest=nearest, quiet=quiet, prefix=prefix)
        end
    end

    invert && (selected = _invert_selected(length(elems), selected))
    _tag_selected!(elems, selected, tag)
    return elems[selected]
end


"""
    select(domain::AbstractDomain, kind::Symbol, selectors...; invert=false, nearest=false, quiet=false, prefix="select", tag=nothing)

Filters entities from a finite element domain (`domain`) by type and selection criteria.
"""
function select(
    domain::AbstractDomain,
    kind::Symbol,
    selectors...;
    invert=false,
    nearest=false,
    quiet=false,
    prefix::AbstractString="select",
    tag::Union{String,Nothing}=nothing,
)
    if kind == :element
        return select(domain.elems, selectors...; invert=invert, nearest=nearest, quiet=quiet, prefix=prefix, tag=tag)
    elseif kind == :face
        return select(domain.faces, selectors...; invert=invert, nearest=nearest, quiet=quiet, prefix=prefix, tag=tag)
    elseif kind == :edge
        return select(domain.edges, selectors...; invert=invert, nearest=nearest, quiet=quiet, prefix=prefix, tag=tag)
    elseif kind == :node
        return select(domain.nodes, selectors...; invert=invert, nearest=nearest, quiet=quiet, prefix=prefix, tag=tag)
    elseif kind == :ip
        length(domain.elems) > 0 && hasfield(typeof(domain.elems[1]), :ips) || return []
        ips = [ip for elem in domain.elems for ip in elem.ips]
        return select(ips, selectors...; invert=invert, nearest=nearest, quiet=quiet, prefix=prefix, tag=tag)
    else
        error("select: unknown kind $(repr(kind))")
    end
end
