# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl


mutable struct Node
    id     ::Int
    coord  ::Vec3
    tag    ::String
    dofs   ::Vector{Dof}
    vals   ::OrderedDict{Symbol,Float64}
    aux    ::Bool
    elems  ::Vector{AbstractCell}

    function Node()
        this         = new()
        this.id      = -1
        this.coord   = Vec3()
        this.dofs    = Dof[]
        this.vals    = OrderedDict{Symbol,Float64}()
        this.aux     = false
        this.elems   = []
        return this
    end
    
    
    function Node(x::Real, y::Real=0.0, z::Real=0.0; tag::String="", id::Int=-1)
        x = round(x, digits=8) + 0.0 # +0.0 required to drop negative bit
        y = round(y, digits=8) + 0.0
        z = round(z, digits=8) + 0.0
        
        this         = new(id, Vec3(x,y,z), tag)
        this.dofs    = Dof[]
        this.vals    = OrderedDict{Symbol,Float64}()
        this.aux     = false
        this.elems   = []
        return this
    end


    function Node(X::AbstractArray; tag::String="", id::Int=-1)
        @assert length(X) in (1,2,3)
        return Node(X...; tag=tag, id=id)
    end
end


const null_Node = Node(NaN, NaN, NaN)
@inline null(::Type{Node}) = null_Node


#Base.hash(n::Node) = hash( (round(n.coord.x, digits=8), round(n.coord.y, digits=8), round(n.coord.z, digits=8)) )
# Base.hash(n::Node) = hash( (n.coord.x, n.coord.y, n.coord.z) )
Base.hash(n::Node) = hash( (n.coord.x+1.0, n.coord.y+2.0, n.coord.z+3.0) ) # 1,2,3 aim to avoid clash in some arrays of nodes.
Base.isequal(n1::Node, n2::Node) = hash(n1)==hash(n2)


function Base.copy(node::Node)
    newnode = Node(node.coord, tag=node.tag, id=node.id)
    for dof in node.dofs
        newdof = copy(dof)
        push!(newnode.dofs, newdof)
    end
    return newnode
end


# The functions below can be used in conjuntion with sort
get_x(node::Node) = node.coord[1]
get_y(node::Node) = node.coord[2]
get_z(node::Node) = node.coord[3]


function add_dof(node::Node, name::Symbol, natname::Symbol)
    for dof in node.dofs
        dof.name == name && return nothing # dof already exists
    end

    dof = Dof(name, natname)
    push!(node.dofs, dof)

    # sets initial values
    dof.vals[name]    = 0.0
    dof.vals[natname] = 0.0
    return nothing
end


function get_dof(node::Node, key::Symbol)
    for dof in node.dofs
        dof.name == key && return dof
        dof.natname == key && return dof
    end
    return nothing
end

get_essential_keys(node::Node) = [ dof.name for dof in node.dofs ]
get_natural_keys(node::Node)  = [ dof.natname for dof in node.dofs ]


# Node collection

function get_ndim(nodes::Vector{Node})
    any( node->node.coord[3] != 0.0, nodes ) && return 3
    any( node->node.coord[2] != 0.0, nodes ) && return 2
    return 1
end


function select(
    nodes::Vector{Node},
    selectors::Union{Symbolic, Expr, String, Vector{Int}, Vector{<:Real}, NTuple{N, Symbolic} where N}...;
    invert = false,
    nearest = true,
    tag = ""
    )

    buffer = []
    for selector in selectors
        if selector isa NTuple
            append!(buffer, selector)
        else
            push!(buffer, selector)
        end
    end
    selectors = buffer

    selected = collect(1:length(nodes))

    for selector in selectors
        if typeof(selector) in (Expr, Symbolic)
            fnodes = nodes[selected]

            T = Bool[]
            for node in fnodes
                x, y, z = node.coord.x, node.coord.y, node.coord.z
                push!(T, evaluate(selector, x=x, y=y, z=z))
            end

            selected = selected[T]
        elseif selector isa String
            selected = [ i for i in selected if nodes[i].tag==selector ]
        elseif selector isa Vector{Int} # selector is a vector of indexes
            selected = intersect(selected, selector)
        elseif selector isa Vector{<:Real}
            X = Vec3(selector)
            if nearest
                node = Serendip.nearest(nodes[selected], X)
                i = findfirst(isequal(node), ips)
                selected = [ i ]
            else
                R = Bool[]
                for i in selected
                    push!(R, norm(nodes[i].coord-X) < 1e-8)
                end
                selected = selected[R]
            end
        else
            error("select: unknown selector type $(typeof(selector))")
        end
    end

    if invert
        selected = setdiff(1:length(nodes), selected)
    end

    # Set tag for selected nodes
    if tag != ""
        for i in selected
            nodes[i].tag = tag
        end
    end

    return nodes[selected]
end


function get_nodes(nodes::Vector{Node}, P::AbstractArray{<:Real})
    R = Node[]
    X = Vec3(P)
    for node in nodes
        if norm(X-node.coord) < 1e-8
            push!(R, node)
        end
    end
    return R
end


# function getfromcoords(nodes::Vector{Node}, P::AbstractArray{<:Real})
#     R = Node[]
#     X = Vec3(P)
#     for node in nodes
#         norm(X-node.coord) < 1e-8 && push!(R, node)
#     end
#     return R
# end


# Get node coordinates for a collection of nodes as a matrix
function get_coords(nodes::Vector{Node}, ndim=3)
    nnodes = length(nodes)
    return [ nodes[i].coord[j] for i in 1:nnodes, j=1:ndim]
end


# Get node values in a dictionary
function get_values(node::Node)
    coords = OrderedDict( :x => node.coord[1], :y => node.coord[2], :z => node.coord[3] )
    all_vals = [ dof.vals for dof in node.dofs ]
    return merge(coords, all_vals...)
end


# function get_values(node::Node)
#     table = DataTable()
#     dict = OrderedDict{Symbol,Real}(:id => node.id)
#     for dof in node.dofs
#         dict = merge(dict, dof.vals)
#     end
#     push!(table, dict)
#     return table
# end


function get_values(nodes::Vector{Node})
    table = DataTable()
    for node in nodes
        dict = OrderedDict{Symbol,Real}(:id => node.id)
        for dof in node.dofs
            dict = merge(dict, dof.vals)
        end
        push!(table, dict)
    end
    return table
end


function nearest(nodes::Vector{Node}, coord)
    n = length(nodes)
    D = zeros(n)
    X = Vec3(coord)

    for (i,node) in enumerate(nodes)
        D[i] = norm(X-node.coord)
    end

    return nodes[sortperm(D)[1]]
end


# function setvalue!(dof::Dof, sym_val::Pair)
#     sym, val = sym_val
#     if haskey(dof.vals, sym)
#         dof.vals[sym] = val
#     end
# end


# function setvalue!(dofs::Vector{Dof}, sym_val::Pair)
#     for dof in dofs
#         setvalue!(dof, sym_val)
#     end
# end