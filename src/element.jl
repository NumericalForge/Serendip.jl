# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

abstract type ElementFormulation end

# abstract type ElementCache{T <: ElementFormulation} end
abstract type ElementCache end

mutable struct  Element{T}<:AbstractCell where T<:ElementFormulation
    id    ::Int
    shape ::CellShape
    role  ::Symbol # :vertex, :line, :bulk, :surface, :contact, :cohesive, :line_interface, :tip
    etype ::T       # ElementFormulation
    tag   ::String
    active::Bool
    nodes ::Vector{Node}
    ips   ::Vector{Ip}
    cmodel::Constitutive
    couplings::Vector{Element}
    cacheV::Vector{FixedSizeVector{Float64}}
    cacheM::Vector{FixedSizeMatrix{Float64}}
    cacheD::Dict{Symbol,Any}
    cache::ElementCache
    ctx  ::Context

    function Element{T}() where T<:ElementFormulation
        return new()
    end
end



# Functions that should be available in all concrete types derived from Element

"""
`elem_config_dofs(elem)`

Sets up the dofs for all nodes in `elem` according with its type.
This function can be specialized by concrete types.
"""
function elem_config_dofs(elem::Element)
    return nothing
end


"""
`elem_init(elem)`

Configures `elem` according to its type.
This function can be specialized by concrete types.
"""
function elem_init(elem::Element)
    for ip in elem.ips
        init_state(elem.cmodel, ip.state)
    end
    return nothing
end


"""
`elem_vals(elem)`

Returns a dictionary with values for the element.
Those values are intended to be constant along the element.
This function can be specialized by concrete types.
"""
function elem_vals(elem::Element)
    return Dict{Symbol, Float64}()
end


"""
`elem_recover_nodal_values(elem)`

Returns a dictionary with nodal values obtained by extrapolation
of values at ip points.
"""
function elem_recover_nodal_values(elem::Element)
    return Dict{Symbol, Float64}()
end



# Auxiliary functions for elements
# ================================

# Get the element coordinates matrix
function get_coords(elem::Element)
    nnodes = length(elem.nodes)
    ndim   = elem.ctx.ndim
    return [ elem.nodes[i].coord[j] for i in 1:nnodes, j=1:ndim]
end

# Set the quadrature points for an element
# Default implementation of bulk line and interface elements
# Other elements such as beams, shells and line_interfaces should specialize this function
function set_quadrature(elem::Element, n::Int=0; state::NamedTuple=NamedTuple())

    if !(n in keys(elem.shape.quadrature))
        alert("set_quadrature: cannot set $n integration points for shape $(elem.shape.name)")
        return
    end

    shape = elem.shape

    ipc = get_ip_coords(shape, n)
    n = size(ipc,1)

    resize!(elem.ips, n)
    for i in 1:n
        R = ipc[i].coord
        w = ipc[i].w
        elem.ips[i] = Ip(R, w)
        elem.ips[i].id = i
        elem.ips[i].state = compat_state_type(typeof(elem.cmodel), typeof(elem.etype))(elem.ctx; state...)
        elem.ips[i].owner = elem
    end

    # finding ips global coordinates
    C = get_coords(elem)

    # fix for interface elements
    if elem.role in (:cohesive, :contact)
        C = C[1:shape.npoints, : ]
    end

    # interpolation
    for ip in elem.ips
        N = shape.func(ip.R)
        ip.coord = extend!(C'*N, 3)
    end
end


function set_quadrature(elems::Array{<:Element,1}, n::Int=0)
    for elem in elems
        set_quadrature(elem, n)
    end
end


function changequadrature!(elems::Array{<:Element,1}, n::Int=0)
    set_quadrature(elems, n)
    foreach(elem_init, elems)
end


function update_material!(elem::Element, mat::Constitutive)
    typeof(elem.cmodel) == typeof(mat) || error("update_material!: The same material type should be used.")
    elem.cmodel = mat
end

"""
`changemat(elems, mat)`

Especifies the material model `mat` to be used to represent the behavior of a set of `Element` objects `elems`.
"""
function update_material!(elems::Array{<:Element,1}, mat::Constitutive)
    length(elems)==0 && notify("update_material!: Defining material model $(typeof(mat)) for an empty array of elements.")

    for elem in elems
        update_material!(elem, mat)
    end
end


# Define the state at all integration points in a collection of elements
# function set_state(elems::Array{<:Element,1}; args...)
function set_state(model, filter; args...)
    found = Set{Symbol}()
    notfound = Set{Symbol}()

    elems = select(model, :element, filter)
    length(elems)==0 && notify("set_state: Setting state to an empty array of elements.")

    for elem in elems
        fields = fieldnames(compat_state_type(typeof(elem.cmodel), typeof(elem.etype)))
        for ip in elem.ips
            for (k,v) in args
                if k in fields
                    setfield!(ip.state, k, v)
                    push!(found, k)
                else
                    push!(notfound, k)
                end
            end
        end
    end

    for k in notfound
        if k in found
            alert("set_state: field '$k' was not found at some elements while setting state values")
        else
            error("set_state: field '$k' was not found at selected elements while setting state values")
        end
    end
end


# Get all nodes from a collection of elements
# function get_nodes(elems::Array{<:Element,1})
#     nodes = Set{Node}()
#     for elem in elems
#         for node in elem.nodes
#             push!(nodes, node)
#         end
#     end
#     return [node for node in nodes]
# end

# Get all dofs from an element
function get_dofs(elem::Element)
    return Dof[ dof for node in elem.nodes for dof in node.dofs ]
end


function get_dofs(elem::Element, dofname::Symbol)
    return Dof[ dof for node in elem.nodes for dof in node.dofs if dof.name==dofname ]
end


# Get all ips from a collection of elements
function get_ips(elems::Array{<:Element,1})
    return Ip[ ip for elem in elems for ip in elem.ips ]
end


# function Base.getproperty(elems::Array{<:Element,1}, s::Symbol)
#     s == :ips   && return get_ips(elems)
#     return invoke(getproperty, Tuple{Array{<:AbstractCell,1}, Symbol}, elems, s)
# end


# General element sorting
function Base.sort!(elems::Array{<:Element,1})
    length(elems)==0 && return

    # General sorting
    sorted = sort(elems, by=elem->sum(get_coords(elem.nodes)))

    # Check type of elements
    shapes = [ elem.shape for elem in elems ]
    shape  = shapes[1]

    if all(shapes.==shape)
        if shape in (LIN2, LIN3)
            node_ids = Set(node.id for elem in sorted for node in elem.nodes)
            for elem in sorted

            end
        elseif shape in (JLIN3, JLIN4)
        end
    end

    # General sorting
    return sorted
end


"""
    get_values(elem::Element)

Returns a `DataTable` containing state variable values at all integration points of the given element.

# Arguments
- `elem::Element`: The element whose integration point values are to be retrieved.

# Returns
- `DataTable`: A table with one row per integration point and columns for each state variable.
"""
function get_values(elem::Element)
    table = DataTable()
    for ip in elem.ips
        D = state_values(elem.cmodel, ip.state)
        push!(table, D)
    end

    return table
end

