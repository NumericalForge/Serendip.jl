# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

import Base.maximum
import Base.minimum
import Base.sort

"""
`IpState`

Abstract type for objects to store the state at integration points.
"""
abstract type IpState
    #ctx::Context
    #other data
end


function init_state(::Constitutive, ::IpState; args...)
end


function Base.copy(src::IpState)
    T = typeof(src)
    dst = ccall(:jl_new_struct_uninit, Any, (Any,), T)
    names = fieldnames(T)
    for name in names
        val = getfield(src, name)
        if hasmethod(copy, (typeof(val),))
            setfield!(dst, name, copy(val))
        else
            setfield!(dst, name, val)
        end
    end
    return dst
end


function Base.copyto!(dst::IpState, src::IpState)
    names = fieldnames(typeof(src))
    for name in names
        val = getfield(src, name)
        if isbits(val)
            setfield!(dst, name, val)
        elseif typeof(val)<:AbstractArray
            copyto!(getfield(dst,name), val)
        else
            setfield!(dst, name, val)
            #error("copyto!(::IpState, ::IpState): unsupported field type: $(typeof(val))")
        end
    end
    return dst
end


"""
`Ip(R, w)`

Creates an `Ip` object that represents an Integration Point in finite element analyses.
`R` is a vector with the integration point local coordinates and `w` is the corresponding integration weight.
"""
mutable struct Ip
    R    ::Vec3
    w    ::Float64
    coord::Vec3
    id   ::Int
    tag  ::String
    owner::AbstractCell  # Element
    state::IpState  # Ip current state

    function Ip(R::AbstractArray{<:Float64}, w::Float64)
        this     = new(Vec3(R), w)
        this.coord = zeros(SVector{3})
        # this.coord = Vec3()
        this.tag = ""
        # this.owner = nothing
        return this
    end
end

# The functions below can be used in conjuntion with sort
# get_x(ip::Ip) = ip.coord[1]
# get_y(ip::Ip) = ip.coord[2]
# get_z(ip::Ip) = ip.coord[3]


function get_values(ip::Ip)
    coords = OrderedDict( :x => ip.coord[1], :y => ip.coord[2], :z => ip.coord[3] )
    vals   = state_values(ip.owner.cmodel, ip.state)
    return merge(coords, vals)
end


function select(
    ips::Vector{Ip},
    selectors::Union{Symbol, Symbolic, Expr, String, Vector{<:Real}, NTuple{N, Symbolic} where N}...;
    invert = false,
    nearest = true,
    tag = ""
    )

    selectors = flatten(selectors)
    selected = collect(1:length(ips))

    for selector in selectors

        if typeof(selector) == Symbol
            if selector == :all
                continue
            else
                error("select: unknown symbol selector $(repr(selector))")
            end
        elseif typeof(selector) in (Expr, Symbolic)
            fips = ips[selected]

            T = Bool[]
            for ip in fips
                x, y, z = ip.coord.x, ip.coord.y, ip.coord.z
                push!(T, evaluate(selector, x=x, y=y, z=z))
            end

            selected = selected[T]
        elseif selector isa String
            selected = [ i for i in selected if ips[i].tag==selector ]
        # elseif selector isa Vector{Int} # selector is a vector of indexes
            # selected = intersect(selected, selector)
        elseif selector isa Vector{<:Real}
            X = Vec3(selector)
            if nearest
                ip = nearest(ips[selected], X)
                i = findfirst(isequal(ip), ips)
                selected = [ i ]
            else
                R = Bool[]
                for i in selected
                    push!(R, norm(ips[i].coord-X) < 1e-8)
                end
                selected = selected[R]
            end
        end
    end

    if invert
        selected = setdiff(1:length(ips), selected)
    end

    # Set tag for selected ips
    if tag != ""
        for i in selected
            ips[i].tag = tag
        end
    end

    return ips[selected]
end

# function Base.getindex(ips::Array{<:Ip,1}, filters::NTuple; kwargs...)
#     return getindex(ips, filters...; kwargs...)
# end


# Index operator for an collection of ips
# function Base.getindex(
#     ips::Array{Ip,1},
#     filters::Union{Expr,Symbolic,String}...;
#     invert = false
#     )

#     filtered = collect(1:length(ips))

#     for filter in filters
#         if typeof(filter) in (Expr, Symbolic)
#             fips = ips[filtered]

#             T = Bool[]
#             for ip in fips
#                 x, y, z = ip.coord.x, ip.coord.y, ip.coord.z
#                 push!(T, evaluate(filter, x=x, y=y, z=z))
#             end

#             filtered = filtered[T]
#         elseif filter isa String
#             filtered = [ i for i in filtered if ips[i].tag==filter ]
#         end
#     end

#     if invert
#         filtered = setdiff(1:length(ips), filtered)
#     end

#     return ips[filtered]
# end


# # Index operator for a ip collection using expression
# function Base.getindex(ips::Array{Ip,1}, filter::Union{Expr,Symbolic})
#     R = Ip[]
#     for ip in ips
#         x, y, z = ip.coord
#         evaluate(filter, x=x, y=y, z=z) && push!(R, ip)
#     end
#     return R
# end


# function Base.getindex(ips::Array{Ip,1}, s::String)
#     return Ip[ ip for ip in ips if ip.tag==s ]
# end


function getfromcoords(ips::Array{Ip,1}, P::AbstractArray{<:Real})
    R = Ip[]
    X = Vec3(P)
    for ip in ips
        norm(X-ip.coord) < 1e-8 && push!(R, ip)
    end
    return R
end


# Get the maximum value of a given coordinate for the whole collection of ips
function maximum(ips::Array{Ip,1}, dir::Symbol)
    idx = findfirst(isequal(dir), (:x, :y, :z))
    _, idx = findmax([ip.coord[idx] for ip in ips])
    return ips[idx]
end


function minimum(ips::Array{Ip,1}, dir::Symbol)
    idx = findfirst(isequal(dir), (:x, :y, :z))
    _, idx = findmin([ip.coord[idx] for ip in ips])
    return ips[idx]
end

export nearest
function nearest(ips::Array{Ip,1}, coord)
    n = length(ips)
    D = zeros(n)
    X = Vec3(coord)

    for (i,ip) in enumerate(ips)
        D[i] = norm(X- ip.coord)
    end

    length(D) == 0 && return nothing
    return ips[sortperm(D)[1]]
end