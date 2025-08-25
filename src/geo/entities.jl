abstract type GeoEntity
end


mutable struct Point<:GeoEntity
    id::Int
    coord::Union{Vec3, Nothing}
    embedded::Bool
    tag::String
    function Point(id::Integer, coord::Union{Vector{Float64}, Vec3, Nothing}=nothing; embedded=false, tag::String="")
        return new(id, coord, embedded, tag)
    end
    function Point(coord::Vector{<:Real}; tag::String="")
        if length(coord)==2
            coord = [ coord; 0.0 ]
        end
        return new(-1, coord, false, tag)
    end
end

function Base.copy(p::Point)
    return Point(-1, copy(p.coord); embedded=p.embedded, tag=p.tag)
end


function get_coords(points::Vector{Point}, ndim=3)
    npoints = length(points)
    return [ points[i].coord[j] for i in 1:npoints, j in 1:ndim ]
end


struct Edge<:GeoEntity
    id::Int
    type::String
    points::Vector{Point}
    n::Int # number of elements in the edge
    tag::String

    function Edge(id::Integer, type::String, points=[]; tag=String="")
        return new(id, type, points, 0, tag)
    end
end


struct Loop<:GeoEntity
    id::Int
end


struct Wire<:GeoEntity
    id::Int
end


struct Surface<:GeoEntity
    id::Int
    transfinite::Bool
    recombine::Bool
    tag::String

    function Surface(id::Integer; transfinite::Bool=false, recombine::Bool=false, tag="")
        return new(id, transfinite, recombine, tag)
    end
end

struct SurfaceLoop<:GeoEntity
    id::Int
end

struct Volume<:GeoEntity
    id::Int
    transfinite::Bool
    tag::String

    function Volume(id::Integer; transfinite::Bool=false, tag="")
        return new(id, transfinite, tag)
    end
end