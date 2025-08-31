# using Gmsh
export add_box, add_cylinder, add_sphere, add_surface_filling, add_volume, add_rectangle, add_disk
export add_point, add_line, add_circle_arc, add_circle, add_bezier, add_loop, add_wire
export add_plane_surface, add_polygon
export add_surface_loop, get_boundary, get_entities, get_points, get_curves
export get_surfaces, get_volumes, get_point, get_curve, get_surface, get_volume
export translate, rotate, extrude, revolve, mirror
export cut, fuse, intersect, fragment, fillet
export set_transfinite_curve, set_transfinite_surface, set_recombine, set_transfinite_volume, set_size
export copy, get_entities, get_points, get_curves, get_surfaces, get_volumes

"""
    add_point(geo, X; size=0, tag="")

Add a point at coordinates `X = [x, y, z]` to the geometry model `geo`. Optionally specify a characteristic mesh size `size` and a tag.

# Returns
- `Point`: Reference to the created point.
"""
function add_point(geo::GeoModel, X::Vector{<:Real}; size=0, embedded::Bool=false, tag::String="")
    id = gmsh.model.occ.addPoint(X[1], X[2], X[3], size)
    gmsh.model.occ.synchronize()
    ent = Point(id, float.(X), embedded=embedded, tag=tag)
    geo.entities[(0,id)] = ent
    return ent
end


"""
    add_line(geo, p1, p2)

Add a straight line between points `p1` and `p2` to the geometry model `geo`.

# Returns
- `Edge`: Reference to the created edge.
"""
function add_line(geo::GeoModel, p1::Point, p2::Point, tag::String="")
    id = gmsh.model.occ.addLine(p1.id, p2.id)
    gmsh.model.occ.synchronize()
    ent = Edge(id, "Line", [p1, p2], tag=tag)
    geo.entities[(1,id)] = ent
    return ent
end


"""
    add_circle_arc(geo, p1, p2, p3, center=true)

Add a circle arc passing through points `p1`, `p2`, and `p3` to the geometry model `geo`.
If center is true, the middle point is the center of the circle; otherwise the circle goes through
the middle point.

# Returns
- `Edge`: Reference to the created edge.
"""
function add_circle_arc(geo::GeoModel, p1::Point, p2::Point, p3::Point; center=true, tag::String="")
    id = gmsh.model.occ.addCircleArc(p1.id, p2.id, p3.id, -1, center)
    gmsh.model.occ.synchronize()
    type = center ? "CircleArcCenter" : "CircleArcNoCenter"
    ent = Edge(id, type, [p1, p2, p3], tag=tag)
    geo.entities[(1,id)] = ent
    return ent
end


"""
    add_circle(geo, X, A, r; angle1=0.0, angle2=2π)

Add a circular curve centered at `X` with normal vector `A` and radius `r`. Define angular range with `angle1` and `angle2`.

# Returns
- `Edge`: Reference to the created edge.
"""
function add_circle(geo::GeoModel, X::Vector{<:Real}, A::Vector{<:Real}, r; angle1=0.0, angle2=2π, tag::String="")
    id = gmsh.model.occ.addCircle(X..., r, -1, angle1, angle2, A)
    gmsh.model.occ.synchronize()
    p = Point(-1, float.(X))
    ent = Edge(id, "Circle", [p], tag=tag)
    geo.entities[(1,id)] = ent
    return ent
end


"""
    add_bezier(geo, points)

Add a Bezier curve through the given sequence of `points`.

# Returns
- `Edge`: Reference to the created edge.
"""
function add_bezier(geo::GeoModel, points::Vector{Point}; tag::String="")
    ids = [ p.id for p in points ]
    id = gmsh.model.occ.addBezier(ids)
    gmsh.model.occ.synchronize()
    ent = Edge(id, "Bezier", points, tag=tag)
    geo.entities[(1,id)] = ent
    return ent
end


"""
    add_wire(geo, edges)

Add a wire (open or closed sequence of edges) to the geometry model `geo`.

# Returns
- `Wire`: Reference to the created wire.
"""
function add_wire(geo::GeoModel, edges::Vector{Edge}; tag::String="")
    ids = [ edge.id for edge in edges ]
    id = gmsh.model.occ.addWire(ids)
    gmsh.model.occ.synchronize()
    ent =  Wire(id, tag=tag)
    geo.entities[(1,id)] = ent
    return ent
end


"""
    add_loop(geo, edges)

Create a closed curve loop from the sequence of `edges`. Used for defining bounded surfaces.

# Returns
- `Loop`: Reference to the created curve loop.
"""
function add_loop(geo::GeoModel, edges::Vector{Edge})
    ids = [ edge.id for edge in edges ]
    id = gmsh.model.occ.addCurveLoop(ids)
    gmsh.model.occ.synchronize()
    return Loop(id)
end


"""
    add_plane_surface(geo, loops...)

Create a planar surface bounded by one or more `loops`.

# Returns
- `Surface`: Reference to the created surface.
"""
function add_plane_surface(geo::GeoModel, loops::Loop...; tag::String="")
    ids = [ loop.id for loop in loops ]
    id = gmsh.model.occ.addPlaneSurface(ids)
    gmsh.model.occ.synchronize()
    ent =  Surface(id, tag=tag)
    geo.entities[(2,id)] = ent
    return ent
end


"""
    add_polygon(geo, points::Vector{Point}; tag::String="")

Create a planar polygonal surface from an ordered list of `points`.

# Arguments
- `geo::GeoModel`: Geometry model where the polygon will be added.
- `points::Vector{Point}`: Ordered vertices of the polygon. Edges are created sequentially and closed automatically.
- `tag::String`: Optional label for identifying the created surface.

# Returns
- `Surface`: Reference to the created surface.
"""
function add_polygon(geo, points::Vector{Point}; tag::String="")
    edges = Edge[]
    np = length(points)
    for i in 1:np
        p1 = points[i]
        p2 = points[ i == np ? 1 : i+1 ]
        e = add_line(geo, p1, p2)
        push!(edges, e)
    end
    loop = add_loop(geo, edges)
    return add_plane_surface(geo, loop; tag=tag)
end


"""
    add_disk(geo, X, A, r1, r2=0)

Create a disk or elliptical surface centered at `X`, oriented along normal vector `A`, with radii `r1` and `r2`. If `r2=0`, creates a circular disk.

# Returns
- `Surface`: Reference to the created surface.
"""
function add_disk(geo::GeoModel, X::Vector{<:Real}, A::Vector{<:Real}, r1::Real, r2::Real=0; tag::String="")
    r2==0 && (r2 = r1)  # If r2 is not provided, use r1 for a full disk
    id = gmsh.model.occ.addDisk(X..., r1, r2, -1, A)
    ent =  Surface(id, tag=tag)
    geo.entities[(2,id)] = ent
    return ent
end


"""
    add_rectangle(geo, X, dx, dy)

Add a rectangular surface starting at corner `X` with dimensions `dx` and `dy`.

# Returns
- `Surface`: Reference to the created surface.
"""
function add_rectangle(geo::GeoModel, X::Vector{<:Real}, dx::Real, dy::Real; tag::String="")
    Y = Vec3(X)
    id = gmsh.model.occ.addRectangle(Y..., dx, dy)
    gmsh.model.occ.synchronize()
    ent = Surface(id, tag=tag)
    geo.entities[(2,id)] = ent
    return ent
end


"""
    add_surface_filling(loop)

Create a filled surface bounded by the closed `loop`. Used for non-planar boundaries.

# Returns
- `Surface`: Reference to the created surface.
"""
function add_surface_filling(loop::Loop)
    id = gmsh.model.occ.addSurfaceFilling(loop.id)
    gmsh.model.occ.synchronize()
    ent = Surface(id, tag=tag)
    geo.entities[(2,id)] = ent
    return ent
end


"""
    add_surface_loop(surfaces)

Create a surface loop from a collection of `surfaces`. Used for defining closed volumes.

# Returns
- `SurfaceLoop`: Reference to the created surface loop.
"""
function add_surface_loop(surfaces::Vector{Surface})
    ids = [ surf.id for surf in surfaces ]
    id = gmsh.model.occ.addSurfaceLoop(ids)
    return SurfaceLoop(id)
end


"""
    add_volume(sloops)

Create a 3D volume bounded by the given `sloops` (surface loops).

# Returns
- `Volume`: Reference to the created volume.
"""
function add_volume(sloops::Vector{SurfaceLoop})
    ids = [ sloop.id for sloop in sloops ]
    id = gmsh.model.occ.addVolume(ids)
    gmsh.model.occ.synchronize()
    ent = Volume(id, tag=tag)
    geo.entities[(3,id)] = ent
    return ent
end


"""
    add_box(geo, X, dx, dy, dz)

Add a box starting at corner `X` with dimensions `dx`, `dy`, and `dz`.

# Returns
- `Volume`: Reference to the created volume.
"""
function add_box(geo::GeoModel, X::Vector{<:Real}, dx::Real, dy::Real, dz::Real; tag::String="")
    id = gmsh.model.occ.add_box(X..., dx, dy, dz)
    gmsh.model.occ.synchronize()
    ent = Volume(id, tag=tag)
    geo.entities[(3,id)] = ent
    return ent
end


"""
    add_cylinder(geo, X, A, r, angle=2π)

Create a cylindrical volume with base center `X`, axis vector `A`, radius `r`, and optional angular span `angle`.

# Returns
- `Volume`: Reference to the created volume.
"""
function add_cylinder(geo::GeoModel, X::Vector{<:Real}, A::Vector{<:Real}, r, angle=2π; tag::String="")
    # X = [x, y, z] center of the base
    # A = [a, b, c] direction vector of the cylinder axis
    id = gmsh.model.occ.add_cylinder(X..., A..., r, -1, angle)
    gmsh.model.occ.synchronize()
    ent = Volume(id, tag=tag)
    geo.entities[(3,id)] = ent
    return ent
end


"""
    add_sphere(geo, X, r)

Create a spherical volume centered at `X` with radius `r`.

# Returns
- `Volume`: Reference to the created volume.
"""
function add_sphere(geo::GeoModel, X::Vector{<:Real}, r::Real; tag::String="")
    # X = [x, y, z] center of the sphere
    id = gmsh.model.occ.addSphere(X..., r)
    gmsh.model.occ.synchronize()
    ent = Volume(id, tag=tag)
    geo.entities[(3,id)] = ent
    return ent
end


# Query functions

function get_boundary(geo::GeoModel, entities; combined=true, oriented=false, recursive=false)
    gmsh.model.occ.synchronize()
    dimids = _get_dimids_from_entities(entities)
    boundary = gmsh.model.get_boundary(dimids, combined, oriented, recursive)

    # Convert to entities
    return _get_entities_from_dimids(geo, boundary)
end


function _get_entity(dim, X; tol=1e-6)
    gmsh.model.occ.synchronize()
    x, y, z = X
    ent_id = nothing
    min_vol = Inf

    for (_, id) in gmsh.model.occ.getEntities(dim)

        # Get bounding box of the entity
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(dim, id)

        # Expand bounding box with tolerance
        xmin -= tol; ymin -= tol; zmin -= tol
        xmax += tol; ymax += tol; zmax += tol

        # Skip if point not in bounding box
        (xmin ≤ x ≤ xmax && ymin ≤ y ≤ ymax && zmin ≤ z ≤ zmax) || continue

        # Compute vol
        vol = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)

        if vol < min_vol
            min_vol = vol
            ent_id = id
        end
    end

    ent_id === nothing && return nothing
    return (dim, ent_id)

end


function get_point(geo::GeoModel, X::Vector{<:Real})
    ent = _get_entity(0, X)
    ent === nothing && return nothing
    return get(geo.entities, ent, Point(ent[2], X))
end


function get_curve(geo::GeoModel, X::Vector{<:Real})
    ent = _get_entity(1, X)
    ent === nothing && return nothing
    type = gmsh.model.getType(1, ent[2])
    return get(geo.entities, ent, Edge(ent[2], type))
end


function get_surface(geo::GeoModel, X::Vector{<:Real})
    ent = _get_entity(2, X)
    ent === nothing && return nothing
    return get(geo.entities, ent, Surface(ent[2]))
end


function get_volume(geo::GeoModel, X::Vector{<:Real})
    ent = _get_entity(3, X)
    ent === nothing && return nothing
    return get(geo.entities, ent, Volume(ent[2]))
end

# function select(geo::GeoModel, kind::Symbol; tag::String="")
#     ents = []
#     if kind == :point
#         ents = get_entities(geo, 0)
#     elseif kind in (:line, :curve)
#         ents = get_entities(geo, 1)
#     elseif kind in (:face, :surface)
#         ents = get_entities(geo, 2)
#     elseif kind == :volume
#         ents = get_entities(geo, 3)
#     else
#         error("select: Unknown kind '$kind'")
#     end

#     if tag != ""
#         ents = [ e for e in ents if e.tag == tag ]
#     end

#     return ents
# end


function select(geo::GeoModel, kind::Symbol, selector::Union{String,Vector{<:Real}}; tag::String="")
    kind in (:point, :line, :curve, :face, :surface, :volume) || error("select: Unknown kind '$kind'")
    dim = kind == :point ? 0 : kind in (:line, :curve) ? 1 : kind in (:face, :surface) ? 2 : 3

    if selector isa Vector
        dim_id = _get_entity(dim, selector)
        dim_id === nothing && return GeoEntity[]

        if kind== :point
            obj = get(geo.entities, dim_id, Point(dim_id[2], X))
        elseif kind in (:line, :curve)
            type = gmsh.model.getType(1, dim_id[2])
            obj = get(geo.entities, dim_id, Edge(dim_id[2], type))
        elseif kind in (:face, :surface)
            obj = get(geo.entities, dim_id, Surface(dim_id[2]))
        elseif kind == :volume
            obj = get(geo.entities, dim_id, Volume(dim_id[2]))
        end
        tag!= "" && (obj.tag = tag)
        geo.entities[dim_id] = obj
        return obj
    else # string tag
        ents = GeoEntity[]
        for (dim_id, obj) in geo.entities
            dim_id[1] == dim && obj.tag == selector && push!(ents, obj)
        end
        return ents
    end
end


# function select(geo::GeoModel, kind::Symbol, selector::String; tag::String="")
# end



# function get_surface2(geo::GeoModel, X::Vector{<:Real})
#     gmsh.model.occ.synchronize()
#     @show X

#     for dim_id in gmsh.model.occ.getEntities(2)
#         @show dim_id
#         U = gmsh.model.getParametrization(dim_id..., X)
#         @show U

#         if gmsh.model.isInside(dim_id..., U, true)==1
#             return dim_id
#         end
#     end
#     return nothing
# end



function get_entities(dim::Int)
    gmsh.model.occ.synchronize()
    dimids = gmsh.model.occ.getEntities(dim)
    return _get_entities_from_dimids(geo, dimids)
end


function get_points(geo, ents=[])
    ents = ents isa Vector ? ents : [ents]

    length(ents) == 0 && return get_entities(0)

    return get_boundary(geo, ents; combined=false, oriented=false, recursive=true)
end


function get_curves(ents=[])
    length(ents) == 0 && return get_entities(1)

    all_edges = []
    for e in ents
        if e isa Volume
            bry = get_boundary([e]; combined=false, oriented=false, recursive=false)
            edges = get_curves(bry)
        elseif e isa Surface
            bry = get_boundary([e]; combined=false, oriented=false, recursive=false)
            edges = [ x for x in bry if x isa Edge ]
        elseif e isa Edge
            edges = [e]
        end
        append!(all_edges, edges)
    end

    return all_edges
end


function get_surfaces(ents=[])
    length(ents) == 0 && return get_entities(2)

    all_surfaces = []
    for e in ents
        if e isa Volume
            bry = get_boundary([e]; combined=false, oriented=false, recursive=false)
            surfaces = [ x for x in bry if x isa Surface ]
        elseif e isa Surface
            surfaces = [e]
        end
        append!(all_surfaces, surfaces)
    end

    return all_surfaces
end


function get_volumes(ents=[])
    length(ents) == 0 && return get_entities(3)

    return [ e for e in ents if e isa Volume ]
end


function _get_entities_from_dimids(geo::GeoModel, dimids)
    entities = []
    for (d, id) in dimids
        ent = get(geo.entities, (d, id), nothing)
        if ent !== nothing
            push!(entities, ent)
            continue
        end

        if d == 3
            push!(entities, Volume(id))
        elseif d == 2
            push!(entities, Surface(id))
        elseif d == 1
            push!(entities, Edge(id))
        elseif d == 0
            push!(entities, Point(id))
        end
    end
    return entities
end


function _get_dimids_from_entities(entities)
    dimids = []
    for e in entities
        if e isa Volume
            push!(dimids, (3, e.id))
        elseif e isa Surface
            push!(dimids, (2, e.id))
        elseif e isa Edge
            push!(dimids, (1, e.id))
        elseif e isa Point
            push!(dimids, (0, e.id))
        end
    end
    return dimids
end


"""
    Base.copy(geo, ents::Vector{<:GeoEntity})

Copy a list of geometric entities in the geometry model `geo`. Synchronizes OCC before and after copying.

# Returns
- `Vector{GeoEntity}`: The copied entities.
"""
function Base.copy(geo::GeoModel, ents::Vector{<:GeoEntity})
    ents = ents isa Vector ? ents : [ents]
    gmsh.model.occ.synchronize()
    dimids = _get_dimids_from_entities(ents)
    dimids = gmsh.model.occ.copy(dimids)
    gmsh.model.occ.synchronize()
    return _get_entities_from_dimids(geo, dimids)
end


"""
    Base.copy(geo, ent::GeoEntity)

Copy a single geometric entity in the geometry model `geo`.

# Returns
- `GeoEntity`: The copied entity.
"""
function Base.copy(geo::GeoModel, ents::GeoEntity)
    ents = [ents]
    return copy(geo, ents)[1]
end


# Transformations

function translate(geo::GeoModel, ents, A)
    ents = ents isa Vector ? ents : [ents]
    gmsh.model.occ.synchronize()
    dimids = _get_dimids_from_entities(ents)
    gmsh.model.occ.translate(dimids, A...)
    gmsh.model.occ.synchronize()
end


function rotate(geo::GeoModel, ents, X, A, angle)
    ents = ents isa Vector ? ents : [ents]
    gmsh.model.occ.synchronize()
    dimids = _get_dimids_from_entities(ents)
    gmsh.model.occ.rotate(dimids, X..., A..., angle)
    gmsh.model.occ.synchronize()
end


function extrude(geo::GeoModel, ents, A; num_elements=[], heights=[], recombine=false)
    ents = ents isa Vector ? ents : [ents]
    gmsh.model.occ.synchronize()
    dimids = _get_dimids_from_entities(ents)
    dimids = gmsh.model.occ.extrude(dimids, A..., num_elements, heights, recombine)
    gmsh.model.occ.synchronize()
    return _get_entities_from_dimids(geo, dimids)
end


function revolve(geo::GeoModel, ents, X, A, angle, num_elements=[], heights=[], recombine=false)
    ents = ents isa Vector ? ents : [ents]
    gmsh.model.occ.synchronize()
    dimids = _get_dimids_from_entities(ents)
    gmsh.model.occ.revolve(dimids, X..., A..., angle, num_elements, heights, recombine)
    gmsh.model.occ.synchronize()
end


function mirror(ents, a, b, c, d)
    ents = ents isa Vector ? ents : [ents]
    gmsh.model.occ.synchronize()
    dimids = _get_dimids_from_entities(ents)
    gmsh.model.occ.mirror(dimids, a, b, c, d)
    gmsh.model.occ.synchronize()
end


# Boolean operations on entities

function cut(geo::GeoModel, ents1, ents2; remove_object=false, remove_tool=true)
    ents1 = ents1 isa Vector ? ents1 : [ents1]
    ents2 = ents2 isa Vector ? ents2 : [ents2]
    gmsh.model.occ.synchronize()

    dim1 = ents1[1] isa Volume ? 3 : ents1[1] isa Surface ? 2 : 1
    dim2 = ents2[1] isa Volume ? 3 : ents2[1] isa Surface ? 2 : 1

    if dim1 == dim2
        ids1 = [ (dim1, e.id) for e in ents1 ]
        ids2 = [ (dim1, e.id) for e in ents2 ]
        dimids, _ = gmsh.model.occ.cut(ids1, ids2, remove_object, remove_tool)
        gmsh.model.occ.synchronize()

        return _get_entities_from_dimids(geo, dimids)
    else
        error("cut: Entities must have the same dimension")
    end
end


function fuse(geo::GeoModel, ents1, ents2, remove_object=true, remove_tool=true)
    ents1 = ents1 isa Vector ? ents1 : [ents1]
    ents2 = ents2 isa Vector ? ents2 : [ents2]
    gmsh.model.occ.synchronize()

    dim1 = ents1[1] isa Volume ? 3 : ents1[1] isa Surface ? 2 : 1
    dim2 = ents2[1] isa Volume ? 3 : ents2[1] isa Surface ? 2 : 1

    if dim1 == dim2
        ids1 = [ (dim1, e.id) for e in ents1 ]
        ids2 = [ (dim2, e.id) for e in ents2 ]
        dimids, _ = gmsh.model.occ.fuse(ids1, ids2, -1, remove_object, remove_tool)
        gmsh.model.occ.synchronize()

        return _get_entities_from_dimids(geo, dimids)

    else
        error("fuse: Entities must have the same dimension")
    end
end


function intersect(geo::GeoModel, ents1, ents2, remove_object=true)
    ents1 = ents1 isa Vector ? ents1 : [ents1]
    ents2 = ents2 isa Vector ? ents2 : [ents2]
    gmsh.model.occ.synchronize()

    dim1 = ents1[1] isa Volume ? 3 : ents1[1] isa Surface ? 2 : 1
    dim2 = ents2[1] isa Volume ? 3 : ents2[1] isa Surface ? 2 : 1

    if dim1 == dim2
        ids1 = [ (dim1, e.id) for e in ents1 ]
        ids2 = [ (dim2, e.id) for e in ents2 ]
        dimids, _ = gmsh.model.occ.intersect(ids1, ids2, -1, remove_object)
        gmsh.model.occ.synchronize()

        return _get_entities_from_dimids(geo, dimids)

    else
        error("intersect: Entities must have the same dimension")
    end
end


function fragment(geo::GeoModel, ents1, ents2, remove_object=false, remove_tool=true)
    ents1 = ents1 isa Vector ? ents1 : [ents1]
    ents2 = ents2 isa Vector ? ents2 : [ents2]
    gmsh.model.occ.synchronize()

    dimids1 = _get_dimids_from_entities(ents1)
    dimids2 = _get_dimids_from_entities(ents2)

    dimids, _ = gmsh.model.occ.fragment(dimids1, dimids2, -1, remove_object, remove_tool)
    gmsh.model.occ.synchronize()

    return _get_entities_from_dimids(geo, dimids)
end


function fillet(geo::GeoModel, volumes, curves, radii, remove_volume=true)
    volumes = volumes isa Vector ? volumes : [volumes]
    curves = curves isa Vector ? curves : [curves]
    radii = radii isa Vector ? radii : [radii]
    gmsh.model.occ.synchronize()

    dimids1 = _get_dimids_from_entities(volumes)
    dimids2 = _get_dimids_from_entities(curves)
    gmsh.model.occ.fillet(dimids1, dimids2, radii, remove_volume)
    gmsh.model.occ.synchronize()

    return _get_entities_from_dimids(geo, dimids)
end

# transfinite

function set_transfinite_curve(geo::GeoModel, ent, num_nodes)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.set_transfinite_curve(ent.id, num_nodes)
end


function set_transfinite_surface(geo::GeoModel, ent)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.set_transfinite_surface(ent.id)
end


function set_recombine(geo::GeoModel, ent)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.set_recombine(2, ent.id)
end


function set_transfinite_volume(ent)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.set_transfinite_volume(ent.id)
end


# Mesh size

# function set_size(geo::GeoModel, point::Point, size::Real)
#     set_size(geo, [point], size)
# end

# function set_size(geo::GeoModel, points::Vector{Point}, size::Real)
#     gmsh.model.occ.synchronize()
#     dim_ids = [ (0, ent.id) for ent in points]
#     gmsh.model.mesh.set_size(dim_ids, size)
# end

"""
    set_size(geo, kind, tag, size)

Assigns a local mesh size to all points on the boundary of entities selected by type (`kind`) and tag within a `GeoModel`.

This function selects entities of a specific geometric type (`:point`, `:edge`, `:surface`, `:volume`) and tag, and sets the mesh size at the corresponding points.

# Arguments
- `geo::GeoModel`: The geometry model where the mesh size is applied.
- `kind::Symbol`: Type of geometric entity to target (`:point, ``:edge`, `:face`, `:volume`).
- `tag::String`: Tag string used to select the target entities.
- `size::Real`: Target mesh size to assign at the boundary points of the selected entities.

# Example
```julia
set_size(geo, :surface, "foundation_zone", 0.05)
```
"""
function set_size(geo::GeoModel, kind::Symbol, tag::String, size::Real)
    gmsh.model.occ.synchronize()
    ents = select(geo, kind, tag)
    dim_ids = _get_dimids_from_entities(ents)
    pts_dim_ids = gmsh.model.get_boundary(dim_ids, true, false, true)

    length(pts_dim_ids) == 0 && return
    gmsh.model.mesh.set_size(pts_dim_ids, size)
end


"""
    set_size(geo, target, size)

Assigns a target mesh size to the point entities on the boundary of the given `target` geometry entity or entities.
If a single `GeoEntity` is passed, it is internally wrapped in a vector.

# Arguments
- `geo::GeoModel`: The geometry model where the mesh size will be assigned.
- `target::Union{GeoEntity,Vector{<:GeoEntity}`: A single `GeoEntity` or a vector of them.
- `size::Real`: The target mesh size to assign to the boundary points of the given entities.

# Example
```julia
set_size(geo, some_line_entity, 0.05)
set_size(geo, [curve1, surface1], 0.1)
```
"""
function set_size(geo::GeoModel, target::Union{GeoEntity, Vector{<:GeoEntity}}, size::Real)
    if target isa GeoEntity
        set_size(geo, [target], size)
    else
        ents = target
        dim_ids = _get_dimids_from_entities(ents)
        pts_dim_ids = gmsh.model.get_boundary(dim_ids, true, false, true)

        length(pts_dim_ids) == 0 && return
        gmsh.model.mesh.set_size(pts_dim_ids, size)
    end
end


"""
    set_size(geo::GeoModel, size::Real)

Sets the global target mesh size for all entities in the Gmsh-based geometry model.
It enforces a uniform mesh size across the model.

# Arguments
- `geo::GeoModel`: The geometry model whose mesh settings will be modified.
- `size::Real`: The desired global mesh size.

# Example
```julia
set_size(geo, 0.05)  # Set global mesh size to 0.05 units
```
"""
function set_size(geo::GeoModel, size::Real)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", size)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", size)
end