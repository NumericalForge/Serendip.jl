# using Gmsh
export add_box, add_cylinder, add_sphere, add_surface_filling, add_volume, add_rectangle, add_disk
export add_point, add_line, add_circle_arc, add_circle, add_bezier, add_loop, add_wire
export add_plane_surface, add_polygon
export add_surface_loop, get_boundary, get_entities, get_points, get_curves
export get_surfaces, get_volumes, get_point, get_curve, get_surface, get_volume
export translate, rotate, extrude, revolve, mirror
export cut, fuse, intersect, fragment, fillet
export copy, get_entities, get_points, get_curves, get_surfaces, get_volumes


"""
    add_point(geo, X; size=0, embedded=false, tag="")

Add a point to the geometric model `geo` at the coordinates `X = [x, y, z]`.

# Arguments
- `geo::GeoModel`: Geometry model where the point will be added.
- `X::Vector{<:Real}`: Coordinates of the point `[x, y, z]`.
- `size::Real=0`: Optional characteristic mesh size at the point.
- `embedded::Bool=false`: Whether the point is embedded into an existing entity.
- `tag::String=""`: Optional identifier for the point.

# Returns
- `Point`: The created point entity.
"""

function add_point(geo::GeoModel, X::Vector{<:Real}; size=0, embedded::Bool=false, tag::String="")
    id = gmsh.model.occ.addPoint(X[1], X[2], X[3], size)
    gmsh.model.occ.synchronize()
    ent = Point(id, float.(X), embedded=embedded, tag=tag)
    geo.entities[(0,id)] = ent
    return ent
end


# """
#     add_line(geo, p1, p2)

# Add a straight line between points `p1` and `p2` to the geometry model `geo`.

# # Returns
# - `Edge`: Reference to the created edge.
# """

"""
    add_line(geo, p1, p2; tag="")

Add a straight line connecting points `p1` and `p2` to the geometric model `geo`.

# Arguments
- `geo::GeoModel`: Geometry model where the line will be added.
- `p1::Point`: Starting point of the line.
- `p2::Point`: Ending point of the line.
- `tag::String=""`: Optional identifier for the line.

# Returns
- `Edge`: The created line entity.
"""
function add_line(geo::GeoModel, p1::Point, p2::Point, tag::String="")
    id = gmsh.model.occ.addLine(p1.id, p2.id)
    gmsh.model.occ.synchronize()
    ent = Edge(id, "Line", [p1, p2], tag=tag)
    geo.entities[(1,id)] = ent
    return ent
end


"""
    add_circle_arc(geo, p1, p2, p3; center=true, tag="")

Add a circular arc to the geometric model `geo`, with points `p1`, `p2`, and `p3`.

# Arguments
- `geo::GeoModel`: Geometry model where the arc will be added.
- `p1::Point`: Starting point of the arc.
- `p2::Point`: Intermediate point or circle center, depending on `center`.
- `p3::Point`: Ending point of the arc.
- `center::Bool=true`: If `true`, `p2` is treated as the circle center; otherwise, it is a point on the arc.
- `tag::String=""`: Optional identifier for the arc.

# Returns
- `Edge`: The created arc entity.
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
    add_circle(geo, X, A, r; angle1=0.0, angle2=2π, tag="")

Add a circular curve to the geometric model `geo`, centered at `X` with normal vector `A` and radius `r`.

# Arguments
- `geo::GeoModel`: Geometry model where the circle will be added.
- `X::Vector{<:Real}`: Center coordinates `[x, y, z]`.
- `A::Vector{<:Real}`: Normal vector to the circle plane.
- `r::Real`: Circle radius.
- `angle1::Real=0.0`: Starting angle in radians.
- `angle2::Real=2π`: Ending angle in radians.
- `tag::String=""`: Optional user-defined identifier for the circle.

# Returns
- `Edge`: The created circle entity.
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
    add_bezier(geo, points; tag="")

Add a Bézier curve to the geometric model `geo` through the given sequence of `points`.

# Arguments
- `geo::GeoModel`: Geometry model where the curve will be added.
- `points::Vector{Point}`: Control points defining the Bézier curve.
- `tag::String=""`: Optional identifier for the curve.

# Returns
- `Edge`: The created Bézier curve entity.
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
    add_wire(geo, edges; tag="")

Add a wire (an open or closed sequence of edges) to the geometric model `geo`.

# Arguments
- `geo::GeoModel`: Geometry model where the wire will be added.
- `edges::Vector{Edge}`: Ordered list of edges forming the wire.
- `tag::String=""`: Optional identifier for the wire.

# Returns
- `Wire`: The created wire entity.
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

Create a closed curve loop from the sequence of `edges`. Used to define surface boundaries.

# Arguments
- `geo::GeoModel`: Geometry model where the loop will be added.
- `edges::Vector{Edge}`: Ordered list of edges forming a closed loop.

# Returns
- `Loop`: The created curve loop entity.
"""
function add_loop(geo::GeoModel, edges::Vector{Edge})
    ids = [ edge.id for edge in edges ]
    id = gmsh.model.occ.addCurveLoop(ids)
    gmsh.model.occ.synchronize()
    return Loop(id)
end


"""
    add_plane_surface(geo, loops...; tag="")

Create a planar surface in the geometric model `geo`, bounded by one or more `loops`.

# Arguments
- `geo::GeoModel`: Geometry model where the surface will be added.
- `loops::Loop...`: One outer loop and optionally inner loops (holes).
- `tag::String=""`: Optional identifier for the surface.

# Returns
- `Surface`: The created surface entity.
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
    add_polygon(geo, points; tag="")

Create a planar polygonal surface in the geometric model `geo` from an ordered list of `points`.

# Arguments
- `geo::GeoModel`: Geometry model where the polygon will be added.
- `points::Vector{Point}`: Ordered vertices of the polygon. Edges are created sequentially and closed automatically.
- `tag::String=""`: Optional user-defined identifier for the surface.

# Returns
- `Surface`: The created polygonal surface entity.
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
    add_disk(geo, X, A, r1, r2=0; tag="")

Create a circular or elliptical surface in the geometric model `geo`, centered at `X`, oriented along the normal vector `A`, with radii `r1` and `r2`.  
If `r2 == 0`, a circular disk is created.

# Arguments
- `geo::GeoModel`: Geometry model where the disk will be added.
- `X::Vector{<:Real}`: Center coordinates `[x, y, z]`.
- `A::Vector{<:Real}`: Normal vector defining the disk orientation.
- `r1::Real`: Major radius.
- `r2::Real=0`: Minor radius (defaults to `r1` if zero).
- `tag::String=""`: Optional user-defined identifier for the surface.

# Returns
- `Surface`: The created disk entity.
"""
function add_disk(geo::GeoModel, X::Vector{<:Real}, A::Vector{<:Real}, r1::Real, r2::Real=0; tag::String="")
    r2==0 && (r2 = r1)  # If r2 is not provided, use r1 for a full disk
    id = gmsh.model.occ.addDisk(X..., r1, r2, -1, A)
    ent =  Surface(id, tag=tag)
    geo.entities[(2,id)] = ent
    return ent
end


"""
    add_rectangle(geo, X, dx, dy; tag="")

Create a rectangular planar surface in the geometric model `geo`, starting at corner `X` with dimensions `dx` and `dy`.

# Arguments
- `geo::GeoModel`: Geometry model where the rectangle will be added.
- `X::Vector{<:Real}`: Corner coordinates `[x, y, z]`.
- `dx::Real`: Rectangle width along the x-direction.
- `dy::Real`: Rectangle height along the y-direction.
- `tag::String=""`: Optional user-defined identifier for the surface.

# Returns
- `Surface`: The created rectangular surface entity.
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
    add_surface_filling(geo, loop; tag="")

Create a filled surface bounded by the closed `loop`. Useful for non-planar or curved boundaries.

# Arguments
- `geo::GeoModel`: Geometry model where the surface will be added.
- `loop::Loop`: Boundary curve loop defining the surface perimeter.
- `tag::String=""`: Optional user-defined identifier for the surface.

# Returns
- `Surface`: The created filled surface entity.
"""
function add_surface_filling(geo::GeoModel, loop::Loop; tag::String="")
    id = gmsh.model.occ.addSurfaceFilling(loop.id)
    gmsh.model.occ.synchronize()
    ent = Surface(id, tag=tag)
    geo.entities[(2,id)] = ent
    return ent
end


"""
    add_surface_loop(surfaces)

Create a surface loop from a collection of `surfaces`. Used to define closed regions for volume creation.

# Arguments
- `surfaces::Vector{Surface}`: List of surfaces forming a closed shell.

# Returns
- `SurfaceLoop`: The created surface loop entity.
"""
function add_surface_loop(surfaces::Vector{Surface})
    ids = [ surf.id for surf in surfaces ]
    id = gmsh.model.occ.addSurfaceLoop(ids)
    return SurfaceLoop(id)
end


"""
    add_volume(geo, sloops; tag="")

Create a 3D volume in the geometric model `geo`, bounded by the given surface loops `sloops`.

# Arguments
- `geo::GeoModel`: Geometry model where the volume will be added.
- `sloops::Vector{SurfaceLoop}`: List of closed surface loops defining the volume boundary.
- `tag::String=""`: Optional user-defined identifier for the volume.

# Returns
- `Volume`: The created volume entity.
"""
function add_volume(geo::GeoModel, sloops::Vector{SurfaceLoop}; tag::String="")
    ids = [ sloop.id for sloop in sloops ]
    id = gmsh.model.occ.addVolume(ids)
    gmsh.model.occ.synchronize()
    ent = Volume(id, tag=tag)
    geo.entities[(3,id)] = ent
    return ent
end


"""
    add_box(geo, X, dx, dy, dz; tag="")

Create a box volume in the geometric model `geo`, starting at corner `X` with dimensions `dx`, `dy`, and `dz`.

# Arguments
- `geo::GeoModel`: Geometry model where the box will be added.
- `X::Vector{<:Real}`: Corner coordinates `[x, y, z]`.
- `dx::Real`: Box size along the x-direction.
- `dy::Real`: Box size along the y-direction.
- `dz::Real`: Box size along the z-direction.
- `tag::String=""`: Optional user-defined identifier for the volume.

# Returns
- `Volume`: The created box entity.
"""
function add_box(geo::GeoModel, X::Vector{<:Real}, dx::Real, dy::Real, dz::Real; tag::String="")
    id = gmsh.model.occ.add_box(X..., dx, dy, dz)
    gmsh.model.occ.synchronize()
    ent = Volume(id, tag=tag)
    geo.entities[(3,id)] = ent
    return ent
end


"""
    add_cylinder(geo, X, A, r; angle=2π, tag="")

Create a cylindrical volume in the geometric model `geo`, with base center `X`, axis vector `A`, radius `r`, and optional angular span `angle`.

# Arguments
- `geo::GeoModel`: Geometry model where the cylinder will be added.
- `X::Vector{<:Real}`: Base center coordinates `[x, y, z]`.
- `A::Vector{<:Real}`: Axis vector defining cylinder height and direction.
- `r::Real`: Cylinder radius.
- `angle::Real=2π`: Angular extent of the cylinder in radians.
- `tag::String=""`: Optional user-defined identifier for the volume.

# Returns
- `Volume`: The created cylindrical volume entity.
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
    add_sphere(geo, X, r; tag="")

Create a spherical volume in the geometric model `geo`, centered at `X` with radius `r`.

# Arguments
- `geo::GeoModel`: Geometry model where the sphere will be added.
- `X::Vector{<:Real}`: Sphere center coordinates `[x, y, z]`.
- `r::Real`: Sphere radius.
- `tag::String=""`: Optional user-defined identifier for the volume.

# Returns
- `Volume`: The created spherical volume entity.
"""
function add_sphere(geo::GeoModel, X::Vector{<:Real}, r::Real; tag::String="")
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
    kind in (:point, :edge, :curve, :face, :surface, :volume) || error("select: Unknown kind '$kind'")
    dim = kind == :point ? 0 : kind in (:edge, :curve) ? 1 : kind in (:face, :surface) ? 2 : 3

    if selector isa Vector
        dim_id = _get_entity(dim, selector)
        dim_id === nothing && return GeoEntity[]

        if kind== :point
            obj = get(geo.entities, dim_id, Point(dim_id[2], X))
        elseif kind in (:edge, :curve)
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
            push!(entities, Edge(id, "Line")) # TODO
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

"""
    translate(geo, ents, dx, dy, dz)

Translate one or more entities `ents` in the geometric model `geo` by the displacement components `(dx, dy, dz)`.

# Arguments
- `geo::GeoModel`: Geometry model containing the entities.
- `ents`: A single entity or a vector of entities to translate.
- `dx::Real`: Translation along the x-axis.
- `dy::Real`: Translation along the y-axis.
- `dz::Real`: Translation along the z-axis.

# Returns
- `Nothing`
"""
function translate(geo::GeoModel, ents, dx, dy, dz)
    ents = ents isa Vector ? ents : [ents]
    gmsh.model.occ.synchronize()
    dimids = _get_dimids_from_entities(ents)
    gmsh.model.occ.translate(dimids, dx, dy, dz)
    gmsh.model.occ.synchronize()
end


"""
    rotate(geo, ents, X, A, angle)

Rotate one or more entities `ents` in the geometric model `geo` by `angle` radians
around the axis defined by point `X` and direction `A`.

# Arguments
- `geo::GeoModel`: Geometry model containing the entities.
- `ents`: A single entity or a vector of entities to rotate.
- `X::AbstractVector{<:Real}`: Point on the rotation axis `[x, y, z]`.
- `A::AbstractVector{<:Real}`: Axis direction vector `[ax, ay, az]`.
- `angle::Real`: Rotation angle in radians.

# Returns
- `Nothing`
"""
function rotate(geo::GeoModel, ents, X, A, angle)
    ents = ents isa Vector ? ents : [ents]
    gmsh.model.occ.synchronize()
    dimids = _get_dimids_from_entities(ents)
    gmsh.model.occ.rotate(dimids, X..., A..., angle)
    gmsh.model.occ.synchronize()
end


"""
    extrude(geo, ents, A; num_elements=[], heights=[], recombine=false)

Extrude one or more entities `ents` in the geometric model `geo` along vector `A`.

# Arguments
- `geo::GeoModel`: Geometry model containing the source entities.
- `ents`: A single entity or a vector of entities to extrude.
- `A::AbstractVector{<:Real}`: Extrusion vector `[dx, dy, dz]`.
- `num_elements::Vector{Int}=[]`: Optional layer counts per step for transfinite extrusion.
- `heights::Vector{<:Real}=[]`: Optional layer heights per step. Must match `num_elements` if provided.
- `recombine::Bool=false`: If `true`, recombine into quads/hexas where possible.

# Returns
- `Vector`: Newly created entities resulting from the extrusion.
"""
function extrude(geo::GeoModel, ents, A; num_elements=[], heights=[], recombine=false)
    ents = ents isa Vector ? ents : [ents]
    gmsh.model.occ.synchronize()
    dimids = _get_dimids_from_entities(ents)
    dimids = gmsh.model.occ.extrude(dimids, A..., num_elements, heights, recombine)
    gmsh.model.occ.synchronize()
    return _get_entities_from_dimids(geo, dimids)
end


"""
    revolve(geo, ents, X, A, angle; num_elements=[], heights=[], recombine=false)

Revolve one or more entities `ents` in the geometric model `geo` by `angle` radians
around the axis defined by point `X` and direction `A`.

# Arguments
- `geo::GeoModel`: Geometry model containing the source entities.
- `ents`: A single entity or a vector of entities to revolve.
- `X::AbstractVector{<:Real}`: Point on the rotation axis `[x, y, z]`.
- `A::AbstractVector{<:Real}`: Axis direction vector `[ax, ay, az]`.
- `angle::Real`: Sweep angle in radians.
- `num_elements::Vector{Int}=[]`: Optional layer counts per step.
- `heights::Vector{<:Real}=[]`: Optional layer heights per step.
- `recombine::Bool=false`: If `true`, recombine into quads/hexas where possible.

# Returns
- `Nothing`
"""
function revolve(geo::GeoModel, ents, X, A, angle, num_elements=[], heights=[], recombine=false)
    ents = ents isa Vector ? ents : [ents]
    gmsh.model.occ.synchronize()
    dimids = _get_dimids_from_entities(ents)
    gmsh.model.occ.revolve(dimids, X..., A..., angle, num_elements, heights, recombine)
    gmsh.model.occ.synchronize()
end


"""
    mirror(ents, a, b, c, d)

Mirror one or more entities `ents` with respect to the plane `a*x + b*y + c*z + d = 0`.

# Arguments
- `ents`: A single entity or a vector of entities to mirror.
- `a::Real`: Plane coefficient for `x`.
- `b::Real`: Plane coefficient for `y`.
- `c::Real`: Plane coefficient for `z`.
- `d::Real`: Plane offset.

# Returns
- `Nothing`
"""
function mirror(ents, a, b, c, d)
    ents = ents isa Vector ? ents : [ents]
    gmsh.model.occ.synchronize()
    dimids = _get_dimids_from_entities(ents)
    gmsh.model.occ.mirror(dimids, a, b, c, d)
    gmsh.model.occ.synchronize()
end


# Boolean operations on entities

"""
    cut(geo, ents1, ents2; remove_object=false, remove_tool=true)

Boolean cut of `ents1` by `ents2` in the geometric model `geo`.  
All entities in each set must have the same topological dimension.

# Arguments
- `geo::GeoModel`: Geometry model where the operation occurs.
- `ents1`: Objects to be cut (all same dimension).
- `ents2`: Tools that cut `ents1` (all same dimension as `ents1`).
- `remove_object::Bool=false`: If `true`, remove original objects.
- `remove_tool::Bool=true`: If `true`, remove tool entities.

# Returns
- `Vector`: Resulting entities from the cut.
"""
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


"""
    fuse(geo, ents1, ents2; remove_object=true, remove_tool=true)

Boolean union of `ents1` and `ents2` in the geometric model `geo`.  
All entities in each set must have the same topological dimension.

# Arguments
- `geo::GeoModel`: Geometry model where the operation occurs.
- `ents1`: First set of entities (same dimension).
- `ents2`: Second set of entities (same dimension).
- `remove_object::Bool=true`: If `true`, remove original objects.
- `remove_tool::Bool=true`: If `true`, remove tool entities.

# Returns
- `Vector`: Resulting fused entities.
"""
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


"""
    intersect(geo, ents1, ents2; remove_object=true)

Boolean intersection of `ents1` and `ents2` in the geometric model `geo`.  
All entities in each set must have the same topological dimension.

# Arguments
- `geo::GeoModel`: Geometry model where the operation occurs.
- `ents1`: First set of entities (same dimension).
- `ents2`: Second set of entities (same dimension).
- `remove_object::Bool=true`: If `true`, remove original inputs.

# Returns
- `Vector`: Resulting intersection entities.
"""
function Base.intersect(geo::GeoModel, ents1, ents2, remove_object=true)
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


"""
    fragment(geo, ents1, ents2; remove_object=true, remove_tool=true)

Mutually fragment `ents1` and `ents2` in the geometric model `geo`, splitting them at intersections.

# Arguments
- `geo::GeoModel`: Geometry model where the operation occurs.
- `ents1`: First set of entities.
- `ents2`: Second set of entities.
- `remove_object::Bool=true`: If `true`, remove original objects.
- `remove_tool::Bool=true`: If `true`, remove tool entities.

# Returns
- `Vector`: All resulting fragmented entities.
"""
function fragment(geo::GeoModel, ents1, ents2, remove_object=true, remove_tool=true)
    ents1 = ents1 isa Vector ? ents1 : [ents1]
    ents2 = ents2 isa Vector ? ents2 : [ents2]
    gmsh.model.occ.synchronize()

    dimids1 = _get_dimids_from_entities(ents1)
    dimids2 = _get_dimids_from_entities(ents2)

    dimids, _ = gmsh.model.occ.fragment(dimids1, dimids2, -1, remove_object, remove_tool)
    gmsh.model.occ.synchronize()

    return _get_entities_from_dimids(geo, dimids)
end


"""
    fillet(geo, volumes, curves, radii; remove_volume=true)

Apply edge fillets of given `radii` to `curves` on `volumes` in the geometric model `geo`.

# Arguments
- `geo::GeoModel`: Geometry model where the operation occurs.
- `volumes`: A volume or vector of volumes to modify.
- `curves`: A curve or vector of curves to fillet.
- `radii`: A radius or vector of radii corresponding to `curves`.
- `remove_volume::Bool=true`: If `true`, remove original volumes.

# Returns
- `Vector`: Resulting modified entities after filleting.
"""
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
