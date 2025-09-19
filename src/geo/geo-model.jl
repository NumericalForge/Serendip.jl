
function Path(edges::Vector{Edge})
    # points = Point[]
    pos = 1
    cmds = [ PathCmd(:M, [1]) ]
    points = [ edges[1].points[1] ] # start with the first point of the first edge

    for edge in edges

        # check if the first point of the edge is the same as the last point of the previous edge
        edge.points[1] == points[end] || error("Path: First point of edge $(edge.id) does not match last point of previous edge")
        append!(points, edge.points[2:end])

        # TODO: fix idxs for the case of repeated points
        if edge.type == "Bezier"
            n = length(edge.points)
            n in (3, 4) || error("Path: Bezier edge must have 3 or 4 points")
            idxs = collect(pos:pos+n-1)
            push!(cmds, PathCmd(:B, idxs))
            pos += n-1
        elseif edge.type == "Line"
            idxs = collect(pos:pos+1)
            push!(cmds, PathCmd(:L, idxs))
            pos += 1
        elseif edge.type == "CircleArcCenter"
            idxs = collect(pos:pos+2)
            push!(cmds, PathCmd(:Ac, idxs))
            pos += 2
        elseif edge.type == "CircleArcNoCenter"
            idxs = collect(pos:pos+2)
            push!(cmds, PathCmd(:A, idxs))
            pos += 2
        else
            error("Path: Unsupported edge type $(edge.type)")
        end
    end

    closed = edges[1].points[1] == points[end]
    return Path(points, cmds, closed=closed)

end

"""
    GPath(path; embedded=false, shape=LIN3, tag="", interface_tag="", tip_tag="", tips=:none)

Creates a geometric path (`GPath`) entity, which represents a curve (typically a sequence of connected edges) embedded or placed in the geometry model.
This structure is typically used to represent linear inclusions such as reinforcements, drains, etc.
If embedded is `false`, interface elements will be created along the path.
If `tips` is not `:none`, tip elements will be created at the endpoints of the path.

# Arguments
- `path::Path`: The geometric path (sequence of edges and points) to be wrapped as a GPath.
- `embedded::Bool=false`: Whether this path should be embedded in the bulk mesh during meshing (e.g. cracks, reinforcements).
- `shape::CellShape=LIN3`: Shape function used for discretizing the path (e.g., quadratic line `LIN3`, cubic, etc.).
- `tag::String=""`: Optional label or name for this path (e.g. to identify or filter later).
- `interface_tag::String=""`: Optional tag for use when interface elements are generated from the path.
- `tips::Symbol=:none`: Specify witch tips are considered (:start, :end, :both, :none) for the generation of tip interface elements.
- `tip_tag::String=""`: Optional tag for tip elements if tip elements are enabled.

Note: Original OCC edges are removed from the geometry after constructing the path.
"""
struct GPath<:GeoEntity
    path::Path
    embedded::Bool
    shape::CellShape
    tag::String
    interface_tag::String
    tip_tag::String
    tips::Symbol

    function GPath(path::Path; embedded::Bool=false, shape::CellShape=LIN3, tag::String="", interface_tag::String="", tip_tag::String="", tips=:none)
        @check tips in (:none, :start, :end, :both)
        this = new(path, embedded, shape, tag, interface_tag, tip_tag, tips)
        return this
    end
end

export add_path, add_array, add_block, set_size, set_refinement, set_transfinite_curve, set_transfinite_surface, set_recombine, set_transfinite_volume


function Base.copy(p::GPath)
    # copy the path
    path = copy(p.path)

    # create a new GPath with the copied path
    return GPath(path; embedded=p.embedded, shape=p.shape, tag=p.tag, interface_tag=p.interface_tag, tip_tag=p.tip_tag, tips=p.tips)
end

function move!(gpath::GPath; dx::Real=0.0, dy::Real=0.0, dz::Real=0.0)
    for p in gpath.path.points
        p.coord = p.coord + Vec3(dx, dy, dz)
    end
end


struct SizeField
    coord::Vec3
    rx::Float64
    ry::Float64
    rz::Float64
    size1::Float64
    size2::Float64
    roundness::Float64
    gradient::Float64
end


"""
    GeoModel(; size=0.0, quiet = false)

Creates a new geometry model (`GeoModel`) using Gmsh's OpenCASCADE (OCC) backend
or blocks for structured meshing.
This struct serves as a container for user-defined geometry entities and geometric paths (e.g., composed by line or arc definitions).

# Arguments
- `size::Real=0.0`: If greater than zero, sets the maximum element size for all entities in the model.
- `quiet::Bool=false`: If `true`, suppresses Gmsh initialization messages in the console.

# Fields Initialized
- `entities`: A dictionary mapping Gmsh entity identifiers `(dim, tag)` to `GeoEntity` objects explicitly added by the user.
- `blocks`: A list of meshing `Block` structures used for meshing control or volume/surface definitions.
- `gpaths`: A list of `GPath` objects used to represent embedded geometric paths (e.g., for reinforcement, interfaces, or spring elements).

# Example
```julia
geo = GeoModel(size=0.5)       # set the maximum element size
geo = GeoModel(quiet=true)     # silent initialization
```
"""
mutable struct GeoModel
    entities::OrderedDict{Tuple{Int,Int},GeoEntity}
    blocks::Vector{Block}
    gpaths::Vector{GPath}
    fields::Vector{SizeField}

    function GeoModel(; size=0.0,quiet::Bool=false)
        quiet || printstyled("Geometry model:\n", bold=true, color=:cyan)
        quiet || println("  gmsh-occ, blocks")

        gmsh.isInitialized()==1 && gmsh.finalize()
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        size>0.0 && gmsh.option.setNumber("Mesh.CharacteristicLengthMax", size)

        return new(
            OrderedDict{Tuple{Int,Int},GeoEntity}(),
            Block[],
            GPath[],
            SizeField[]
        )
    end
end

# Save gmsh geometry as step
function save(geometry::GeoModel, filename::String, quiet=false)
    gmsh.model.occ.synchronize()

    # Supress output
    redirect_stdout(devnull) do
        redirect_stderr(devnull) do
            gmsh.write(filename)
        end
    end
    quiet || printstyled("  file $filename written\n", color=:cyan)
    return nothing
end


function add_block(geometry::GeoModel, X1, X2;
    nx::Int=0, ny::Int=0, nz::Int=0, n::Int=0,
    rx::Real=1.0, ry::Real=1.0, rz::Real=1.0, r::Real=0.0,
    shape=nothing, tag="")

    bl = Block(X1, X2; nx=nx, ny=ny, nz=nz, n=n, rx=rx, ry=ry, rz=rz, r=r, shape=shape, tag=tag)
    push!(geometry.blocks, bl)
    return bl
end


"""
    add_path(geometry, edges; embedded=false, shape=LIN3, tag="", interface_tag="", tip_tag="", tips=:none)

Adds a logical path structure (`GPath`) to the geometric model from a sequence of connected `Edge` objects.
This is useful for modeling discrete and embedded 1D elements such as reinforcement bars, drains, or inclusions.

# Arguments
- `geometry::GeoModel`: Target geometric model.
- `edges::Vector{Edge}`: Sequence of connected edges that define the path.
- `embedded::Bool=false`: Whether the path is embedded into a solid domain. If embedded is `false`, interface elements will be created.
- `shape::CellShape=LIN3`: Finite element shape for discretization of the path (LIN2 or LIN3).
- `tag::String=""`: Identifier tag for the path.
- `interface_tag::String=""`: Optional tag used for interface elements.
- `tip_tag::String=""`: Optional tag for used for tip elements at endpoints.
- `tips::Symbol=:none`: Specifies which endpoint elements are generated (e.g., `:start`, `:end`, `:both`). If tips is `:none`, no tips elements are created.

# Returns
- `GPath`: The path structure added to the model.

Note: Original OCC edges are removed from the geometry after constructing the path.

# Example
```julia
geo = GeoModel()
point1 = add_point(geo, [0,0,0])
point2 = add_point(geo, [1,0,0])
point3 = add_point(geo, [2,0,0])
edge1 = add_line(geo, p1, p2)
edge2 = add_line(geo, p2, p3)

path = add_path(geo, [edge1, edge2]; tag="reinforcement", interface_tag="contact")
```
"""
function add_path(geometry::GeoModel, edges::Vector{Edge}; embedded::Bool=false, shape::CellShape=LIN3, tag::String="", interface_tag::String="", tip_tag::String="", tips=:none)
    @check tips in (:none, :start, :end, :both) "'tips' must be one of :none, :start, :end, :both"

    path = Path(edges)
    gpath = GPath(path; tag=tag, embedded=embedded, shape=shape, interface_tag=interface_tag, tip_tag=tip_tag, tips=tips)
    push!(geometry.gpaths, gpath)

    # remove edges from the geometry model
    dim_ids = [ (1, edge.id) for edge in edges ]

    for dim_id in dim_ids
        gmsh.model.occ.remove([dim_id], true)
        gmsh.model.occ.synchronize()
    end

    return gpath
end


function add_path(geometry::GeoModel, coords::Vector{<:Real}...; embedded::Bool=false, shape::CellShape=LIN3, tag::String="", interface_tag::String="", tip_tag::String="", tips=:none)
    @check tips in (:none, :start, :end, :both)
    embedded && interface_tag!="" && warn("add_path: 'interface_tag' is ignored when 'embedded=true'")

    path = Path(coords...)

    gpath = GPath(path; tag=tag, embedded=embedded, shape=shape, interface_tag=interface_tag, tip_tag=tip_tag, tips=tips)
    push!(geometry.gpaths, gpath)

    return gpath
end


"""
    add_array(geometry, path; nx=1, ny=1, nz=1, dx=0.0, dy=0.0, dz=0.0)

Creates a regular array of copies of a given `GPath` and adds them to the model.

# Arguments
- `geometry::GeoModel`: Target geometric model.
- `gpath::GPath`: The original path to replicate.
- `nx, ny, nz`: Number of copies in each spatial direction.
- `dx, dy, dz`: Spacing between copies in each direction.

The original path is keeped in the model as part of the array.
"""
function add_array(geometry::GeoModel, gpath::GPath; nx::Int=1, ny::Int=1, nz::Int=1, dx::Float64=0.0, dy::Float64=0.0, dz::Float64=0.0)
    for k in 0:nz-1
        for j in 0:ny-1
            for i in 0:nx-1
                i==j==k==0 && continue
                cp = copy(gpath)
                move!(cp, dx=i*dx, dy=j*dy, dz=k*dz)
                push!(geometry.gpaths, cp)
            end
        end
    end
end





# ❱❱❱ Transfinite

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


# ❱❱❱ Mesh size


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

Sets the maximum mesh size for all entities in the Gmsh-based geometry model.

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
end


# ❱❱❱ Mesh refinement


"""
    set_refinement(geo, X, rx, ry, rz, size1, size2; gradient=0.1, roundness=0.5)

Register a mesh-refinement field centered at `X`.
The field targets element size `size1` near `X`, transitioning to `size2` away from it.
The radii `rx`, `ry`, and `rz` set the extents along the x/y/z axes. If `ry` or `rz`
are omitted, they default to `rx`.

# Arguments
- `geo`::GeoModel            : geometry model.
- `X`::AbstractVector{<:Real}: 3D center `[x, y, z]`.
- `rx`::Real                 : extent along x.
- `ry`::Real                 : extent along y.
- `rz`::Real                 : extent along z.
- `size1`::Real              : inner target size.
- `size2`::Real              : outer target size.

# Keywords
- `gradient`::Real ∈ [0,1] = 0.1 : transition gradient from `size1` to `size2`
  (`0` = sharp, `1` = linear).
- `roundness`::Real  ∈ [0,1] = 0.5 : shape of the refinement region
  (`0` = box-like, `1` = ellipsoidal).

# Example
```julia
set_refinement(geo, [0.5, 0.5, 0.5], 10.0, 20.0, 15.0, 0.1, 0.5;
               gradient=0.5, roundness=0.8)
```
"""
function set_refinement(geo::GeoModel, X::Vector{<:Real}, rx::Real, ry::Real, rz::Real, size1::Real, size2::Real; gradient::Real=0.1, roundness::Real=0.5)
    gmsh.model.occ.synchronize()

    @check length(X) == 3 "set_refinement: coordinate 'X' must be a vector of length 3"

    @check rx>0 "set_refinement: 'rx' must be positive"
    @check ry>0 "set_refinement: 'ry' must be positive"
    @check rz>0 "set_refinement: 'rz' must be positive"
    @check size1>0 "set_refinement: 'size1' must be positive"
    @check size2>0 "set_refinement: 'size2' must be positive"
    @check 0.0<=gradient<=1.0 "set_refinement: 'gradient' must be in [0, 1]"
    @check 0.0<=roundness<=1.0 "set_refinement: 'roundness' must be in [0, 1]"

    push!(geo.fields, SizeField(X, rx, ry, rz, size1, size2, roundness, gradient))

end