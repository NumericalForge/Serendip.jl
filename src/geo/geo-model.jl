
function Path(edges::Vector{Edge}; closed::Union{Bool,Symbol}=:auto)
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

    return _build_path(points, cmds; closed=closed)

end


"""
    GPath(path; mode=:interface, quadratic=true, tag="", interface_tag="", tip_tag="", tips=:none, n=1)

Creates a geometric path (`GPath`) entity, which represents a curve (typically a sequence of connected edges) embedded or placed in the geometry model.
This structure is typically used to represent linear inclusions such as reinforcements, drains, etc.
The `mode` option controls how line elements are generated along the path.
If `mode` is `:interface`, interface elements will be created along the path.
If `tips` is not `:none`, tip elements will be created at the endpoints of the path in `:interface` mode.

# Arguments
- `path::Path`: The geometric path (sequence of edges and points) to be wrapped as a GPath.
- `mode::Symbol=:interface`: Path generation mode. Use `:interface`, `:embedded`, `:conforming`, or `:free`.
- `quadratic::Bool=true`: If `true`, discretize the path with quadratic 3-node line elements. If `false`, use linear 2-node line elements.
- `tag::String=""`: Optional label or name for this path (e.g. to identify or filter later).
- `interface_tag::String=""`: Optional tag for use when interface elements are generated from the path.
- `tips::Symbol=:none`: Specify witch tips are considered (:start, :end, :both, :none) for the generation of tip interface elements.
- `tip_tag::String=""`: Optional tag for tip elements if tip elements are enabled.
- `n::Int=1`: Number of free-path line elements created per path command. Only used when `mode == :free`.

Note: Original OCC edges are removed from the geometry after constructing the path.
"""
struct GPath<:GeoEntity
    path::Path
    mode::Symbol
    quadratic::Bool
    tag::String
    interface_tag::String
    tip_tag::String
    tips::Symbol
    n::Int

    function GPath(path::Path; mode::Symbol=:interface, quadratic::Bool=true, tag::String="", interface_tag::String="", tip_tag::String="", tips=:none, n::Int=1)
        @check mode in (:interface, :embedded, :conforming, :free) "'mode' must be one of :interface, :embedded, :conforming, :free"
        @check tips in (:none, :start, :end, :both)
        @check n >= 1 "'n' must be greater than or equal to 1"
        if mode != :interface
            interface_tag != "" && warn("GPath: 'interface_tag' is ignored when mode=$(repr(mode))")
            tip_tag != "" && warn("GPath: 'tip_tag' is ignored when mode=$(repr(mode))")
            tips != :none && warn("GPath: 'tips' is ignored when mode=$(repr(mode))")
        end
        mode != :free && n != 1 && error("GPath: 'n' is only supported when mode=:free")
        this = new(path, mode, quadratic, tag, interface_tag, tip_tag, tips, n)
        return this
    end
end

export add_path, add_array, add_block, set_size, set_refinement, set_transfinite_curve, set_transfinite_surface, set_recombine, set_transfinite_volume

_resolve_geo_shape(shape::Nothing) = nothing
_resolve_geo_shape(shape::Union{Symbol,CellShape}) = get_shape(shape)


function Base.copy(p::GPath)
    # copy the path
    path = copy(p.path)

    # create a new GPath with the copied path
    return GPath(path; mode=p.mode, quadratic=p.quadratic, tag=p.tag, interface_tag=p.interface_tag, tip_tag=p.tip_tag, tips=p.tips, n=p.n)
end

function move(gpath::GPath; dx::Real=0.0, dy::Real=0.0, dz::Real=0.0)
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
    transition::Float64
end


"""
    GeoModel(; size=0.0, quiet = false)
    GeoModel(filename; size=0.0, quiet = false)

Creates a new geometry model (`GeoModel`) using Gmsh's OpenCASCADE (OCC) backend
or blocks for structured meshing.
This struct serves as a container for user-defined geometry entities and geometric paths (e.g., composed by line or arc definitions).

# Arguments
- `filename::AbstractString`: Optional STEP, IGES, or BREP file imported through Gmsh OCC.
- `size::Real=0.0`: If greater than zero, sets the maximum element size for all entities in the model.
- `quiet::Bool=false`: If `true`, suppresses Gmsh initialization messages in the console.

# Fields Initialized
- `entities`: A dictionary mapping Gmsh entity identifiers `(dim, tag)` to `GeoEntity` objects explicitly added by the user.
- `blocks`: A list of meshing `Block` structures used for meshing control or volume/surface definitions.
- `gpaths`: A list of `GPath` objects used to represent path-generated line elements (e.g., for reinforcement, interfaces, or spring elements).
- `meshes`: A list of existing meshes to include in generated meshes and use as boundary constraints.

# Example
```julia
geo = GeoModel(size=0.5)       # set the maximum element size
geo = GeoModel(quiet=true)     # silent initialization
geo = GeoModel("part.step")    # import an OpenCASCADE CAD file
```
"""
mutable struct GeoModel
    entities::OrderedDict{Tuple{Int,Int},GeoEntity}
    blocks::Vector{Block}
    gpaths::Vector{GPath}
    fields::Vector{SizeField}
    meshes::Vector{AbstractDomain}

    function GeoModel(; size=0.0, min_size=0.0, max_size=0.0, quiet::Bool=false)
        quiet || printstyled("Geometry model:\n", bold=true, color=:cyan)
        quiet || println("  gmsh-occ, blocks")

        gmsh.isInitialized()==1 && gmsh.finalize()
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        if size>0.0
            max_size = 1.2*size
            min_size = 0.83*size
        end

        max_size>0.0 && gmsh.option.setNumber("Mesh.CharacteristicLengthMax", max_size)
        min_size>0.0 && gmsh.option.setNumber("Mesh.CharacteristicLengthMin", min_size)

        return new(
            OrderedDict{Tuple{Int,Int},GeoEntity}(),
            Block[],
            GPath[],
            SizeField[],
            AbstractDomain[]
        )
    end
end


const _OCC_IMPORT_EXTENSIONS = (".step", ".stp", ".iges", ".igs", ".brep", ".brp")

function _check_occ_import_extension(filename::AbstractString)
    ext = lowercase(splitext(filename)[2])
    ext in _OCC_IMPORT_EXTENSIONS ||
        error("GeoModel: unsupported OpenCASCADE CAD file extension '$ext'. Supported extensions are: $(join(_OCC_IMPORT_EXTENSIONS, ", "))")

    return nothing
end


function GeoModel(
    filename::AbstractString;
    size=0.0,
    min_size=0.0,
    max_size=0.0,
    quiet::Bool=false,
    highest_dim_only::Bool=true,
)
    isfile(filename) || error("GeoModel: file not found: $filename")
    _check_occ_import_extension(filename)

    geo = GeoModel(; size=size, min_size=min_size, max_size=max_size, quiet=quiet)
    dimids = gmsh.model.occ.importShapes(filename, highest_dim_only)
    gmsh.model.occ.synchronize()
    isempty(dimids) && error("GeoModel: no OpenCASCADE entities imported from $filename")

    for ((dim, id), entity) in zip(dimids, _get_entities_from_dimids(geo, dimids))
        geo.entities[(Int(dim), Int(id))] = entity
    end

    quiet || printstyled("  file $filename imported\n", color=:cyan)
    return geo
end


# Save gmsh geometry as step
function save(geo::GeoModel, filename::String, quiet=false)
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


"""
    add_mesh(geo::GeoModel, mesh)

Adds an existing mesh to a geometry model.

Stored meshes are included in the mesh returned by `Mesh(geo)`. When OCC
entities are meshed, the boundary vertex nodes of the stored meshes are also
used as Gmsh point constraints so adjacent generated regions can conform to the
stored mesh boundary.

The mesh is stored by reference, following the same convention as other geometry
model contents. Later mutations to the mesh object are therefore reflected when
`Mesh(geo)` is generated.
"""
function add_mesh(geo::GeoModel, mesh::AbstractDomain)
    push!(geo.meshes, mesh)
    return mesh
end


"""
    add_block(geometry, X, dx, dy, dz;
        nx=0, ny=0, nz=0, n=0,
        rx=1.0, ry=1.0, rz=1.0, r=0.0,
        quadratic=false, shape=nothing, tag="")

Adds a structured 1D, 2D, or 3D block (`Block`) to the geometry model from an
origin point `X` and edge lengths `dx`, `dy`, and `dz`.

# Arguments
- `geometry::GeoModel`: Target geometric model to which the block is appended.
- `X::Vector{<:Real}`: Block origin coordinates. Shorter vectors are padded with zeros.
- `dx::Real, dy::Real, dz::Real`: Block lengths along the x, y, and z directions.
- `nx::Int=0, ny::Int=0, nz::Int=0`: Number of divisions along each active direction.
- `n::Int=0`: Number of divisions for 1D blocks; when provided, the block is treated as 1D and `nx` is ignored.
- `rx::Real=1.0, ry::Real=1.0, rz::Real=1.0`: Grading ratios along x, y, and z.
- `r::Real=0.0`: Grading ratio for 1D blocks; when positive, it overrides `rx`.
- `quadratic::Bool=false`: If `true` and `shape` is omitted, choose a quadratic default shape (`:lin3`, `:quad8`, or `:hex20`).
- `shape`: Optional cell shape. If omitted, it is inferred from the block dimension and `quadratic`.
- `tag::String=""`: Optional tag identifier for the block.

# Behavior
- The block dimension is inferred from `dx`, `dy`, and `dz`:
  - 1D when only `dx` is nonzero, or when `n > 0`
  - 2D when `dx` and `dy` are nonzero and `dz == 0`
  - 3D when `dz != 0`
- The chosen `shape` must be compatible with the inferred dimension.

# Returns
- `Block`: The block object appended to `geometry.blocks`.

# Example
```julia
geo = GeoModel()

# Add a hexahedral block with graded mesh along z
blk = add_block(geo, [0.0, 0.0, 0.0], 2.0, 1.0, 1.0;
    nx=8, ny=4, nz=4, rz=1.3,
    shape=:hex8,
    tag="foundation")
```
"""
function add_block(geometry::GeoModel, 
    X, dx::Real, dy::Real, dz::Real;
    nx::Int=0, ny::Int=0, nz::Int=0, n::Int=0,
    rx::Real=1.0, ry::Real=1.0, rz::Real=1.0, r::Real=0.0, 
    quadratic=false,
    shape=nothing, tag="")

    bl = Block(X, float(dx), float(dy), float(dz);
        nx=nx, ny=ny, nz=nz, n=n,
        rx=rx, ry=ry, rz=rz, r=r,
        quadratic=quadratic, shape=_resolve_geo_shape(shape), tag=tag )
    # bl = Block(X1, X2; nx=nx, ny=ny, nz=nz, n=n, rx=rx, ry=ry, rz=rz, r=r, shape=shape, tag=tag)
    push!(geometry.blocks, bl)
    return bl
end


"""
    add_block(geometry, coords;
        nx=0, ny=0,
        rx=1.0, ry=1.0,
        shape=nothing, quadratic=false, tag="")

Adds a 2D structured block from explicit quadrilateral control-point
coordinates. `coords` must have either 4 corner points or 8 points including
midside points, ordered as `[1, 2, 3, 4, 5, 6, 7, 8]` around the boundary.
"""
function add_block(geometry::GeoModel,
    coords::Matrix{<:Real};
    nx::Int=0, ny::Int=0,
    rx::Real=1.0, ry::Real=1.0,
    quadratic=false,
    shape=nothing, tag="")

    npoints = size(coords, 1)
    npoints in (4, 8) || error("add_block: coordinate matrix must have 4 or 8 points for a quadrilateral block.")

    bl = Block(coords;
        nx=nx, ny=ny,
        rx=rx, ry=ry,
        quadratic=quadratic, shape=_resolve_geo_shape(shape), tag=tag)
    push!(geometry.blocks, bl)
    return bl
end


"""
    add_path(geometry, edges; mode=:interface, quadratic=true, tag="", interface_tag="", tip_tag="", tips=:none, n=1, closed=:auto)

Adds a logical path structure (`GPath`) to the geometric model from a sequence of connected `Edge` objects.
This is useful for modeling discrete and embedded 1D elements such as reinforcement bars, drains, or inclusions.

The `mode` argument defines how the path is converted into line cells when `Mesh(geo)` is built:
- `:interface`: traces the path through host cells, creates line elements along the path, and also creates `:line_interface` coupling elements between each line segment and its host cell. Optional `tips`/`tip_tag` are only meaningful in this mode.
- `:embedded`: traces the path through host cells and creates line elements coupled directly to their host cell as embedded elements (`embedded=true`), without creating interface cells.
- `:conforming`: reuses existing mesh edges that already lie on the path. It does not create new path geometry inside cells; instead it walks the existing mesh topology from one existing path endpoint to the other.
- `:free`: creates standalone line elements directly from the path geometry, independent of host-cell crossings. Existing endpoint nodes are reused when available, otherwise new nodes are created. The optional `n` keyword subdivides each path edge into `n` line elements in this mode.

# Arguments
- `geometry::GeoModel`: Target geometric model.
- `edges::Vector{Edge}`: Sequence of connected edges that define the path.
- `mode::Symbol=:interface`: Path generation mode. Use `:interface`, `:embedded`, `:conforming`, or `:free`.
- `quadratic::Bool=true`: If `true`, discretize the path with quadratic 3-node line elements. If `false`, use linear 2-node line elements.
- `tag::String=""`: Identifier tag for the path.
- `interface_tag::String=""`: Optional tag used for interface elements.
- `tip_tag::String=""`: Optional tag for used for tip elements at endpoints.
- `tips::Symbol=:none`: Specifies which endpoint elements are generated (e.g., `:start`, `:end`, `:both`). If tips is `:none`, no tips elements are created.
- `n::Int=1`: Number of line elements created per path edge in `:free` mode.
- `closed=:auto`: Path closure policy. Use `:auto` to close only when the start and end points coincide, `false` to keep the path open, or `true` to force closure with an implicit final straight segment when needed.

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
function add_path(geometry::GeoModel, edges::Vector{Edge}; mode::Symbol=:interface, quadratic::Bool=true, tag::String="", interface_tag::String="", tip_tag::String="", tips=:none, n::Int=1, closed::Union{Bool,Symbol}=:auto)
    @check tips in (:none, :start, :end, :both) "'tips' must be one of :none, :start, :end, :both"

    path = Path(edges; closed=closed)
    gpath = GPath(path; tag=tag, mode=mode, quadratic=quadratic, interface_tag=interface_tag, tip_tag=tip_tag, tips=tips, n=n)
    push!(geometry.gpaths, gpath)

    # remove source OCC curves but keep their endpoint points available for
    # later OCC operations such as rebuilding surfaces from the same vertices.
    dim_ids = [ (1, edge.id) for edge in edges ]
    gmsh.model.occ.remove(dim_ids, false)
    gmsh.model.occ.synchronize()
    foreach(dim_id -> delete!(geometry.entities, dim_id), dim_ids)

    return gpath
end


"""
    add_path(geometry, coords...; mode=:interface, quadratic=true, tag="", interface_tag="", tip_tag="", tips=:none, n=1, closed=:auto)

Adds a logical path structure (`GPath`) to the geometric model from an ordered sequence of point coordinates.
Consecutive coordinates are connected by straight line segments.

Keyword arguments follow the same rules as [`add_path(geometry, edges; ...)`](@ref), including `closed`:
- `:auto`: close only when the first and last coordinates coincide.
- `false`: keep the path open even when the endpoint coordinates coincide.
- `true`: force closure by adding an implicit final straight segment when needed.
"""
function add_path(geometry::GeoModel, coords::Vector{<:Real}...; mode::Symbol=:interface, quadratic::Bool=true, tag::String="", interface_tag::String="", tip_tag::String="", tips=:none, n::Int=1, closed::Union{Bool,Symbol}=:auto)
    @check tips in (:none, :start, :end, :both)

    path = Path(coords...; closed=closed)

    gpath = GPath(path; tag=tag, mode=mode, quadratic=quadratic, interface_tag=interface_tag, tip_tag=tip_tag, tips=tips, n=n)
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
                move(cp, dx=i*dx, dy=j*dy, dz=k*dz)
                push!(geometry.gpaths, cp)
            end
        end
    end
end


# ❱❱❱ Selection


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
    entities = _flatten(entities, GeoEntity)
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


function get_entities(geo::GeoModel, dim::Int)
    gmsh.model.occ.synchronize()
    dimids = gmsh.model.occ.getEntities(dim)
    return _get_entities_from_dimids(geo, dimids)
end


function _get_geo_dim(kind::Symbol)
    kind in (:point, :edge, :curve, :face, :surface, :volume) || error("select: Unknown kind '$kind'")
    return kind == :point ? 0 : kind in (:edge, :curve) ? 1 : kind in (:face, :surface) ? 2 : 3
end


function _get_bbox_points(dim::Int, id::Int)
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(dim, id)
    xs = unique((xmin, xmax))
    ys = unique((ymin, ymax))
    zs = unique((zmin, zmax))
    points = Vec3[Vec3(x, y, z) for x in xs for y in ys for z in zs]
    center = Vec3((xmin + xmax)/2, (ymin + ymax)/2, (zmin + zmax)/2)
    any(point -> norm(point - center) < 1e-12, points) || push!(points, center)
    return points
end


"""
    select(geo::GeoModel, kind::Symbol, selectors...; invert=false, tag=nothing)

Select OCC geometry entities from `geo` using one or more filters. Filters are
applied sequentially (logical AND), starting from all entities of the requested
kind. Tuple selectors are flattened and applied in order.

Supported entity kinds are `:point`, `:edge`, `:curve`, `:face`,
`:surface`, and `:volume`. The symbols `:edge`/`:curve` and
`:face`/`:surface` are treated as aliases.

Supported selectors:
- `:all` keeps the current selection unchanged.
- `:none` clears the selection.
- `String` matches `entity.tag`.
- `Expr` or `Symbolic` evaluates a coordinate condition using `x`, `y`, `z`
  at the entity bounding-box corners and center. The entity is kept only if the
  condition is true at all tested bounding-box points.
- `AbstractVector{<:Real}` selects by point coordinates, keeping the smallest
  candidate entity whose bounding box contains the point.

# Arguments
- `geo::GeoModel`: Geometry model containing the OCC entities.
- `kind::Symbol`: Entity kind to select (`:point`, `:edge`, `:curve`,
  `:face`, `:surface`, or `:volume`).
- `selectors`: One or more selectors applied in sequence.

# Keyword Arguments
- `invert::Bool=false`: return the complement of the final selection.
- `tag::Union{String,Nothing}=nothing`: if a string is provided, assigns it to
  all selected entities and persists them in `geo.entities`.

# Returns
- `Vector{GeoEntity}` with the selected entities.
"""
function select(
    geo::GeoModel,
    kind::Symbol,
    selectors...;
    invert::Bool = false,
    tag::Union{String, Nothing} = nothing,
    )
    dim = _get_geo_dim(kind)

    selectors = _flatten_selectors(selectors)
    length(selectors) == 1 && selectors[1] === nothing && (selectors = ())

    entities = get_entities(geo, dim)
    selected = collect(1:length(entities))

    for selector in selectors
        if selector isa Symbol
            if selector == :all
                continue
            elseif selector == :none
                selected = Int[]
            else
                error("select: unknown symbol selector $(repr(selector))")
            end
        elseif selector isa String
            selected = Int[i for i in selected if entities[i].tag == selector]
        elseif selector isa Expr || selector isa Symbolic
            selected = Int[
                i for i in selected if all(
                    evaluate(selector, x=corner.x, y=corner.y, z=corner.z)
                    for corner in _get_bbox_points(dim, entities[i].id)
                )
            ]
        elseif selector isa AbstractVector{<:Real}
            ids = [entities[i].id for i in selected]
            dim_id = _get_entity(dim, selector; ids=ids)
            if dim_id === nothing
                selected = Int[]
            else
                i = findfirst(j -> entities[j].id == dim_id[2], selected)
                selected = i === nothing ? Int[] : Int[selected[i]]
            end
        else
            error("select: unknown selector type $(typeof(selector))")
        end
    end

    if invert
        selected = setdiff(1:length(entities), selected)
    end

    selection = entities[selected]

    if tag !== nothing
        dimids = [(dim, entity.id) for entity in selection]
        selection = _ensure_entities(geo, dimids)
        for entity in selection
            entity.tag = tag
        end
    end

    return selection
end


# ❱❱❱ Transfinite

"""
    set_transfinite_curve(geo, curve, num_nodes)

Define a transfinite discretization along a curve entity `curve` in the geometric model `geo`,  
specifying the number of mesh nodes to distribute along it.

# Arguments
- `geo::GeoModel`: Geometry model containing the curve.
- `curve`: Curve entity (typically an `Edge`) to be meshed with transfinite spacing.
- `num_nodes::Int`: Number of nodes along the curve, including endpoints.

# Returns
- `Nothing`
"""
function set_transfinite_curve(geo::GeoModel, curves::AbstractVector, num_nodes::Int)
    gmsh.model.occ.synchronize()
    for curve in curves
        gmsh.model.mesh.set_transfinite_curve(curve.id, num_nodes)
    end
end

function set_transfinite_curve(geo::GeoModel, curve::GeoEntity, num_nodes::Int)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.set_transfinite_curve(curve.id, num_nodes)
end


"""
    set_transfinite_surface(geo, surfaces)

Assign a transfinite mesh distribution to surface entities `surfaces` in the geometric model `geo`.  
Used to ensure structured meshing consistent with adjoining transfinite curves.

# Arguments
- `geo::GeoModel`: Geometry model containing the surface.
- `surfaces`: Surface entities to be meshed transfinetely.

# Returns
- `Nothing`
"""
function set_transfinite_surface(geo::GeoModel, surfaces)
    gmsh.model.occ.synchronize()
    for s in surfaces
        gmsh.model.mesh.set_transfinite_surface(s.id)
    end
end


"""
    set_recombine(geo, ents)

Enable element recombination for surface entities `ents` in the geometric model `geo`.  
This converts triangular subdivisions into quadrilateral elements during meshing.

# Arguments
- `geo::GeoModel`: Geometry model containing the surface.
- `ents`: Surface entities where recombination is applied.

# Returns
- `Nothing`
"""
function set_recombine(geo::GeoModel, ents)
    gmsh.model.occ.synchronize()
    for ent in ents
        gmsh.model.mesh.set_recombine(2, ent.id)
    end
end


"""
    set_transfinite_volume(geo, volume)

Assign a transfinite mesh distribution to a volume entity `volume` in the geometric model `geo`.  
Ensures structured, block-like hexahedral meshing compatible with transfinite boundaries.

# Arguments
- `geo::GeoModel`: Geometry model containing the volume.
- `volume`: Volume entity to be meshed using transfinite interpolation.

# Returns
- `Nothing`
"""
function set_transfinite_volume(geo::GeoModel, volume)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.set_transfinite_volume(volume.id)
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
    set_refinement(geo, X, rx, ry, rz, size1, size2; transition=0.1, roundness=0.5)

Register a mesh-refinement field centered at `X`.
The field targets element size `size1` near `X`, transitioning to `size2` away from it.
The radii `rx`, `ry`, and `rz` set the extents along the x/y/z axes.

# Arguments
- `geo`::GeoModel            : geometry model.
- `X`::AbstractVector{<:Real}: 3D center `[x, y, z]`.
- `rx`::Real                 : extent along x.
- `ry`::Real                 : extent along y.
- `rz`::Real                 : extent along z.
- `size1`::Real              : inner target size.
- `size2`::Real              : outer target size.

# Keywords
- `transition`::Real (>=0) = 0.1 : transition smoothness from `size1` to `size2`
  (`0` = sharp, `1` = linear).
- `roundness`::Real ∈ [0,1] = 0.5 : shape of the refinement region
  (`0` = box-like, `1` = ellipsoidal).

# Example
```julia
set_refinement(geo, [0.5, 0.5, 0.5], 10.0, 20.0, 15.0, 0.1, 0.5;
               transition=0.5, roundness=0.8)
```
"""
function set_refinement(
    geo::GeoModel,
    X::Vector{<:Real},
    rx::Real,
    ry::Real,
    rz::Real,
    size1::Real,
    size2::Real;
    transition::Real=0.1,
    gradient=nothing,
    roundness::Real=0.5
)
    gmsh.model.occ.synchronize()

    @check length(X) == 3 "set_refinement: coordinate 'X' must be a vector of length 3"

    @check rx>0 "set_refinement: 'rx' must be positive"
    @check ry>0 "set_refinement: 'ry' must be positive"
    @check rz>0 "set_refinement: 'rz' must be positive"
    @check size1>0 "set_refinement: 'size1' must be positive"
    @check size2>0 "set_refinement: 'size2' must be positive"

    if gradient !== nothing
        gradient isa Real || error("set_refinement: deprecated keyword 'gradient' must be a Real")
        @check gradient>=0 "set_refinement: deprecated 'gradient' must be non-negative"
        (transition != 0.1 && transition != gradient) &&
            error("set_refinement: pass only one of 'transition' or deprecated 'gradient'")
        transition = gradient
    end

    @check transition>=0 "set_refinement: 'transition' must be non-negative"
    @check 0.0<=roundness<=1.0 "set_refinement: 'roundness' must be in [0, 1]"

    min_size = min(size1, gmsh.option.getNumber("Mesh.CharacteristicLengthMin"))
    max_size = max(size2, gmsh.option.getNumber("Mesh.CharacteristicLengthMax"))
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", min_size)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", max_size)

    push!(geo.fields, SizeField(X, rx, ry, rz, size1, size2, roundness, transition))

end
