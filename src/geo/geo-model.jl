struct GPath<:GeoEntity
    path::Path
    embedded::Bool
    shape::CellShape
    tag::String
    interface_tag::String
    tip_tag::String
    tips::Symbol

    function GPath(path::Path; embedded::Bool=false, shape::CellShape=LIN3, tag::String="", interface_tag::String="", tip_tag::String="", tips=:none)
        this = new(path, embedded, shape, tag, interface_tag, tip_tag, tips)
        return this
    end
end

export add_path, add_array, add_block


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


mutable struct GeoModel
    added_entities::OrderedDict{Tuple{Int,Int},GeoEntity}
    blocks::Vector{AbstractBlock}
    gpaths::Vector{GPath}

    function GeoModel(; quiet::Bool=false)
        quiet || printstyled("Geometry model:\n", bold=true, color=:cyan)
        quiet || println("  gmsh-occ, blocks")

        gmsh.isInitialized()==1 && gmsh.finalize()
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        this = new()
        this.blocks = AbstractBlock[]
        this.added_entities = OrderedDict{Tuple{Int,Int},GeoEntity}()
        this.gpaths = GPath[]
        return this
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
"""
function add_path(geometry::GeoModel, edges::Vector{Edge}; embedded::Bool=false, shape::CellShape=LIN3, tag::String="", interface_tag::String="", tip_tag::String="", tips=:none)
    @check tips in (:none, :start, :end, :both)

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
