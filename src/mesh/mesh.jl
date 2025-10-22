# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl


mutable struct Mesh<:AbstractDomain
    nodes::Vector{Node}
    elems::Vector{Cell}
    faces::Vector{Cell}
    edges::Vector{Cell}
    node_data::OrderedDict{String,Array}
    elem_data ::OrderedDict{String,Array}
    ctx::MeshContext

    _pointdict::Dict{UInt64,Node}
    _elempartition::ElemPartition

    function Mesh(ndim::Integer=0)
        this = new()
        this.nodes = []
        this.elems  = []
        this.faces  = []
        this.edges  = []
        this.node_data = OrderedDict()
        this.elem_data = OrderedDict()
        this.ctx = MeshContext(ndim)
        this._pointdict = Dict{UInt64, Node}()
        this._elempartition = ElemPartition()
        return this
    end
end


function get_node(nodes::Dict{UInt64,Node}, C::AbstractArray{<:Real})
    hs = hash(Node(C))
    return get(nodes, hs, nothing)
end


function Base.copy(mesh::AbstractDomain)
    ndim = mesh.ctx.ndim
    newmesh = Mesh(ndim)
    newmesh.nodes = copy.(mesh.nodes)

    for elem in mesh.elems
        idxs = [ node.id for node in elem.nodes ]
        newelemnodes = newmesh.nodes[idxs]
        newelem = Cell(elem.shape, elem.role, newelemnodes, tag=elem.tag, id=elem.id, active=elem.active)
        push!(newmesh.elems, newelem)
    end

    compute_facets(newmesh)

    newmesh._pointdict = Dict( hash(node) => node for node in newmesh.nodes )
    newmesh.node_data = copy(mesh.node_data)
    newmesh.elem_data = copy(mesh.elem_data)

    # fixing references for linked elements #todo check
    for elem in newmesh.elems
        if length(elem.couplings)>0
            idxs = [ e.id for e in elem.couplings ]
            elem.couplings = newmesh.elems[idxs]
        end
    end

    return newmesh
end


function get_outer_facets(cells::Array{<:AbstractCell,1})
    face_d = OrderedDict{UInt64, Cell}()

    # Get only unique faces. If dup, original and dup are deleted
    for cell in cells
        for face in get_facets(cell) 
            hs = hash(face)
            if haskey(face_d, hs)
                delete!(face_d, hs)
            else
                face_d[hs] = face
            end
        end
    end

    return CellFace[ face for face in values(face_d) ]
end


function get_outer_facets_by_id(cells::Array{<:AbstractCell,1})
    # face_d = OrderedDict{UInt64, Cell}()
    face_d = Dict()
    hash1(edge) = sort([ n.id for n in edge.nodes ])

    # Get only unique faces. If dup, original and dup are deleted
    for cell in cells
        for face in get_facets(cell)
            hs = hash1(face)
            if haskey(face_d, hs)
                delete!(face_d, hs)
            else
                face_d[hs] = face
            end
        end
    end

    return collect(values(face_d))
end


function get_edges(surf_cells::Array{<:AbstractCell,1})
    edges_dict = Dict{UInt64, Cell}()

    # Get only unique edges
    for cell in surf_cells
        for edge in get_edges(cell)
            edge.owner = cell.owner
            edges_dict[hash(edge)] = edge
        end
    end

    return collect(values(edges_dict))
end


# Return a list of neighbors for each cell
function get_neighbors(cells::Vector{Cell})
    faces_dict = Dict{UInt64, Cell}()
    neighbors = [ Cell[] for i in 1:length(cells) ]

    # Get cell faces. If dup, original and dup are deleted but neigh info saved
    for cell in cells
        for face in get_facets(cell)
            hs = hash(face)
            other = get(faces_dict, hs, nothing)
            if other === nothing
                faces_dict[hs] = face
            else
                push!(neighbors[face.owner.id], other.owner)
                push!(neighbors[other.owner.id], face.owner)
                delete!(faces_dict, hs)
            end
        end
    end

    return neighbors

end

# Return the cell patch for each point
function get_patches(mesh::Mesh)
    patches = [ Cell[] for i in 1:length(mesh.nodes) ]
    for cell in mesh.elems
        for pt in cell.nodes
            push!(patches[pt.id], cell)
        end
    end

    return mesh.nodes, patches
end


"""
    sort_mesh(mesh; reverse=true) -> Mesh

Reorder the nodes and elements of a finite element `mesh` to improve 
sparsity structure and cache locality in FEM analysis.
Nodes are renumbered using the **Reverse Cuthill–McKee (RCM) algorithm.
Elements are subsequently sorted by the minimum node id they contain.

# Arguments
- `mesh::Mesh`: The mesh to be reordered.
- `reverse::Bool=true`: If true, apply Reverse Cuthill–McKee; if false, 
  use standard Cuthill–McKee.

# Returns
The same `mesh` object, modified in place, with reordered nodes and elements.
"""
function sort_mesh(mesh::Mesh; reverse::Bool=true)
    nnodes = length(mesh.nodes)

    # Build node–node adjacency from cell connectivities
    adjacents = [Int[] for _ in 1:nnodes]
    for cell in mesh.elems
        ids = getfield.(cell.nodes, :id)
        ncellnodes = length(ids)
        for i in 1:ncellnodes-1, j in i+1:ncellnodes
            u, v = ids[i], ids[j]
            push!(adjacents[u], v)
            push!(adjacents[v], u)
        end
    end

    unique!.(adjacents)
    deg = length.(adjacents)

    # RCM: BFS starting at min-degree nodes, neighbors sorted by degree
    visited = falses(nnodes)
    perm = Int[]

    while length(perm) < nnodes
        # start at unvisited node with minimal degree
        start = argmin([(visited[i] ? typemax(Int) : deg[i]) for i in 1:nnodes])
        queue = Int[]
        push!(queue, start)
        visited[start] = true
        head_idx = 1

        while head_idx <= length(queue)
            u = queue[head_idx]
            head_idx += 1
            push!(perm, u)
            # enqueue neighbors by increasing degree
            
            buffer = filter(v -> !visited[v], adjacents[u])
            sorted = sort!(buffer, by = v -> deg[v])
            for v in sorted
                visited[v] = true
                push!(queue, v)
            end
        end
    end

    reverse && reverse!(perm)
    mesh.nodes = mesh.nodes[perm]

    # Update nodal fields
    for (key, data) in mesh.node_data
        key == "node-id" && continue
        sz = size(data)
        if length(sz)==1
            mesh.node_data[key] = data[perm]
        else
            mesh.node_data[key] = data[perm,:]
        end
    end

    # Renumbering node ids
    for (i, node) in enumerate(mesh.nodes)
        node.id = i
    end
    mesh.node_data["node-id"]   = collect(1:length(mesh.nodes))

    # Sorting elements
    keys = [ minimum(getfield.(cell.nodes, :id)) for cell in mesh.elems]
    perm = sortperm(keys)
    mesh.elems = mesh.elems[perm]

    # Update cell fields
    for (key, data) in mesh.elem_data
        key == "elem-id" && continue
        sz = size(data)
        if length(sz)==1
            mesh.elem_data[key] = data[perm]
        else
            mesh.elem_data[key] = data[perm,:]
        end
    end

    # Renumbering element ids
    for (i, elem) in enumerate(mesh.elems)
        elem.id = i
    end
    mesh.elem_data["elem-id"]   = collect(1:length(mesh.elems))

    return mesh
end


function compute_facets(mesh::Mesh)

    if mesh.ctx.ndim==2
        mesh.edges = get_outer_facets(mesh.elems)
        mesh.faces = mesh.edges
    elseif mesh.ctx.ndim==3
        solids = filter( elem->elem.shape.ndim==3, mesh.elems )
        surfaces = filter( elem-> elem.role==:surface || (elem.shape.ndim==2 && elem.role==:bulk), mesh.elems )
        for cell in surfaces
            cell.owner = cell # set itself as owner
        end
        mesh.faces = [ get_outer_facets(solids); surfaces ]
        mesh.edges = get_edges(mesh.faces)
    end

    # add line elements as edges
    lines = filter( elem->elem.role==:line, mesh.elems )
    for line in lines
        line.owner = line # set itself as owner
    end
    append!(mesh.edges, lines)
end


# Syncs the mesh data
function synchronize!(mesh::Mesh; sort=false, cleandata=false)

    ndim = mesh.ctx.ndim
    if ndim!=3
        ndim = any( node.coord[3] != 0.0 for node in mesh.nodes ) ? 3 : 2
        if ndim == 2
            ndim = any( node.coord[2] != 0.0 for node in mesh.nodes ) ? 2 : 1
        end
    end

    mesh.ctx.ndim = max(ndim, mesh.ctx.ndim)

    # check if there is a 3d surface
    if mesh.ctx.ndim==2 && any( elem.role==:surface for elem in mesh.elems )
        notify("synchronize!: 2D mesh with surface elements detected. Converting to 3D mesh.")
        mesh.ctx.ndim = 3
    end

    # Numberig nodes
    for (i,p) in enumerate(mesh.nodes)
        p.id = i
    end

    # Numberig cells and setting ctx
    for (i,elem) in enumerate(mesh.elems)
        elem.id = i
        elem.ctx = mesh.ctx
    end

    # Faces and edges
    compute_facets(mesh)

    # Update node adjacents
    for node in mesh.nodes
        resize!(node.elems, 0)
    end
    for elem in mesh.elems
        for node in elem.nodes
            push!(node.elems, elem)
        end
    end
    for node in mesh.nodes
        unique!(node.elems)
    end

    # Quality
    Q = Float64[]
    for c in mesh.elems
        c.quality = cell_quality(c)
        push!(Q, c.quality)
    end
    
    # update data
    if cleandata
        empty!(mesh.node_data)
        empty!(mesh.elem_data)
    end
    
    mesh.node_data["node-id"]   = collect(1:length(mesh.nodes))
    mesh.elem_data["quality"]   = Q
    mesh.elem_data["elem-id"]   = collect(1:length(mesh.elems))
    mesh.elem_data["cell-type"] = [ cell.role in (:interface, :line_interface) ? Int32(VTK_POLY_VERTEX) : Int32(cell.shape.vtk_type) for cell in mesh.elems ]

    # Ordering
    sort && sort_mesh(mesh)
end



# Mesh quality
function quality!(mesh::Mesh)
    for c in mesh.elems
        c.quality = cell_quality(c)
    end
    Q = Float64[ c.quality for c in mesh.elems]
    mesh.elem_data["quality"] = Q
    return nothing
end


function join_mesh!(mesh::Mesh, m2::Mesh)

    # mesh.ctx.ndim = max(mesh.ctx.ndim, m2.ctx.ndim)

    pointdict = Dict{UInt, Node}()
    for m in (mesh, m2)
        for node in m.nodes
            pointdict[hash(node)] = node
        end
    end
    nodes = collect(values(pointdict))

    elems = Cell[]
    for m in (mesh, m2)
        for elem in m.elems
            newelemnodes = Node[pointdict[hash(node)] for node in elem.nodes ]
            newelem = Cell(elem.shape, newelemnodes, tag=elem.tag)
            push!(elems, newelem)
        end
    end

    mesh.nodes = nodes
    mesh.elems = elems
    mesh._pointdict = pointdict

    synchronize!(mesh, sort=false)

    return nothing
end


function join_meshes(m1::Mesh, m2::Mesh)
    mesh = copy(m1)
    join_mesh!(mesh, m2)
    return mesh
end


"""
    Mesh(coordinates, connectivities, cellshapes=CellShape[]; tag="", quiet=true)

Creates a `Mesh` from a nodal `coordinates` matrix,
an array of `connectivities` and a list of `cellshapes`.
If `cellshapes` are not provided, they are guessed based on the geometry.
A `tag` string for all generated cells can be provided optionally.

# Examples

```jldoctest
julia> using Serendip;
julia> C = [ 0 0; 1 0; 1 1; 0 1 ];
julia> L = [ [ 1, 2, 3 ], [ 1, 3, 4 ] ];
julia> S = [ TRI3, TRI3 ];
julia> Mesh(C, L, S, tag="triangle")

Mesh
  ndim: 2
  nodes: 4-element Vector{Node}:
    1: Node  id=1
    2: Node  id=2
    3: Node  id=3
    4: Node  id=4
  elems: 2-element Vector{Cell}:
    1: Cell  id=1  tag="triangle"
    2: Cell  id=2  tag="triangle"
  faces: 4-element Vector{Cell}:
    1: Cell  id=-1  tag="triangle"
    2: Cell  id=-1  tag="triangle"
    3: Cell  id=-1  tag="triangle"
    4: Cell  id=-1  tag="triangle"
  edges: 4-element Vector{Cell}:
    1: Cell  id=-1  tag="triangle"
    2: Cell  id=-1  tag="triangle"
    3: Cell  id=-1  tag="triangle"
    4: Cell  id=-1  tag="triangle"
  node_data: OrderedDict{String, Array} with 1 entry
    "node-id" => [1, 2, 3, 4]
  elem_data: OrderedDict{String, Array} with 3 entries
    "quality" => [0.8915188114208271, 0.8915188114208271]
    "elem-id" => [1, 2]
    "cell-type" => [5, 5]
```
"""


"""
    Mesh(coordinates, connectivities, cellshapes=CellShape[]; tag="", quiet=true) -> Mesh

Construct a `Mesh` object directly from nodal coordinates and element connectivities.

# Arguments
- `coordinates::Matrix{<:Real}`: Node coordinate matrix of size `(nnodes, ndim)`.
- `connectivities::Vector{Vector{Int}}`: Connectivity list, where each entry is a vector
  of node indices (1-based) defining an element.
- `cellshapes::Vector{CellShape}=CellShape[]`: Optional vector of element shapes (LIN2, QUAD4, etc.).  
  If not provided, shapes are inferred automatically from the number of nodes and
  spatial dimension.
- `tag::String=""`: Optional tag applied to all created elements.
- `quiet::Bool=false`: If `true`, suppresses console output during construction.

# Returns
- `mesh::Mesh`: A finite element mesh with nodes and elements created from
  the input data.

# Example
```julia
using Serendip

# square domain with 4 points and 2 triangular elements
coordinates = [ 0.0 0.0;
                1.0 0.0;
                1.0 1.0;
                0.0 1.0 ]

connectivities = [ [1, 2, 3], [1, 3, 4] ]
cellshapes     = [ TRI3, TRI3 ]

mesh = Mesh(coordinates, connectivities, cellshapes; tag="tri")
println(mesh)
```
"""
function Mesh(
              coordinates::Matrix{<:Real},
              connectivities::Vector{Vector{Int64}},
              cellshapes ::Vector{CellShape}=CellShape[];
              tag        ::String="",
              quiet      ::Bool=false,
             )


    n = size(coordinates, 1) # number of nodes
    m = size(connectivities , 1) # number of cells

    nodes = Node[]
    for i in 1:n
        C = coordinates[i,:]
        push!(nodes, Node(C))
    end

    # Get ndim
    ndim = size(coordinates,2)
    ctx  = MeshContext(ndim)

    cells = Cell[]
    for i in 1:m
        pts = nodes[connectivities[i]]
        if length(cellshapes)>0
            shape = cellshapes[i]
        else
            shape = get_shape_from_ndim_npoints(length(pts), ndim)
        end
        role = shape in (LIN2, LIN3, LIN4) ? :line : :bulk
        cell = Cell(shape, role, pts, tag=tag, ctx=ctx)
        push!(cells, cell)
    end

    mesh = Mesh()
    mesh.ctx.ndim = ndim
    mesh.nodes = nodes
    mesh.elems = cells
    synchronize!(mesh, sort=false) # no node ordering

    return mesh
end



export stats
function stats(mesh::Mesh)
    printstyled("Mesh stats:\n", bold=true, color=:cyan)

    npoints = length(mesh.nodes)
    ncells  = length(mesh.elems)
    @printf "  %3dd mesh                             \n" mesh.ctx.ndim
    @printf "  %4d nodes\n" npoints
    @printf "  %4d cells\n" ncells

    L = Float64[]
    for elem in mesh.elems.solids
        l = cell_extent(elem)^(1/mesh.ctx.ndim)
        push!(L, l)
    end
    lavg = mean(L)
    lmdn = quantile(L, 0.5)
    minl = minimum(L)
    maxl = maximum(L)

    bin = (maxl-minl)/10
    hist  = fit(Histogram, L, minl:bin:maxl, closed=:right).weights
    # @show hist
    lmod = (findmax(hist)[2]-1)*bin + bin/2

    @printf "  lmin = %7.5f\n" minl
    @printf "  lmax = %7.5f\n" maxl
    @printf "  lavg = %7.5f\n" lavg
    @printf "  lmdn = %7.5f\n" lmdn
    @printf "  lmod = %7.5f\n" lmod
end


# Construct a mesh given a set of elements
function Mesh(elems::Vector{Cell})
    length(elems)==0 && return Mesh()

    nodes = [ node for elem in elems for node in elem.nodes ]

    # check if all nodes have id
    all( node.id > 0 for node in nodes ) || error("Mesh: all input nodes must have a valid id")

    # newnodes
    nodes_d  = Dict{Int,Node}( node.id => copy(node) for node in nodes )
    newnodes = collect(values(nodes_d))
    node_ids = collect(keys(nodes_d))

    # copy elements
    newelems = Cell[]
    for elem in elems
        elem_nodes = Node[ nodes_d[node.id] for node in elem.nodes ]
        newelem    = Cell(elem.shape, elem_nodes, tag=elem.tag)
        push!(newelems, newelem)
    end

    ndim = any( node.coord[3] != 0.0 for node in newnodes ) ? 3 : 2

    newmesh       = Mesh(ndim)
    newmesh.nodes = newnodes
    newmesh.elems = newelems

    synchronize!(newmesh, sort=true)
    return newmesh
end


function threshold(mesh::Mesh, field::Union{Symbol,String}, minval::Float64, maxval::Float64)
    field = string(field)

    # check if field exists
    found = haskey(mesh.elem_data, field)
    if found
        vals = mesh.elem_data[field]
    else
        found = haskey(mesh.node_data, field)
        found || error("threshold: field $field not found")
        data  = mesh.node_data[field]
        vals = [ mean(data[i]) for i in 1:ncells ]
    end

    # filter cells
    cells = Cell[]
    for (cell, val) in zip(mesh.elems, vals)
        if minval <= val <= maxval
            push!(cells, cell)
        end
    end

    # get nodes
    nodes = get_nodes(cells)

    # ids from selected cells and nodes
    cids = [ c.id for c in cells ]
    pids = [ p.id for p in nodes ]

    # new mesh object
    new_mesh = Mesh()
    new_mesh.nodes = nodes
    new_mesh.elems = cells

    # select relevant data
    for (key,vals) in mesh.node_data
        new_mesh.node_data[key] = vals[pids,:]
    end

    for (key,vals) in mesh.elem_data
        new_mesh.elem_data[key] = vals[cids]
    end

    # update node numbering, facets and edges
    synchronize!(new_mesh, sort=false)

    return new_mesh

end


export get_segment_data
function get_segment_data(msh::AbstractDomain, X1::Array{<:Real,1}, X2::Array{<:Real,1}, filename::String=""; n=200)
    data = msh.node_data
    table = DataTable(["s"; collect(keys(data))])
    X1 = [X1; 0.0][1:3]
    X2 = [X2; 0.0][1:3]
    Δ = (X2-X1)/(n-1)
    Δs = norm(Δ)
    s1 = 0.0

    for i in 1:n
        X = X1 + Δ*(i-1)
        s = s1 + Δs*(i-1)
        cell = find_elem(X, msh.elems, msh._elempartition, 1e-7, Cell[])
        coords =get_coords(cell)
        R = inverse_map(cell.shape, coords, X)
        N = cell.shape.func(R)
        map = [ p.id for p in cell.nodes ]
        vals = [ s ]
        for (k,V) in data
            val = dot(V[map], N)
            push!(vals, val)
        end
        push!(table, vals)
    end

    filename != "" && save(table, filename)

    return table
end

export randmesh
function randmesh(n::Int...; shape=nothing)
    ndim = length(n)
    geo = GeoModel()
    if ndim==2
        lx, ly = (1.0, 1.0)
        nx, ny = n
        isnothing(shape) && (shape=rand((TRI3, TRI6, QUAD4, QUAD8)))
        add_block(geo, [0.0, 0.0], [lx, ly], nx=nx, ny=ny, shape=shape)
        # return Mesh(geo; quiet=true)
        # m = Mesh(Block([0.0 0.0; lx ly], nx=nx, ny=ny, shape=shape), quiet=true)
    else
        lx, ly, lz = (1.0, 1.0, 1.0)
        nx, ny, nz = n
        isnothing(shape) && (shape=rand((TET4, TET10, HEX8, HEX20)))
        add_block(geo, [0.0, 0.0, 0.0], [lx, ly, lz], nx=nx, ny=ny, nz=nz, shape=shape)
        # m = Mesh(Block([0.0 0.0 0.0; lx ly lz], nx=nx, ny=ny, nz=nz, shape=shape), quiet=true)
    end
    return Mesh(geo; quiet=true)
end





