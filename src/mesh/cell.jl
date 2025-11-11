# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl


# Cell
mutable struct Cell<:AbstractCell
    id     ::Integer
    shape  ::CellShape
    role::Symbol  # :vertex, :line, :bulk, :surface, :contact, :cohesive, :line_interface, :tip
    nodes  ::Vector{Node}
    tag    ::String
    active ::Bool
    quality::Float64              # quality index: surf/(reg_surf)
    embedded::Bool                # flag for embedded cells
    crossed::Bool                 # flag if cell crossed by linear inclusion
    owner  ::Union{AbstractCell,Nothing}  # owner cell if this cell is a face/edge
    couplings::Vector{AbstractCell}   # neighbor cells in case of interface cells
    ctx::MeshContext                 # mesh environment variables


    function Cell(shape::CellShape, role::Symbol, nodes::Vector{Node}; ctx::MeshContext=MeshContext(0), tag::String="", owner=nothing, id::Int=-1, active=true)
        this = new()
        this.id = id
        this.shape = shape
        this.role = role
        this.nodes = copy(nodes)
        this.tag = tag
        this.active  = active
        this.quality = 0.0
        this.embedded= false
        this.crossed = false
        this.owner   = owner
        this.couplings = []
        this.ctx     = ctx
        return this
    end
end

const CellFace=Cell
const CellEdge=Cell
const Facet=Cell


### Cell methods

Base.hash(cell::Cell) = sum(hash(node) for node in cell.nodes)
Base.isequal(c1::Cell, c2::Cell) = hash(c1)==hash(c2)


Base.copy(cell::Cell)  = Cell(cell.shape, cell.nodes, tag=cell.tag, owner=cell.owner)


function get_coords(c::AbstractCell, ndim=3)
    n = length(c.nodes)
    C = Array{Float64}(undef, n, ndim)
    for (i,p) in enumerate(c.nodes)
        C[i,1] = p.coord.x
        if ndim>1 C[i,2] = p.coord.y end
        if ndim>2 C[i,3] = p.coord.z end
    end
    return C
end


# Return all nodes in cells
# function get_nodes(cells::Array{<:AbstractCell,1})
#     return collect(Set(node for cell in cells for node in cell.nodes))
# end



function get_nodes(elems::Vector{<:AbstractCell})
    # get all nodes using object id
    points_d = Dict{UInt64, Node}()
    for elem in elems
        for node in elem.nodes
            points_d[objectid(node)] = node
        end
    end

    return collect(values(points_d))
end


"""
    select(elems::Vector{<:AbstractCell}, selectors...; invert=false, tag="")

Filters a list of finite element cells (`elems`) based on one or more `selectors`.

Selectors can be:
- `:all`             → select all elements (no filtering)
- `:bulk`, `:line`, `:contact`, `:cohesive`, `:line_interface`, `:tip` → select by element role
- `:active`          → select only active elements
- `:embedded`        → select embedded line elements (with couplings)
- `String`           → match element tag
- `Expr` or `Symbolic` → spatial condition using coordinates `x`, `y`, `z`
- `Vector{Int}`      → list of element indices to select
- `NTuple{N, Symbolic}` → multiple symbolic coordinate conditions

# Keyword Arguments
- `invert::Bool`: If `true`, returns the complement of the selected set.
- `tag::String`: If non-empty, assigns this tag to all selected elements.

# Returns
- A filtered list of elements matching the criteria.
"""
function select(
    elems::Vector{<:AbstractCell},
    selectors::Union{Symbol,Expr,Symbolic,String,Vector{Int},NTuple{N, Symbolic} where N}...;
    invert = false,
    tag::String = ""
    )

    selectors = flatten(selectors)
    selected = collect(1:length(elems)) # selected indexes

    for (i,selector) in enumerate(selectors)

        if isa(selector, Symbol)
            if selector == :all
                # do nothing (don't filter)
            elseif selector in (:bulk, :line, :cohesive, :contact, :line_interface, :tip)
                selected = Int[ i for i in selected if elems[i].role==selector ]
            elseif selector == :active
                selected = Int[ i for i in selected if elems[i].active ]
            elseif selector == :embedded
                selected = Int[ i for i in selected if elems[i].role==:line && length(elems[i].couplings)>0 ]
            else
                error("select: cannot filter array of Cell with symbol $(repr(selector))")
            end
        elseif isa(selector, String)
            selected = Int[ i for i in selected if elems[i].tag==selector ]
        elseif isa(selector, Expr) || isa(selector, Symbolic)
            nodes = [ node for cell in elems for node in cell.nodes ]
            max_id = maximum( n->n.id, nodes )
            pointmap = zeros(Int, max_id) # points and pointmap may have different sizes

            T = Bool[]
            for (i,node) in enumerate(nodes)
                pointmap[node.id] = i
                x, y, z = node.coord.x, node.coord.y, node.coord.z
                push!(T, evaluate(selector, x=x, y=y, z=z))
            end

            selected = Int[ i for i in selected if all( T[pointmap[node.id]] for node in elems[i].nodes ) ]
        elseif isa(selector, Vector{Int}) # selector is a vector of indexes
            selected = intersect(selected, selector)
        end
    end

    if invert
        selected = setdiff(1:length(elems), selected)
    end

    selection = elems[selected]

    # Set tag for selected elements
    if tag != ""
        for elem in selection
            elem.tag = tag
        end
    end

    return elems[selected]
end

"""
    select(domain::AbstractDomain, kind::Symbol, selectors...; invert=false, tag="")

Filters entities from a finite element domain (`domain`) by type and selection criteria.

# Arguments
- `domain::AbstractDomain`: The mesh or domain containing entities.
- `kind::Symbol`: One of `:element`, `:face`, `:edge`, or `:node` to specify which entity to select.
- `selectors`: One or more selectors.

Selectors can be:
- `:all`             → select all elements (no filtering)
- `:bulk`, `:line`, `:contact`, `:cohesive`, `:line_interface`, `:tip` → select by element role
- `:active`          → select only active elements
- `:embedded`        → select embedded line elements (with couplings)
- `String`           → match element tag
- `Expr` or `Symbolic` → spatial condition using coordinates `x`, `y`, `z`
- `Vector{Int}`      → list of element indices to select
- `NTuple{N, Symbolic}` → multiple symbolic coordinate conditions

# Keyword Arguments
- `invert::Bool`: If `true`, returns the entities not matching the selectors.
- `tag::String`: If non-empty, assigns this tag to all selected entities.

# Returns
- A list of selected entities of the specified kind.
"""
function select(
    domain::AbstractDomain,
    kind::Symbol,
    selectors::Union{Symbol,Expr,Symbolic,String,CellShape,Vector{<:Real},NTuple{N, Symbolic} where N}...;
    invert = false,
    nearest = true,
    tag::String = ""
    )

    if kind == :element
        return select(domain.elems, selectors...; invert=invert, tag=tag)
    elseif kind == :face
        return select(domain.faces, selectors...; invert=invert, tag=tag)
    elseif kind == :edge
        return select(domain.edges, selectors...; invert=invert, tag=tag)
    elseif kind == :node
        return select(domain.nodes, selectors...; invert=invert, nearest=nearest, tag=tag)
    elseif kind == :ip
        ips = [ ip for elem in domain.elems for ip in elem.ips ]
        return select(ips, selectors...; invert=invert, nearest=nearest, tag=tag)
    else
        error("select: unknown kind $(repr(kind))")
    end

end


function cellnormal(cell::AbstractCell)
    iscoplanar(cell) || return nothing

    tol = 1e-8
    ndim = cell.shape.ndim+1
    C = get_coords(cell, ndim)
    I = ones(size(C,1))
    N = pinv(C.+tol/100)*I # best fit normal
    normalize!(N)
    return N
end


function isparallelto(A,B)
    tol = 1e-8

    dotAB = dot(A,B)
    normAB = norm(A)*norm(B)
    abs(dotAB-normAB) < tol && return true
    abs(dotAB+normAB) < tol && return true
    return false
end


function iscoplanar(cell::AbstractCell)
    tol = 1e-8

    coords = get_coords(cell, 3)

    # find a plane
    X1 = coords[1,:]
    X2 = coords[2,:]
    X1X2 = X2-X1

    # look for a non-collinear point
    local X, N
    for i in 3:length(cell.nodes)
        X = coords[i,:]
        X1X = X-X1
        N = cross(X1X2, X1X)
        norm(N) > tol && break
    end

    # test the plane at each point
    for i in 3:length(cell.nodes)
        X = coords[i,:]

        if dot(X-X1, N) > tol
            return false
        end
    end

    return true
end


function nearest(cells::Vector{Cell}, coord)
    n = length(cells)
    D = zeros(n)
    X = vec(coord)

    for (i,cell) in enumerate(cells)
        C = vec(mean(get_coords(cell), dims=1))
        D[i] = norm(X-C)
    end

    return cells[sortperm(D)[1]]
end



# Gets the coordinates of a bounding box for an array of nodes
function bounding_box(nodes::Vector{Node})
    minx = miny = minz =  Inf
    maxx = maxy = maxz = -Inf
    for node in nodes
        node.coord.x < minx && (minx = node.coord.x)
        node.coord.y < miny && (miny = node.coord.y)
        node.coord.z < minz && (minz = node.coord.z)
        node.coord.x > maxx && (maxx = node.coord.x)
        node.coord.y > maxy && (maxy = node.coord.y)
        node.coord.z > maxz && (maxz = node.coord.z)
    end
    return [ minx miny minz; maxx maxy maxz ]
end


# Gets the coordinates of a bounding box for a cell
function bounding_box(cell::AbstractCell)
    return bounding_box(cell.nodes)
end


# Gets the coordinates of a bounding box for an array of cells
function bounding_box(cells::Array{<:AbstractCell,1})
    nodes = unique( Node[ p for c in cells for p in c.nodes ] )
    return bounding_box(nodes)
end


# gets all facets of a cell
function get_facets(cell::AbstractCell)
    faces  = Cell[]
    all_facets_idxs = cell.shape.facet_idxs
    facet_shape    = cell.shape.facet_shape

    facet_shape==() && return faces

    sameshape = typeof(facet_shape) == CellShape # check if all facets have the same shape
    role = cell.shape.ndim ==3 ? :surface : :line # 2D cells have only lines as facets

    # Iteration for each facet
    for (i, face_idxs) in enumerate(all_facets_idxs)
        nodes = cell.nodes[face_idxs]
        shape  = sameshape ? facet_shape : facet_shape[i]
        face   = Cell(shape, role, nodes, tag=cell.tag, owner=cell)
        face.nodes = nodes # update nodes since Cell creates a copy
        push!(faces, face)
    end

    return faces
end


# gets all edges of a cell
function get_edges(cell::AbstractCell)
    edges  = Cell[]
    all_edge_idxs = cell.shape.edge_idxs

    for edge_idx in all_edge_idxs
        nodes = cell.nodes[edge_idx]
        shape  = (LIN2, LIN3, LIN4)[length(nodes)-1]
        edge   = Cell(shape, :line, nodes, tag=cell.tag, owner=cell)
        push!(edges, edge)
    end

    return edges
end


# Returns the volume/area/length of a cell
function cell_extent(c::AbstractCell)
    IP = get_ip_coords(c.shape)
    nip = size(IP,1)
    nldim = c.shape.ndim # cell basic dimension

    # get coordinates matrix
    C = get_coords(c)
    J = Array{Float64}(undef, size(C,2), nldim)

    # calc metric
    vol = 0.0
    for i in 1:nip
        R    = IP[i].coord
        w    = IP[i].w
        dNdR = c.shape.deriv(R)
        @mul J = C'*dNdR
        normJ = norm2(J)
        #if normJ<0
            #error("cell_extent: Negative Jacobian while calculating cell volume/area/length id=$(c.id) shape=$(c.shape.name) ")
        #end
        vol += normJ*w
    end
    return vol
end


# Returns the surface/perimeter of a regular element given the volume/area of a cell
function regular_surface(metric::Float64, shape::CellShape)
    if shape in [ TRI3, TRI6, TRI9, TRI10 ]
        A = metric
        a = 2.0*√(A/√3.0)
        return 3*a
    end
    if shape in [ QUAD4, QUAD8, QUAD9, QUAD12, QUAD16 ]
        A = metric
        a = √A
        return 4*a
    end
    if shape in [ PYR5, PYR13 ]
        V = metric
        a = ( 3.0*√2.0*V )^(1.0/3.0)
        return (1 + √3.0)*a^2
    end
    if shape in [ TET4, TET10 ]
        V = metric
        a = ( 6.0*√2.0*V )^(1.0/3.0)
        return √3.0*a^2
    end
    if shape in [ HEX8, HEX20, HEX27 ]
        V = metric
        a = V^(1.0/3.0)
        return 6.0*a^2.0
    end
    if shape in [ WED6, WED15 ]
        V = metric
        a2 = (16.0/3.0*V^2)^(1.0/3.0)
        return (3.0 + √3.0/2.0)*a2
    end
    error("No regular surface/perimeter value for shape $(shape.name)")
end


# Returns the area/volume of a regular element given the perimeter/surface
function regular_volume(metric::Float64, shape::CellShape)
    if shape in [ TRI3, TRI6, TRI9, TRI10 ]
        p = metric
        a = p/3
        return a^2/4*√3.0
    end
    if shape in [ QUAD4, QUAD8, QUAD9, QUAD12, QUAD16 ]
        p = metric
        a = p/4
        return a^2
    end
    #if shape in [ PYR5 ]
        #V = metric
        #a = ( 3.0*√2.0*V )^(1.0/3.0)
        #return (1 + √3.0)*a^2
    #end
    if shape in [ TET4, TET10 ]
        s = metric
        A = s/4
        a = 2.0*√(A/√3.0)
        return a^3/(6*√2.0)
    end
    if shape in [ HEX8, HEX20 ]
        s = metric
        A = s/6
        a = √A
        return a^3
    end
    #if shape in [ WED6, WED15 ]
        #s = metric
        #return (3.0 + √3.0/2.0)*a2
    #end
    error("No regular area/volume value for shape $(get_name(shape))")
end



# Returns the cell quality ratio as vol/reg_vol
function cell_quality_2(c::AbstractCell)::Float64
    # get faces
    faces = get_facets(c)
    length(faces)==0 && return 1.0

    # cell surface
    surf = sum( cell_extent(f) for f in faces )

    # quality calculation
    vol = cell_extent(c) # volume or area
    rvol = regular_volume(surf, c.shape)
    return min(vol/rvol, 1.0)
end


# Returns the cell quality ratio as reg_surf/surf
function cell_quality(c::AbstractCell)::Float64
    # get faces
    c.role != :bulk && return 1.0

    faces = get_facets(c)
    length(faces)==0 && return 1.0

    # cell surface
    surf = 0.0
    for f in faces
        surf += cell_extent(f)
    end

    # quality calculation
    extent = cell_extent(c) # volume or area
    rsurf  = regular_surface(abs(extent), c.shape)

    k = 2
    q = sign(extent)*(rsurf/surf)^k
    if q>1
        q = 2 - q
    end
    return q
end


function cell_aspect_ratio(c::AbstractCell)::Float64
    edges = get_edges(c)
    L = [ cell_extent(e) for e in edges ]
    return maximum(L)/minimum(L)
end


# Get an array with shares for all nodes
function get_patches(cells::Array{<:AbstractCell,1})
    # get all nodes from cells if needed
    pointsd = Dict{UInt64, Node}()
    for cell in cells
        for node in cell.nodes
            pointsd[hash(node)] = node
        end
    end

    nodes = collect(values(pointsd))
    np     = length(nodes)

    # backup nodes ids
    bk_pt_id = [ pt.id for pt in nodes ]
    for i in 1:np
        nodes[i].id = i
    end

    # get incidence array
    patches  = [ AbstractCell[] for i in 1:np ]
    for cell in cells
        for pt in cell.nodes
            push!(patches[pt.id], cell)
        end
    end

    # restore nodes ids
    for i in 1:np
        nodes[i].id = bk_pt_id[i]
    end

    return nodes, patches
end


function inverse_map(cell::AbstractCell, X::AbstractVector{Float64}, tol=1.0e-7)
    return inverse_map(cell.shape, get_coords(cell), X, tol)
end


# function get_point(s::Float64, coords::Array{Float64,2})
#     #  Interpolates coordinates for s between 0 and 1
#     #
#     #  0               +1  -->s
#     #  1---2---3---..---n  -->idx

#     @assert 0<=s<=1
#     n = size(coords,1)
#     m = n - 1 # number of segments
#     δ = 1/m   # parametric length per segment

#     # segment index and parametric coordinates
#     i  = floor(Int, s/δ) + 1 # index for current segment
#     i  = min(i, m)
#     si = (i-1)*δ     # global parametric coordinate at the beginning of segment i
#     t  = (s - si)/δ  # local parametric coordiante

#     return coords[i,:]*(1-t) + coords[i+1,:]*t
# end
