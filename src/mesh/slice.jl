"""
    getpolygon(nodes)

Orders an array of `nodes` to get a closed polygon.
"""
function getpolygon(nodes::Vector{Node})
    # the nodes might be placed an arbitrary 3D plane

    # find the best fit plane
    C = get_coords(nodes, 3)
    nnodes = length(nodes)

    # move the coordinates to avoid singular case
    # when the regression line/planes crosses the origin
    Cm =  C .+ [pi 2*pi 3*pi]

    I = ones(nnodes)
    N = pinv(Cm)*I # best fit normal
    N = round.(N, digits=15) # to avoid almost zero values

    # rotation matrix
    if N[1]==0.0
        V1 = [ 1.0, 0.0, 0.0 ]
    elseif N[2]==0.0
        V1 = [ 0.0, 1.0, 0.0 ]
    else
        V1 = [ 0.0, 0.0, 1.0 ]
    end

    V2 = cross(N, V1)
    V3 = cross(V1, V2)

    normalize!(V2)
    normalize!(V3)

    R = [ V1 V2 V3 ]
    C = C*R

    # get center
    Xc = vec(sum(C, dims=1))./nnodes

    # get polar angles
    angles = zeros(nnodes)
    for i in 1:nnodes
        x, y, _ = C[i,:]-Xc
        angles[i] = atan(y,x)
    end

    perm = sortperm(angles)
    C = C[perm,:]

    # get internal angles
    minα = 2*pi
    mini = 0
    for i in 1:nnodes
        prev = i==1 ? nnodes : i-1
        next = i==nnodes ? 1 : i+1
        X = C[i,:]
        V1 = normalize(C[prev,:]-X)
        V2 = normalize(C[next,:]-X)
        α = acos(clamp(dot(V1,V2), -1, 1))
        if α < minα
            minα=α
            mini = i
        end
    end

    # put point with smaller internal angle as first
    circshift(perm, 1-mini)

    # Serendip.@showm get_coords(nodes[perm])*R

    # build contour
    return nodes[perm]
end



function sortnodes(nodes::Vector{Node}, Vx´, Vy´)
    # find the best fit plane
    C = get_coords(nodes, 3)
    nnodes = length(nodes)

    # get center
    Xc = vec(sum(C, dims=1))./nnodes
    angles = zeros(nnodes) # polar angles

    # get polar angles and update C
    for i in 1:nnodes
        X = C[i,:]
        x´ = dot(X - Xc, Vx´)
        y´ = dot(X - Xc, Vy´)
        C[i,:] = [x´, y´, 0.0]
        angles[i] = atan(y´, x´)
    end

    perm = sortperm(angles)
    C = C[perm,:]

    angles = zeros(nnodes) # polygon internal angles
    for i in 1:nnodes
        prev = mod1(i-1, nnodes)
        next = mod1(i+1, nnodes)
        X    = C[i,:]
        V1   = normalize(C[prev,:]-X)
        V2   = normalize(C[next,:]-X)
        α    = acos(clamp(dot(V1,V2), -1, 1))
        angles[i] = α
    end
    _, min_idx = findmin(angles)

    # put point with smaller internal angle as first
    circshift(perm, 1-min_idx)

    # build contour
    return nodes[perm]
end


function makecells(nodes::Vector{Node}; quadratic=false)
    cells = Cell[]
    nnodes = length(nodes)
    # @assert nnodes==shape.npoints

        # nodes = getpolygon(nodes)
    if nnodes==3
        cell = Cell(TRI3, :solid, nodes)
        push!(cells, cell)
    elseif nnodes==4
        cell = Cell(QUAD4, :solid, nodes)
        push!(cells, cell)
    elseif nnodes==6
        cell1 = Cell(QUAD4, :solid, nodes[1:4])
        cell2 = Cell(QUAD4, :solid, nodes[[4,5,6,1]])
        push!(cells, cell1)
        push!(cells, cell2)
    elseif nnodes==5
        cell1 = Cell(TRI3, :solid, nodes[1:3])
        cell2 = Cell(QUAD4, :solid, nodes[[3,4,5,1]])
        push!(cells, cell1)
        push!(cells, cell2)
    end

    return cells

end


"""
    slice(mesh; base, axis)

Generates a planar mesh by slicing a 3D `mesh` using a plane defined by a `base` point and an `axis`.
The original nodal data is interpolated to the nodes of the resulting mesh.
"""
function slice(
    mesh::AbstractDomain;
    base::AbstractArray{<:Real,1} = Float64[],
    axis::AbstractArray{<:Real,1} = Float64[],
    project=true # makes a x-y projection of the slice
    )

    ndim = mesh.ctx.ndim
    @check ndim==3
    @check length(base)==3
    @check length(axis)==3

    axis = normalize(Vec3(axis))
    base = Vec3(base)

    # find projection system
    if axis[3]==1.0
        Vx´ = Vec3(1,0,0)
        Vy´ = Vec3(0,1,0)
    elseif axis[3]==-1
        Vx´ = -Vec3(1,0,0)
        Vy´ = Vec3(0,1,0)
    elseif axis[3]==0.0
        Vx´ = cross(Vec3(0,0,1), axis)
        Vy´ = Vec3(0,0,1)
    else
        pr = Vec3(axis[1], axis[2], 0)
        Vx´ = normalize(cross(axis, pr))
        Vy´ = normalize(cross(axis, Vx´))
    end

    # get all mesh edges
    edgedict = Dict{Tuple{Vararg{UInt64}}, CellEdge}()
    for cell in mesh.elems
        cell.role == :solid || continue
        for edge in get_edges(cell)
            edgedict[_topology_key(edge)] = edge
        end
    end

    edges = collect(values(edgedict))

    # find intersection points between edges and the plane
    int_nodes_d = NodePosMap()
    edge_node_d = Dict{Tuple{Vararg{UInt64}}, Node}()
    tolf = 1e-8
    maxits = 100
    for edge in edges

        Xa = edge.nodes[1].coord
        Xb = edge.nodes[2].coord

        fa = dot(Xa-base, axis)
        fb = dot(Xb-base, axis)
        fa*fb > 0.0 && continue

        Xi = nothing

        if fa==0.0
           Xi = Xa
        end
        if fb==0.0
            Xi = Xb
        end

        if Xi===nothing
            Xi = Vec3(0,0,0)
            C  = get_coords(edge)
            a  = -1.0
            b  = +1.0
            fi = 0.0
            for i in 1:maxits
                ξ = (a + b)/2
                N = edge.shape.func([ξ])
                Xi = C'*N

                fi = dot(Xi - base, axis)
                if fa*fi < 0.0
                    fb = fi
                    b  = ξ
                else
                    fa = fi
                    a  = ξ
                end
                abs(fi) < tolf && break
            end
        end

        Xi = round.(Xi, digits=6)
        int_node = Node(Xi)
        int_node = get!(int_nodes_d, int_node, int_node)
        edge_node_d[_topology_key(edge)] = int_node
    end

    # mount slice elements
    bulk_elems = Cell[]
    for cell in mesh.elems
        cell.role == :solid || continue

        nodedict = IdDict{Node,Nothing}()
        for edge in get_edges(cell)
            Xi = get(edge_node_d, _topology_key(edge), nothing)
            Xi === nothing && continue
            nodedict[Xi] = nothing
        end

        nodes = collect(keys(nodedict))
        nnodes = length(nodes)
        nnodes > 2 || continue

        nodes = sortnodes(nodes, Vx´, Vy´)
        cells = makecells(nodes)
        for c in cells
            c.tag = cell.tag
            c.owner = cell
            push!(bulk_elems, c)
        end
    end

    bulk_elems = unique(_topology_key, bulk_elems)
    new_nodes = get_nodes(bulk_elems)

    # non bulk elements
    d = dot(base, axis)
    l, m, n = axis
    # other_elems = mesh.elems[:($l*x + $m*y + $n*z == $d) ]
    other_elems = select(mesh, :element, :( $l*x + $m*y + $n*z == $d ) )
    new_nodes = [ new_nodes; get_nodes(other_elems) ]

    # Numberig nodes
    for (i, node) in enumerate(new_nodes)
        node.id = i
    end

    # Numberig cells and setting ctx
    for (i, elem) in enumerate(bulk_elems)
        elem.id = i
    end

    # interpolate node data
    labels = [ key for key in keys(mesh.node_fields) if key != "node-id" ]
    nnodes = length(new_nodes)
    node_fields = OrderedDict{String,Array}()


    for key in labels
        data = mesh.node_fields[key]
        sz   = size(data)
        dim  = length(sz)
        if dim==1
            newdata = zeros(eltype(data), nnodes)
        else
            newdata = zeros(eltype(data), nnodes, sz[2])
        end
        node_fields[key] = newdata
    end

    count = zeros(Int, nnodes)

    for cell in bulk_elems
        ocell = cell.owner
        ids = [ node.id for node in ocell.nodes ]
        coords = get_coords(ocell)
        for node in cell.nodes
            count[node.id] > 0 && continue
            count[node.id] += 1
            Ξ = inverse_map(ocell.shape, coords, node.coord)
            N = ocell.shape.func(Ξ)
            for key in labels
                V = mesh.node_fields[key][ids,:]
                val = V'*N
                node_fields[key][node.id, :] .= val
            end
        end
    end


    for key in labels
        for cell in other_elems
            for node in cell.nodes
                node_fields[key][node.id, :] .= mesh.node_fields[key][node.id, :]
            end
        end
    end


    # update cell data
    elem_fields=OrderedDict{String,Array}()
    nelems = length(bulk_elems)
    for (key,data) in mesh.elem_fields
        key in ("elem-id", "quality", "cell-type") && continue
        contains(key, "tag") && continue

        sz = size(data)
        dim = length(sz)
        if dim==1
            newdata = zeros(eltype(data), nelems)
        else
            newdata = zeros(eltype(data), nelems, sz[2])
        end

        for cell in bulk_elems
            ocell = cell.owner
            val = data[ocell.id, :]
            newdata[cell.id, :] .= val
        end

        elem_fields[key] = newdata
    end

    if project
        for node in new_nodes
            X = node.coord
            x´ = dot(X-base, Vx´)
            y´ = dot(X-base, Vy´)
            node.coord = Vec3(x´, y´, 0.0)
        end

        # update U so the slice can be plotted using warp
        if haskey(node_fields, "U")
            U = node_fields["U"]
            R = [ Vx´ Vy´ zeros(3) ]'
            for i in 1:size(U, 1)
                U[i, :] = R*U[i, :]
            end
        end
    end

    ndim = project ? 2 : 3
    newmesh = Mesh(ndim)

    newmesh.nodes = new_nodes
    newmesh.elems = [ bulk_elems; other_elems ]
    newmesh.node_fields = node_fields
    newmesh.elem_fields = elem_fields
    synchronize(newmesh, sort=false)

    return newmesh
end
