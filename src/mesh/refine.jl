function _refine_node!(nodes::Vector{Node}, nodes_d::NodePosMap, coord, shared::Bool)
    node = Node(coord)

    if shared
        existing = get(nodes_d, node, nothing)
        if existing === nothing
            nodes_d[node] = node
            push!(nodes, node)
        else
            node = existing
        end
    else
        push!(nodes, node)
    end

    return node
end


function _refine_cell(cell::Cell, shape::CellShape, nodes::Vector{Node})
    newcell = Cell(shape, cell.role, nodes, tag=cell.tag, active=cell.active)
    newcell.embedded = cell.embedded
    newcell.crossed = cell.crossed
    return newcell
end


function refine(mesh::Mesh; mode::Symbol=:h, n=2, quiet=true)
    if mode == :h
        return hrefine(mesh, n=n, quiet=quiet)
    elseif mode == :p
        return prefine(mesh, n=n, quiet=quiet)
    end

    error("refine: unknown refinement mode $(repr(mode)). Expected :h or :p")
end


function hrefine(mesh::Mesh; n=2, quiet=true)
    n == 1 && return copy(mesh)
    n > 0 || error("hrefine: n must be positive")

    newmesh = Mesh(mesh.ctx.ndim)
    nodes = Node[]
    nodes_d = NodePosMap()
    cells = Cell[]

    for cell in mesh.elems
        coords = get_coords(cell.nodes)

        if cell.shape == TRI3
            p_arr = Array{Node}(undef, n+1, n+1)
            for j = 1:n+1
                for i = 1:n+1
                    i + j > n + 2 && continue

                    r = (i - 1) / n
                    s = (j - 1) / n
                    C = cell.shape.func([r, s])' * coords
                    shared = i == 1 || j == 1 || i + j == n + 2
                    p_arr[i, j] = _refine_node!(nodes, nodes_d, C, shared)
                end
            end

            for j = 1:n
                for i = 1:n
                    i + j >= n + 2 && continue

                    p1 = p_arr[i  , j  ]
                    p2 = p_arr[i+1, j  ]
                    p3 = p_arr[i  , j+1]
                    push!(cells, _refine_cell(cell, TRI3, [p1, p2, p3]))

                    if i + j < n + 1
                        p4 = p_arr[i+1, j+1]
                        push!(cells, _refine_cell(cell, TRI3, [p2, p4, p3]))
                    end
                end
            end
            continue
        end

        if cell.shape == QUAD4
            p_arr = Array{Node}(undef, n+1, n+1)
            for j = 1:n+1
                for i = 1:n+1
                    r = -1.0 + 2.0 * (i - 1) / n
                    s = -1.0 + 2.0 * (j - 1) / n
                    C = cell.shape.func([r, s])' * coords
                    shared = i in (1, n+1) || j in (1, n+1)
                    p_arr[i, j] = _refine_node!(nodes, nodes_d, C, shared)
                end
            end

            for j = 1:n
                for i = 1:n
                    p1 = p_arr[i  , j  ]
                    p2 = p_arr[i+1, j  ]
                    p3 = p_arr[i+1, j+1]
                    p4 = p_arr[i  , j+1]
                    push!(cells, _refine_cell(cell, QUAD4, [p1, p2, p3, p4]))
                end
            end
            continue
        end

        if cell.shape == HEX8
            p_arr = Array{Node}(undef, n+1, n+1, n+1)
            for k = 1:n+1
                for j = 1:n+1
                    for i = 1:n+1
                        r = -1.0 + 2.0 * (i - 1) / n
                        s = -1.0 + 2.0 * (j - 1) / n
                        t = -1.0 + 2.0 * (k - 1) / n
                        C = cell.shape.func([r, s, t])' * coords
                        shared = i in (1, n+1) || j in (1, n+1) || k in (1, n+1)
                        p_arr[i, j, k] = _refine_node!(nodes, nodes_d, C, shared)
                    end
                end
            end

            for k = 1:n
                for j = 1:n
                    for i = 1:n
                        p1 = p_arr[i  , j  , k  ]
                        p2 = p_arr[i+1, j  , k  ]
                        p3 = p_arr[i+1, j+1, k  ]
                        p4 = p_arr[i  , j+1, k  ]
                        p5 = p_arr[i  , j  , k+1]
                        p6 = p_arr[i+1, j  , k+1]
                        p7 = p_arr[i+1, j+1, k+1]
                        p8 = p_arr[i  , j+1, k+1]
                        push!(cells, _refine_cell(cell, HEX8, [p1, p2, p3, p4, p5, p6, p7, p8]))
                    end
                end
            end
            continue
        end

        error("hrefine: Cannot refine mesh containing elements of type $(cell.shape.kind)")
    end

    newmesh.nodes = nodes
    newmesh.elems = cells
    synchronize(newmesh, sort=true)

    return newmesh
end


function prefine(mesh::Mesh; n=2, quiet=true)
    _ = (n, quiet)

    newshape_d = Dict{CellShape, CellShape}(
        TRI3  => TRI6,
        QUAD4 => QUAD8,
        TET4  => TET10,
        HEX8  => HEX20,
    )

    cells = Cell[]
    nodes = Node[]
    nodes_d = NodePosMap()

    for cell in mesh.elems
        newshape = get(newshape_d, cell.shape, nothing)
        newshape === nothing && error("prefine: Cannot refine mesh containing elements of type $(cell.shape.kind)")

        coords = get_coords(cell.nodes)
        points = Node[]
        for i in 1:newshape.npoints
            R = newshape.nat_coords[i, :]
            C = coords' * cell.shape.func(R)
            node = Node(C)
            node = get!(nodes_d, node, node)
            push!(points, node)
        end

        push!(cells, _refine_cell(cell, newshape, points))
    end

    nodes = collect(values(nodes_d.store))

    newmesh = Mesh(mesh.ctx.ndim)
    newmesh.nodes = nodes
    newmesh.elems = cells
    synchronize(newmesh, sort=true)

    return newmesh
end
