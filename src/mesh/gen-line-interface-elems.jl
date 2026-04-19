# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

const PATH_NODE_TOL = 1e-8
const HOST_BISECT_TOL = 1e-9
const HOST_PATH_TOL = 1e-4
const HOST_CELL_TOL = 1e-9
const HOST_PROBE_STEP = 1.0

"""
    gen_path_insets(mesh::Mesh, subpath::GPath)
    gen_path_insets(mesh::Mesh, subpaths::Vector{GPath})

Discretize one or more geometric paths into line elements inside `mesh`.

The generation strategy is controlled by `subpath.mode`:
- `:interface` traces the path through host cells and creates line-interface
  cells plus optional tip cells.
- `:embedded` traces the path through host cells and embeds the generated line
  cells in their hosts.
- `:conforming` walks existing mesh topology and creates standalone line cells
  that reuse matching mesh edges.
- `:free` creates standalone line cells from path commands, reusing existing
  endpoint nodes when present and creating missing nodes otherwise.
"""
function gen_path_insets(mesh::Mesh, subpaths::Vector{GPath})
    for subpath in subpaths
        gen_path_insets(mesh, subpath)
    end
end


function gen_path_insets(mesh::Mesh, subpath::GPath)
    if subpath.mode in (:interface, :embedded)
        gen_hosted_path_insets(mesh, subpath)
    elseif subpath.mode == :conforming
        gen_conforming_path_insets(mesh, subpath)
    elseif subpath.mode == :free
        gen_free_path_insets(mesh, subpath)
    else
        error("gen_path_insets: invalid path mode $(repr(subpath.mode))")
    end
end


function _find_mesh_node(mesh::Mesh, X; tol=PATH_NODE_TOL, required=false, label="point")
    matches = Node[]
    for node in mesh.nodes
        norm(node.coord - X) <= tol && push!(matches, node)
    end

    if isempty(matches)
        required && error("gen_path_insets: no existing mesh node found for $label at $(X)")
        return nothing
    end
    length(matches) > 1 && error("gen_path_insets: multiple mesh nodes found for $label at $(X)")
    return matches[1]
end


function _find_or_create_mesh_node!(mesh::Mesh, X; tol=PATH_NODE_TOL)
    node = _find_mesh_node(mesh, X, tol=tol)
    node !== nothing && return node

    node = Node(X)
    push!(mesh.nodes, node)
    return node
end


function _path_command_endpoints(path::Path, cmd::PathCmd)
    return evaluate(path, cmd, 0.0), evaluate(path, cmd, 1.0)
end


function _oriented_line_nodes(edge::AbstractCell, forward::Bool)
    forward && return copy(edge.nodes)

    n = length(edge.nodes)
    if n == 2
        return edge.nodes[[2, 1]]
    elseif n == 3
        return edge.nodes[[2, 1, 3]]
    elseif n == 4
        return edge.nodes[[2, 1, 4, 3]]
    else
        error("gen_path_insets: cannot orient edge with $n nodes")
    end
end


function _matching_oriented_edge(path::Path, cmd::PathCmd, edge::AbstractCell, current_node::Node, target_t::Float64; tol=PATH_NODE_TOL)
    n1, n2 = edge.nodes[1], edge.nodes[2]

    if current_node === n1
        forward = true
        next_node = n2
    elseif current_node === n2
        forward = false
        next_node = n1
    else
        return nothing
    end

    t_current = command_parameter(path, cmd, current_node.coord, tol=tol)
    t_next = command_parameter(path, cmd, next_node.coord, tol=tol)
    t_current === nothing && return nothing
    t_next === nothing && return nothing
    t_next > t_current + tol || return nothing
    t_next <= target_t + tol || return nothing

    for node in edge.nodes
        t = command_parameter(path, cmd, node.coord, tol=tol)
        t === nothing && return nothing
        t_current - tol <= t <= t_next + tol || return nothing
    end

    return (_oriented_line_nodes(edge, forward), next_node)
end


_line_shape_for_nodes(nodes::Vector{Node}) = (LIN2, LIN3, LIN4)[length(nodes)-1]


function _path_line_nodes!(mesh::Mesh, path::Path, cmd::PathCmd, line_shape::CellShape, P1::Node, P2::Node; t1=0.0, t2=1.0)
    line_shape in (LIN2, LIN3, LIN4) || error("gen_path_insets: unsupported line shape $(line_shape.kind)")
    line_shape == LIN2 && return [P1, P2]

    points = [P1, P2]
    for i in 3:line_shape.npoints
        r = line_shape.nat_coords[i, 1]
        t = t1 + 0.5 * (r + 1.0) * (t2 - t1)
        P = Node(evaluate(path, cmd, t))
        push!(mesh.nodes, P)
        push!(points, P)
    end
    return points
end


function _add_line_cell!(mesh::Mesh, shape::CellShape, nodes::Vector{Node}, tag::String)
    line_cell = Cell(shape, :line, nodes, tag=tag)
    push!(mesh.elems, line_cell)
    return line_cell
end


function _add_line_interface_cell!(mesh::Mesh, host_cell::AbstractCell, line_cell::AbstractCell, tag::String)
    joint_cell = Cell(POLYVERTEX, :line_interface, vcat(host_cell.nodes, line_cell.nodes), tag=tag)
    joint_cell.couplings = [host_cell, line_cell]
    push!(mesh.elems, joint_cell)
    return joint_cell
end


function _add_tip_cell!(mesh::Mesh, host_cell::AbstractCell, line_cell::AbstractCell, tip::Node, tag::String)
    tip_cell = Cell(POLYVERTEX, :tip, vcat(host_cell.nodes, tip), tag=tag)
    tip_cell.couplings = [host_cell, line_cell]
    push!(mesh.elems, tip_cell)
    return tip_cell
end


function _attach_hosted_line!(
    mesh::Mesh,
    host_cell::AbstractCell,
    line_cell::AbstractCell,
    subpath::GPath,
    P1::Node,
    P2::Node;
    is_path_start::Bool,
    is_path_end::Bool,
    is_closed::Bool,
)
    if subpath.mode == :embedded
        line_cell.embedded = true
        line_cell.couplings = [host_cell]
    else
        _add_line_interface_cell!(mesh, host_cell, line_cell, subpath.interface_tag)

        if !is_closed && is_path_start && subpath.tips in (:start, :both)
            _add_tip_cell!(mesh, host_cell, line_cell, P1, subpath.tip_tag)
        end
        if !is_closed && is_path_end && subpath.tips in (:end, :both)
            _add_tip_cell!(mesh, host_cell, line_cell, P2, subpath.tip_tag)
        end
    end

    host_cell.crossed = true
    return line_cell
end


function gen_hosted_path_insets(mesh::Mesh, subpath::GPath)
    path      = subpath.path
    commands  = path.cmds
    is_closed = path.closed
    line_tag  = subpath.tag
    line_shape    = subpath.shape

    segment_end = 1.0
    first_path_node = Node(path.points[1].coord)
    last_path_node = nothing
    is_path_start = true
    is_path_end = false
    for (cmd_idx, cmd) in enumerate(commands)
        cmd.key == :M && continue

        # find the initial element
        X_start = evaluate(path, cmd, HOST_PATH_TOL) # a little bit ahead from 0.0
        current_cell = find_elem(X_start, mesh.elems, mesh._elempartition, tol=HOST_CELL_TOL) # first crossed cell

        is_last_segment = false
        is_last_command = cmd_idx == length(commands)
        ξ = 0.0
        ξ_prev = 0.0

        # Splitting cmd
        while true
            if current_cell !== nothing
                current_cell_coords = get_coords(current_cell)
                # Default step
                step = 0.5*(1.0 - ξ)

                # Finding step (introduced to catch elements in curved paths)
                ξ_trial = ξ
                n_probe_its = round(Int, 1.0 / HOST_PROBE_STEP)
                for _ in 1:n_probe_its
                    ξ_trial += HOST_PROBE_STEP
                    ξ_trial > 1.0 && break
                    X = evaluate(path, cmd, ξ_trial)
                    if !is_inside(current_cell.shape, current_cell_coords, X, tol=HOST_CELL_TOL)
                        step = 0.5 * (ξ_trial - ξ)
                        break
                    end
                end

                ξ += step
                X = evaluate(path, cmd, ξ)
                n_bisect_its = floor(Int, log(2, step / HOST_BISECT_TOL)) + 1

                for _ in 1:n_bisect_its
                    step *= 0.5
                    if is_inside(current_cell.shape, current_cell_coords, X, tol=HOST_CELL_TOL)
                        ξ += step
                    else
                        ξ -= step
                    end

                    X = evaluate(path, cmd, ξ)
                end
            else # current_cell is nothing (hole or gap)
                step = 0.5 * (1.0 - ξ)
                ξ += step
                X = evaluate(path, cmd, ξ)
                n_bisect_its = floor(Int, log(2, step / HOST_BISECT_TOL)) + 1

                for _ in 1:n_bisect_its
                    step *= 0.5
                    current_cell = find_elem(X, mesh.elems, mesh._elempartition, tol=HOST_CELL_TOL)
                    if current_cell === nothing
                        ξ += step
                    else
                        ξ -= step
                    end

                    X = evaluate(path, cmd, ξ)
                end
            end

            # Check if end was reached
            if ξ > segment_end - HOST_PATH_TOL
                is_last_segment = true
                is_path_end = is_last_command
            end

            # Getting line cell nodes
            if last_path_node === nothing
                P1 = first_path_node
                push!(mesh.nodes, P1)
            else
                P1 = last_path_node
            end

            if is_last_command && is_last_segment && is_closed
                P2 = first_path_node
            else
                P2 = Node(X)
                push!(mesh.nodes, P2)
            end

            points = _path_line_nodes!(mesh, path, cmd, line_shape, P1, P2, t1=ξ_prev, t2=ξ)

            line_cell = _add_line_cell!(mesh, line_shape, points, line_tag)

            if current_cell !== nothing
                _attach_hosted_line!(
                    mesh,
                    current_cell,
                    line_cell,
                    subpath,
                    P1,
                    P2,
                    is_path_start=is_path_start,
                    is_path_end=is_path_end,
                    is_closed=is_closed,
                )
            end

            # update first and last points
            last_path_node = P2
            is_path_start = false

            is_last_segment && break

            # Preparing for the next iteration
            X = evaluate(path, cmd, ξ + HOST_PATH_TOL)
            current_cell = find_elem(X, mesh.elems, mesh._elempartition, tol=HOST_CELL_TOL, exclude=[current_cell])
            ξ_prev = ξ
            ξ = ξ + HOST_PATH_TOL
        end

    end

end


function gen_conforming_path_insets(mesh::Mesh, subpath::GPath)
    path = subpath.path
    tol = PATH_NODE_TOL

    for cmd in path.cmds
        cmd.key == :M && continue

        X_start, X_end = _path_command_endpoints(path, cmd)
        current_node = _find_mesh_node(mesh, X_start, tol=tol, required=true, label="conforming path start")
        end_node = _find_mesh_node(mesh, X_end, tol=tol, required=true, label="conforming path end")
        target_t = command_parameter(path, cmd, end_node.coord, tol=tol)
        target_t === nothing && error("gen_path_insets: conforming path end does not lie on command")

        while current_node !== end_node
            matches = Tuple{Vector{Node}, Node}[]
            for edge in get_edges(current_node)
                match = _matching_oriented_edge(path, cmd, edge, current_node, target_t, tol=tol)
                match === nothing || push!(matches, match)
            end

            isempty(matches) && error("gen_path_insets: conforming path cannot continue from node $(current_node.id)")
            length(matches) > 1 && error("gen_path_insets: ambiguous conforming path at node $(current_node.id)")

            line_nodes, next_node = matches[1]
            _add_line_cell!(mesh, _line_shape_for_nodes(line_nodes), line_nodes, subpath.tag)
            current_node = next_node
        end
    end
end


function gen_free_path_insets(mesh::Mesh, subpath::GPath)
    path = subpath.path
    line_shape = subpath.shape
    tol = PATH_NODE_TOL

    for cmd in path.cmds
        cmd.key == :M && continue

        X_start, X_end = _path_command_endpoints(path, cmd)
        P1 = _find_or_create_mesh_node!(mesh, X_start, tol=tol)
        P2 = _find_or_create_mesh_node!(mesh, X_end, tol=tol)
        points = _path_line_nodes!(mesh, path, cmd, line_shape, P1, P2)

        _add_line_cell!(mesh, line_shape, points, subpath.tag)
    end
end
