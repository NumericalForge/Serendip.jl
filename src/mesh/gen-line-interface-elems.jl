# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

function gen_insets!(mesh::Mesh, subpaths::Vector{GPath})
    for subpath in subpaths
        gen_insets!(mesh, subpath)
    end
end


function gen_insets!(mesh::Mesh, subpath::GPath)
    path      = subpath.path
    commands  = path.cmds
    is_closed = path.closed
    line_tag  = subpath.tag
    interface_tag = subpath.interface_tag
    tip_tag       = subpath.tip_tag
    tips          = subpath.tips
    is_embedded   = subpath.embedded
    line_shape    = subpath.shape

    # tolerances
    eps_bisect = 1e-9  # for bisection
    eps_path   = 1e-4  # for checking first cell and end of path
    eps_cell   = 1e-9  # for checking if point is inside a cell
    probe_step = 1.0   # step for checking multiple intersections in one cell (0 < λ < 1)

    nline_points = line_shape == LIN2 ? 2 : 3
    joint_shape = POLYVERTEX
    # jntshape = subpath.shape
    # jntshape = subpath.shape==LIN2 ? JLINK2 : JLINK3

    segment_end = 1.0
    first_path_node = Node(path.points[1].coord)
    last_path_node = nothing
    is_path_start = true
    is_path_end = false
    for (cmd_idx, cmd) in enumerate(commands)
        cmd.key == :M && continue

        # find the initial element
        X_start = evaluate(path, cmd, eps_path) # a little bit ahead from 0.0
        current_cell = find_elem(X_start, mesh.elems, mesh._elempartition, tol=eps_cell) # first crossed cell

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
                n_probe_its = round(Int, 1.0 / probe_step)
                for _ in 1:n_probe_its
                    ξ_trial += probe_step
                    ξ_trial > 1.0 && break
                    X = evaluate(path, cmd, ξ_trial)
                    if !is_inside(current_cell.shape, current_cell_coords, X, tol=eps_cell)
                        step = 0.5 * (ξ_trial - ξ)
                        break
                    end
                end

                ξ += step
                X = evaluate(path, cmd, ξ)
                n_bisect_its = floor(Int, log(2, step / eps_bisect)) + 1

                for _ in 1:n_bisect_its
                    step *= 0.5
                    if is_inside(current_cell.shape, current_cell_coords, X, tol=eps_cell)
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
                n_bisect_its = floor(Int, log(2, step / eps_bisect)) + 1

                for _ in 1:n_bisect_its
                    step *= 0.5
                    current_cell = find_elem(X, mesh.elems, mesh._elempartition, tol=eps_cell)
                    if current_cell === nothing
                        ξ += step
                    else
                        ξ -= step
                    end

                    X = evaluate(path, cmd, ξ)
                end
            end

            # Check if end was reached
            if ξ > segment_end - eps_path
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

            if nline_points == 2
                points = [P1, P2]
            else
                X = evaluate(path, cmd, (ξ_prev + ξ) / 2)
                P3 = Node(X)
                push!(mesh.nodes, P3)
                points = [P1, P2, P3]
            end

            # Saving line cell
            line_cell = Cell(line_shape, :line, points, tag=line_tag)
            push!(mesh.elems, line_cell)

            if current_cell !== nothing
                if is_embedded
                    # Set line as embedded
                    line_cell.embedded = true
                    line_cell.couplings = [current_cell]
                else
                    # Generate a continuous joint element
                    joint_points = vcat(current_cell.nodes, line_cell.nodes)
                    joint_cell = Cell(joint_shape, :line_interface, joint_points, tag=interface_tag)
                    push!(mesh.elems, joint_cell)
                    joint_cell.couplings = [current_cell, line_cell]

                    # generate tip joints
                    if !is_closed && is_path_start && tips in (:start, :both)
                        tip = P1
                        tip_joint_nodes = vcat(current_cell.nodes, tip)
                        tip_joint_cell = Cell(POLYVERTEX, :tip, tip_joint_nodes, tag=tip_tag)
                        tip_joint_cell.couplings = joint_cell.couplings
                        push!(mesh.elems, tip_joint_cell)
                    end
                    if !is_closed && is_path_end && tips in (:end, :both)
                        tip = P2
                        tip_joint_nodes = vcat(current_cell.nodes, tip)
                        tip_joint_cell = Cell(POLYVERTEX, :tip, tip_joint_nodes, tag=tip_tag)
                        tip_joint_cell.couplings = joint_cell.couplings
                        push!(mesh.elems, tip_joint_cell)
                    end
                end
                current_cell.crossed = true
            end

            # update first and last points
            last_path_node = P2
            is_path_start = false

            is_last_segment && break

            # Preparing for the next iteration
            X = evaluate(path, cmd, ξ + eps_path)
            current_cell = find_elem(X, mesh.elems, mesh._elempartition, tol=eps_cell, exclude=[current_cell])
            ξ_prev = ξ
            ξ = ξ + eps_path
        end

    end

end
