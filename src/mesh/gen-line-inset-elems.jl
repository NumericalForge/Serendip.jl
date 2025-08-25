# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

function gen_insets!(mesh::Mesh, subpaths::Vector{GPath})
    for subpath in subpaths
        gen_insets!(mesh, subpath)
    end
end


function gen_insets!(mesh::Mesh, subpath::GPath)
    cmds      = subpath.path.cmds
    closed    = subpath.path.closed
    tag       = subpath.tag
    interface_tag = subpath.interface_tag
    tip_tag   = subpath.tip_tag
    tips      = subpath.tips
    embedded  = subpath.embedded
    shape     = subpath.shape
    path      = subpath.path

    # @show closed

    # TODO: add option: merge_points
    # TODO: check closed when hole

    # tolerances
    ε  = 1e-9  # for bissection
    εn = 1e-4  # for checking first cell and end of path
    εc = 1e-9  # for checking if point is inside a cell
    λ  = 1.0   # step for checking multiple intersections in one cell (0 < λ < 1)

    npoints = subpath.shape==LIN2 ? 2 : 3
    jntshape = POLYVERTEX
    # jntshape = subpath.shape
    # jntshape = subpath.shape==LIN2 ? JLINK2 : JLINK3

    # Initial conditions
    len = 1.0

    firstnode = Node(path.points[1].coord)
    # firstnode = Node(cmds[1].points[1].coord)
    lastnode = nothing
    last_cmd = false
    for (i,cmd) in enumerate(cmds)
        cmd.key==:M && continue

        # find the initial element
        X0    = evaluate(path, cmd, εn) # a little bit ahead from 0.0
        ccell = find_elem(X0, mesh.elems, mesh._elempartition, tol=εc) # the first tresspased cell

        first_segment = true
        last_segment  = false
        last_cmd = i==length(cmds)
        s  = 0.0
        s1 = 0.0

        # Splitting cmd
        k = 0
        while true
            k +=1

            if ccell !== nothing
                ccell_coords = get_coords(ccell)
                # Default step
                step  = 0.5*(1.0-s)

                # Finding step (introduced to catch elements in curved paths)
                str = s     # trial point
                nits = round(Int, 1.0/λ)
                for i in 1:nits
                    str += λ
                    str>1.0 && break
                    X = evaluate(path, cmd, str)
                    if !is_inside(ccell.shape, ccell_coords, X, tol=εc)
                        step = 0.5*(str-s)
                        break
                    end
                end

                s += step
                X  = evaluate(path, cmd, s)
                n  = floor(Int, log(2, step/ε)) + 1  # number of required iterations to find intersection

                for i in 1:n
                    step *= 0.5
                    if is_inside(ccell.shape, ccell_coords, X, tol=εc)
                        s += step
                    else
                        s -= step
                    end

                    X = evaluate(path, cmd, s)
                end
            else # ccell is nothing (hole or gap)
                step  = 0.5*(1.0-s)
                s += step
                X  = evaluate(path, cmd, s)
                n  = floor(Int, log(2, step/ε)) + 1  # number of required iterations to find intersection

                for i in 1:n
                    step *= 0.5
                    ccell = find_elem(X, mesh.elems, mesh._elempartition, tol=εc)
                    if ccell === nothing
                        s += step
                    else
                        s -= step
                    end

                    X = evaluate(path, cmd, s)
                end
            end

            # Check if end was reached
            if s > len - εn
                last_segment = true
            end

            # Getting line cell nodes
            if lastnode===nothing
                P1 = firstnode
                push!(mesh.nodes, P1)
            else
                P1 = lastnode
            end

            if last_cmd && last_segment && closed
                P2 = firstnode
            else
                P2 = Node(X)
                push!(mesh.nodes, P2)
            end

            if npoints==2
                points = [P1, P2]
            else
                X = evaluate(path, cmd, (s1+s)/2)
                P3 = Node(X)
                push!(mesh.nodes, P3)
                points = [P1, P2, P3]
            end

            # Saving line cell
            lcell = Cell(shape, :line, points, tag=tag)
            push!(mesh.elems, lcell)

            if ccell!== nothing
                if embedded
                    # Set line as embedded
                    lcell.embedded = true
                    lcell.couplings = [ ccell ]
                else
                    # Generate a continuous joint element
                    jntpts  = vcat( ccell.nodes, lcell.nodes )
                    jntcell = Cell(jntshape, :line_interface, jntpts, tag=interface_tag)
                    push!(mesh.elems, jntcell)
                    jntcell.couplings = [ccell, lcell]

                    # generate tip joints
                    if first_segment && tips in (:front, :both)
                        tip = P1
                        tipjointnodes = vcat(ccell.nodes, tip)
                        tipjointcell = Cell(POLYVERTEX, :tip, tipjointnodes, tag=tip_tag)
                        tipjointcell.couplings = jntcell.couplings
                        push!(mesh.elems, tipjointcell)
                    end
                    if last_segment && tips in (:end, :both)
                        tip = P2
                        tipjointnodes = vcat(ccell.nodes, tip)
                        tipjointcell = Cell(POLYVERTEX, :tip, tipjointnodes, tag=tip_tag)
                        tipjointcell.couplings = jntcell.couplings
                        push!(mesh.elems, tipjointcell)
                    end
                end
                ccell.crossed = true
            end

            # update first and last points
            if firstnode === nothing;
                firstnode = P1
            end
            lastnode = P2
            first_segment = false

            last_segment && break

            # Preparing for the next iteration
            X = evaluate(path, cmd, s + εn)
            ccell = find_elem(X, mesh.elems, mesh._elempartition, tol=εc, exclude=[ccell])
            s1    = s
            s     = s + εn
        end

    end

end

