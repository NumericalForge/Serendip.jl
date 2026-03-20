# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl


function faces_normal(faces::Vector{Cell}, facetol)
    ndim = 1 + faces[1].shape.ndim
    normals = Vector{Float64}[]

    for face in faces
        C = get_coords(face, ndim)

        # Shift coordinates to avoid singular regressions when the best-fit
        # line/plane crosses the origin.
        if ndim == 2
            C .+= [pi pi^1.1]
        else
            C .+= [pi pi^1.1 pi^1.2]
        end

        I = ones(size(C, 1))
        N = pinv(C) * I
        normalize!(N)

        if all(norm(N - NN) > facetol for NN in normals)
            push!(normals, N)
        end
    end

    return normals
end


mutable struct sNode
    node::Node
    faces::Array{Cell}
    normals
end


function str_histogram(hist::Vector{Int64})
    m = maximum(hist)
    H = m == 0 ? zeros(Int, length(hist)) : round.(Int, hist ./ m * 7)
    chars = [" ", "_", "▁", "▂", "▃", "▄", "▅", "▆", "▇", "█"]
    return "[" * join(hist[i] == 0 ? " " : chars[H[i] + 2] for i in 1:length(H)) * "]"
end


"""
    smooth(mesh; maxit=20, quiet=true, fixed_boundary=false, mintol=1e-2, tol=1e-4,
           facetol=1e-5, binsize=0.05, smart=false, weighted=false)

Smooth a finite element mesh in place using Laplacian smoothing.

By default, the smoother updates each node from the average position of its
patch neighbors. With `weighted=true`, centroidal Laplacian smoothing is used.
With `smart=true`, node moves that worsen the minimum quality in the local
patch are rejected.

The input mesh is mutated and returned.
"""
function smooth(
    mesh::Mesh;
    maxit::Int64 = 20,
    quiet = true,
    fixed_boundary = false,
    mintol::Float64 = 1e-2,
    tol::Float64 = 1e-4,
    facetol::Float64 = 1e-5,
    binsize::Float64 = 0.05,
    smart = false,
    weighted = false,
)
    quiet || printstyled(
        "Mesh $(smart ? "smart-" : "")$(weighted ? "weighted-" : "")Laplacian smoothing:\n",
        bold = true,
        color = :cyan,
    )

    ndim = mesh.ctx.ndim

    nodes, patches = get_patches(mesh)

    patch_nodes = Vector{Node}[]
    for (node, patch) in zip(nodes, patches)
        curr_patch_nodes = get_nodes(patch)
        idx = findfirst(p -> hash(p) == hash(node), curr_patch_nodes)
        idx !== nothing && splice!(curr_patch_nodes, idx)
        push!(patch_nodes, curr_patch_nodes)
    end

    surf_cells = get_outer_facets(mesh.elems)
    surf_nodes, surf_patches = get_patches(surf_cells)
    border_nodes = [sNode(node, patch, nothing) for (node, patch) in zip(surf_nodes, surf_patches)]

    for snode in border_nodes
        snode.normals = faces_normal(snode.faces, facetol)
    end

    in_border = falses(length(mesh.nodes))
    border_idxs = [snode.node.id for snode in border_nodes]
    in_border[border_idxs] .= true

    map_pn = zeros(Int, length(mesh.nodes))
    for (i, snode) in enumerate(border_nodes)
        map_pn[snode.node.id] = i
    end

    Q = Float64[cell_quality(c) for c in mesh.elems]
    q = mean(Q)
    qmin = minimum(Q)
    qmax = maximum(Q)
    dev = stdm(Q, q)
    q1, q2, q3 = quantile(Q, [0.25, 0.5, 0.75])
    hist = fit(Histogram, Q, 0.0:binsize:1.0, closed = :right).weights

    mesh.elem_fields["quality"] = Q

    quiet || @printf(
        "%4s  %5s  %5s  %5s  %5s  %5s  %5s  %7s  %9s  %10s\n",
        "it", "qmin", "q1", "q2", "q3", "qmax", "qavg", "sdev", "time", "histogram (0:$binsize:1]",
    )
    quiet || @printf(
        "%4d  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %7.5f  %9s",
        0, qmin, q1, q2, q3, qmax, q, dev, "-",
    )
    quiet || println("  ", str_histogram(hist))

    for i in 1:maxit
        sw = StopWatch()

        for (node, patch, local_nodes) in zip(nodes, patches, patch_nodes)
            id = node.id
            X0 = [node.coord.x, node.coord.y, node.coord.z][1:ndim]

            if in_border[id]
                fixed_boundary && continue
                snode = border_nodes[map_pn[id]]
                normals = snode.normals
                nnorm = length(normals)

                (nnorm == 1 && ndim == 2) || (nnorm in (1, 2) && ndim == 3) || continue
            end

            if weighted
                Xcs = [mean(get_coords(c, ndim), dims = 1)[1:ndim] for c in patch]
                Vs = [abs(cell_extent(c)) for c in patch]
                X = X0 + sum(Vs[j] * (Xcs[j] - X0) for j in 1:length(patch)) / sum(Vs)
            else
                X = mean(get_coords(local_nodes), dims = 1)[1:ndim]
            end

            if in_border[id]
                ΔX = X - X0
                if nnorm == 1
                    n1 = normals[1]
                    X = X0 + ΔX - dot(ΔX, n1) * n1
                else
                    n3 = normalize(cross(normals[1], normals[2]))
                    X = X0 + dot(ΔX, n3) * n3
                end
            end

            node.coord = X

            if smart
                patch_qmin0 = minimum(c.quality for c in patch)
                patch_q = [cell_quality(c) for c in patch]
                patch_qmin = minimum(patch_q)

                if patch_qmin < patch_qmin0
                    node.coord = [X0; zeros(3 - length(X0))]
                else
                    for (cell, quality) in zip(patch, patch_q)
                        cell.quality = quality
                    end
                end
            end
        end

        if !smart
            for c in mesh.elems
                c.quality = cell_quality(c)
            end
        end

        Q = Float64[c.quality for c in mesh.elems]
        new_q = mean(Q)
        new_qmin = minimum(Q)

        mesh.elem_fields["quality"] = Q

        Δq = abs(q - new_q)
        Δqmin = new_qmin - qmin

        q = new_q
        qmin = new_qmin
        qmax = maximum(Q)
        dev = stdm(Q, q)
        q1, q2, q3 = quantile(Q, [0.25, 0.5, 0.75])
        hist = fit(Histogram, Q, 0.0:binsize:1.0, closed = :right).weights

        quiet || @printf(
            "%4d  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %7.5f  %9s",
            i, qmin, q1, q2, q3, qmax, q, dev, see(sw, format = :ms),
        )
        quiet || println("  ", str_histogram(hist))

        if Δq < tol && Δqmin < mintol && i > 2
            break
        end
    end

    n_bad_cells = count(quality -> quality <= 0, Q)
    n_bad_cells > 0 && @warn "smooth: $n_bad_cells invalid cells obtained"

    mesh.elem_fields["quality"] = Q

    return mesh
end
