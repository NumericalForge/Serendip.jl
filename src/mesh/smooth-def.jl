# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl


function _def_basic_coords(shape::CellShape)
    if shape == TRI3
        a = 2.0 / 3.0^0.25
        return [0.0 0.0; a 0.0; a / 2.0 a / 2.0 * sqrt(3.0)]
    elseif shape == TRI6
        a = 2.0 / 3.0^0.25
        h = a / 2.0 * sqrt(3.0)
        return [0.0 0.0; a 0.0; a / 2.0 h; a / 2.0 0.0; 0.75 * a h / 2.0; 0.25 * a h / 2.0]
    elseif shape == QUAD4
        return [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
    elseif shape == QUAD8
        return [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0; 0.5 0.0; 1.0 0.5; 0.5 1.0; 0.0 0.5]
    elseif shape == QUAD9
        return [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0; 0.5 0.0; 1.0 0.5; 0.5 1.0; 0.0 0.5; 0.5 0.5]
    elseif shape == QUAD12
        return [
            0.0 0.0
            1.0 0.0
            1.0 1.0
            0.0 1.0
            1.0 / 3.0 0.0
            1.0 1.0 / 3.0
            2.0 / 3.0 1.0
            0.0 2.0 / 3.0
            2.0 / 3.0 0.0
            1.0 2.0 / 3.0
            1.0 / 3.0 1.0
            0.0 1.0 / 3.0
        ]
    elseif shape == TET4
        a = (6.0 * sqrt(2.0))^(1.0 / 3.0)
        return [
            0.0 0.0 0.0
            a 0.0 0.0
            a / 2.0 sqrt(3.0) / 2.0 * a 0.0
            a / 2.0 sqrt(3.0) / 6.0 * a sqrt(6.0) / 3.0 * a
        ]
    elseif shape == HEX8
        return [
            0.0 0.0 0.0
            1.0 0.0 0.0
            1.0 1.0 0.0
            0.0 1.0 0.0
            0.0 0.0 1.0
            1.0 0.0 1.0
            1.0 1.0 1.0
            0.0 1.0 1.0
        ]
    elseif shape == HEX20
        return [
            0.0 0.0 0.0
            1.0 0.0 0.0
            1.0 1.0 0.0
            0.0 1.0 0.0
            0.0 0.0 1.0
            1.0 0.0 1.0
            1.0 1.0 1.0
            0.0 1.0 1.0
            0.5 0.0 0.0
            1.0 0.5 0.0
            0.5 1.0 0.0
            0.0 0.5 0.0
            0.5 0.0 1.0
            1.0 0.5 1.0
            0.5 1.0 1.0
            0.0 0.5 1.0
            0.0 0.0 0.5
            1.0 0.0 0.5
            1.0 1.0 0.5
            0.0 1.0 0.5
        ]
    elseif shape == WED6
        a = (4.0 / sqrt(3.0))^(1.0 / 3.0)
        return [
            0.0 0.0 0.0
            a 0.0 0.0
            a / 2.0 a / 2.0 * sqrt(3.0) 0.0
            0.0 0.0 a
            a 0.0 a
            a / 2.0 a / 2.0 * sqrt(3.0) a
        ]
    elseif shape == WED15
        a = (4.0 / sqrt(3.0))^(1.0 / 3.0)
        b = a / 2.0 * sqrt(3.0)
        return [
            0.0 0.0 0.0
            a 0.0 0.0
            a / 2.0 b 0.0
            0.0 0.0 a
            a 0.0 a
            a / 2.0 b a
            a / 2.0 0.0 0.0
            0.75 * a b / 2.0 0.0
            0.25 * a b / 2.0 0.0
            a / 2.0 0.0 a
            0.75 * a b / 2.0 a
            0.25 * a b / 2.0 a
            0.0 0.0 a / 2.0
            a 0.0 a / 2.0
            a / 2.0 b a / 2.0
        ]
    elseif shape == PYR5
        a = sqrt(2.0) / 9.0
        return [
            0.0 0.0 0.0
            a 0.0 0.0
            a a 0.0
            0.0 a 0.0
            a / 2.0 a / 2.0 sqrt(2.0) / 2.0 * a
        ]
    end

    error("deformation_smooth: no reference coordinates for shape $(shape.kind)")
end


function _def_matrix_D(E::Float64, nu::Float64)
    c = E / ((1.0 + nu) * (1.0 - 2.0 * nu))
    return [
        c * (1.0 - nu) c * nu c * nu 0.0 0.0 0.0
        c * nu c * (1.0 - nu) c * nu 0.0 0.0 0.0
        c * nu c * nu c * (1.0 - nu) 0.0 0.0 0.0
        0.0 0.0 0.0 c * (1.0 - 2.0 * nu) 0.0 0.0
        0.0 0.0 0.0 0.0 c * (1.0 - 2.0 * nu) 0.0
        0.0 0.0 0.0 0.0 0.0 c * (1.0 - 2.0 * nu)
    ]
end


function _def_matrix_B(ndim::Int, dNdX::Matx, detJ::Float64, B::Matx)
    nnodes = size(dNdX, 1)
    sqr2 = sqrt(2.0)
    B .= 0.0

    if ndim == 2
        for i in 1:nnodes
            j = i - 1
            B[1, 1 + j * ndim] = dNdX[i, 1]
            B[2, 2 + j * ndim] = dNdX[i, 2]
            B[4, 1 + j * ndim] = dNdX[i, 2] / sqr2
            B[4, 2 + j * ndim] = dNdX[i, 1] / sqr2
        end
    else
        for i in 1:nnodes
            dNdx = dNdX[i, 1]
            dNdy = dNdX[i, 2]
            dNdz = dNdX[i, 3]
            j = i - 1
            B[1, 1 + j * ndim] = dNdx
            B[2, 2 + j * ndim] = dNdy
            B[3, 3 + j * ndim] = dNdz
            B[4, 1 + j * ndim] = dNdy / sqr2
            B[4, 2 + j * ndim] = dNdx / sqr2
            B[5, 2 + j * ndim] = dNdz / sqr2
            B[5, 3 + j * ndim] = dNdy / sqr2
            B[6, 1 + j * ndim] = dNdz / sqr2
            B[6, 3 + j * ndim] = dNdx / sqr2
        end
    end

    return detJ
end


function _def_matrix_K(cell::Cell, ndim::Int, E::Float64, nu::Float64)
    nnodes = length(cell.nodes)
    C = get_coords(cell.nodes, ndim)
    K = zeros(nnodes * ndim, nnodes * ndim)
    B = zeros(6, nnodes * ndim)
    DB = Array{Float64}(undef, 6, nnodes * ndim)
    J = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    D = _def_matrix_D(E, nu)

    for ip in get_ip_coords(cell.shape)
        dNdR = cell.shape.deriv(ip.coord)
        @mul J = C' * dNdR
        @mul dNdX = dNdR * inv(J)
        detJ = det(J)
        _def_matrix_B(ndim, dNdX, detJ, B)

        coef = detJ * ip.w
        @mul DB = D * B
        @mul K += coef * B' * DB
    end

    return K
end


function _def_dof_map(cell::Cell)
    ndim = cell.shape.ndim
    map = Int[]
    for node in cell.nodes
        for i in 1:ndim
            push!(map, (node.id - 1) * ndim + i)
        end
    end
    return map
end


function _def_global_stiffness(mesh::Mesh, E::Float64, nu::Float64, A)
    ndim = mesh.ctx.ndim
    ndof = ndim * length(mesh.nodes)
    R, C, V = Int[], Int[], Float64[]

    for cell in mesh.elems
        Ke = _def_matrix_K(cell, ndim, E, nu)
        map = _def_dof_map(cell)
        nr, nc = size(Ke)
        for i in 1:nr
            for j in 1:nc
                push!(R, map[i])
                push!(C, map[j])
                push!(V, Ke[i, j])
            end
        end
    end

    nbc = size(A, 1)
    for (i, j, val) in zip(findnz(A)...)
        push!(R, ndof + i)
        push!(C, j)
        push!(V, val)
        push!(R, j)
        push!(C, ndof + i)
        push!(V, val)
    end

    return sparse(R, C, V, ndof + nbc, ndof + nbc)
end


function _def_condition_functions(conds)
    conds === nothing && return Function[]
    items = conds isa Tuple || conds isa AbstractVector ? conds : (conds,)
    functions = Function[]

    for cond in items
        fn = eval(quote
            (x, y, z) -> ($cond)
        end)
        push!(functions, fn)
    end

    return functions
end


function _def_constraints(mesh::Mesh, fixed_boundary::Bool, conds, facetol)
    ndim = mesh.ctx.ndim
    nnodes = length(mesh.nodes)
    surf_cells = get_outer_facets(mesh.elems)
    surf_nodes, surf_patches = get_patches(surf_cells)
    border_nodes = [sNode(node, patch, nothing) for (node, patch) in zip(surf_nodes, surf_patches)]
    cond_functions = _def_condition_functions(conds)

    nconstraints = 0
    for snode in border_nodes
        if !isempty(cond_functions)
            p = snode.node
            if any(fn(p.coord.x, p.coord.y, p.coord.z) for fn in cond_functions)
                snode.normals = Vector{Float64}[]
                nconstraints += ndim
                continue
            end
        end

        snode.normals = faces_normal(snode.faces, facetol)
        nnorm = length(snode.normals)
        nconstraints += nnorm in (1, 2) ? nnorm : ndim
    end

    R, C, V = Int[], Int[], Float64[]
    if fixed_boundary
        for (i, snode) in enumerate(border_nodes)
            for j in 1:ndim
                push!(R, (i - 1) * ndim + j)
                push!(C, (snode.node.id - 1) * ndim + j)
                push!(V, 1.0)
            end
        end
        return sparse(R, C, V, length(border_nodes) * ndim, nnodes * ndim)
    end

    row = 0
    for snode in border_nodes
        basecol = (snode.node.id - 1) * ndim
        normals = snode.normals
        nnorm = length(normals)

        if nnorm in (1, 2)
            for normal in normals
                row += 1
                for j in 1:ndim
                    push!(R, row)
                    push!(C, basecol + j)
                    push!(V, normal[j])
                end
            end
        else
            for j in 1:ndim
                push!(R, row + j)
                push!(C, basecol + j)
                push!(V, 1.0)
            end
            row += ndim
        end
    end

    return sparse(R, C, V, nconstraints, nnodes * ndim)
end


function _def_rigid_transform(source::Matrix{Float64}, target::Matrix{Float64}, pindexes::Vector{Int}=Int[])
    A = copy(source)
    B = copy(target)
    cA = mean(A, dims=1)
    cB = mean(B, dims=1)

    A .-= cA
    B .-= cB

    for i in pindexes
        A[i, :] .*= 10.0
    end

    U, _, V = svd(A' * B)
    R = V * U'
    if det(R) < 0.0
        V[:, end] .*= -1.0
        R = V * U'
    end

    T = cB - cA * R'
    return R, T
end


function _def_force_bc(mesh::Mesh, E::Float64, nu::Float64, alpha::Float64, extended::Bool)
    nnodes = length(mesh.nodes)
    ndim = mesh.ctx.ndim
    Fbc = zeros(nnodes * ndim)
    qmin = extended ? minimum(cell.quality for cell in mesh.elems) : 0.0

    for cell in mesh.elems
        ncell_nodes = length(cell.nodes)
        C0 = get_coords(cell.nodes, ndim)
        extent = abs(cell_extent(cell))
        scale = extent^(1.0 / ndim)
        extended && (scale *= alpha)

        reference = _def_basic_coords(cell.shape) * scale
        R, d = _def_rigid_transform(reference, C0)
        C1 = reference * R' .+ repeat(d, ncell_nodes, 1)
        U = vec((C1 - C0)')
        K = _def_matrix_K(cell, ndim, E, nu)

        if extended
            F = qmin < 1.0 ? K * U * ((1.0 - cell.quality) / (1.0 - qmin)) : zeros(length(U))
        else
            F = K * U
        end

        for (i, node) in enumerate(cell.nodes)
            for j in 1:ndim
                Fbc[(node.id - 1) * ndim + j] += F[(i - 1) * ndim + j]
            end
        end
    end

    return Fbc
end


"""
    deformation_smooth(mesh; maxit=10, quiet=true, fixed_boundary=false, mintol=2e-2,
                       tol=1e-3, facetol=1e-4, binsize=0.05, smart=false,
                       alpha=1.0, extended=false, conds=nothing)

Smooth a finite element mesh in place using a deformation-based physics
analogy. Boundary constraints are imposed with Lagrange multipliers.
"""
function deformation_smooth(
    mesh::Mesh;
    maxit::Int64 = 10,
    quiet = true,
    fixed_boundary = false,
    mintol::Float64 = 2e-2,
    tol::Float64 = 1e-3,
    facetol::Float64 = 1e-4,
    binsize::Float64 = 0.05,
    smart = false,
    alpha::Float64 = 1.0,
    extended = false,
    conds = nothing,
)

    quiet || printstyled(
        "Mesh $(smart ? "smart-" : "")deformation smoothing:\n",
        bold = true,
        color = :cyan,
    )

    for cell in mesh.elems
        cell.role == :solid ||
            error("deformation_smooth: cells with role $(repr(cell.role)) are not allowed for smoothing: $(cell.shape.kind)")
    end

    E = 1.0
    nu = 0.0
    ndim = mesh.ctx.ndim
    nnodes = length(mesh.nodes)
    _, patches = get_patches(mesh)

    for cell in mesh.elems
        cell.quality = cell_quality(cell)
    end

    Q = Float64[cell.quality for cell in mesh.elems]
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

    A = _def_constraints(mesh, fixed_boundary, conds, facetol)
    nbc = size(A, 1)

    for i in 1:maxit
        sw = StopWatch()
        F = _def_force_bc(mesh, E, nu, alpha, extended)

        if ndim == 2
            mesh.node_fields["forces"] = [reshape(F, ndim, nnodes)' zeros(nnodes)]
        else
            mesh.node_fields["forces"] = reshape(F, ndim, nnodes)'
        end

        F = vcat(F, zeros(nbc))
        K = _def_global_stiffness(mesh, E, nu, A)
        U = lu(K) \ F

        for node in mesh.nodes
            X0 = node.coord
            pos = (node.id - 1) * ndim + 1
            X = X0 + extend!(vec(U[pos:pos + ndim - 1]), 3)
            node.coord = X

            if smart
                patch = patches[node.id]
                patch_qmin0 = minimum(cell.quality for cell in patch)
                patch_q = [cell_quality(cell) for cell in patch]
                patch_qmin = minimum(patch_q)

                if patch_qmin < 0.8 * patch_qmin0
                    node.coord = X0
                else
                    for (cell, quality) in zip(patch, patch_q)
                        cell.quality = quality
                    end
                end
            end
        end

        if !smart
            for cell in mesh.elems
                cell.quality = cell_quality(cell)
            end
        end

        Q = Float64[cell.quality for cell in mesh.elems]
        new_q = mean(Q)
        new_qmin = minimum(Q)
        mesh.elem_fields["quality"] = Q

        delta_q = abs(q - new_q)
        delta_qmin = new_qmin - qmin

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

        if delta_q < tol && delta_qmin < mintol && i > 1
            break
        end
    end

    n_bad_cells = count(quality -> quality <= 0, Q)
    n_bad_cells > 0 && @warn "deformation_smooth: $n_bad_cells invalid cells obtained"

    mesh.node_fields["forces"] = zeros(length(mesh.nodes), 3)
    mesh.elem_fields["quality"] = Q

    return mesh
end
