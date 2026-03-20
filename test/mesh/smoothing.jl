using Serendip
using Test, Statistics

printstyled("\nMesh smoothing\n", color = :blue, bold = true)


function tri_mesh()
    coordinates = [
        0.0 0.0;
        0.5 0.0;
        1.0 0.0;
        0.0 0.5;
        0.5 0.5;
        1.0 0.5;
        0.0 1.0;
        0.5 1.0;
        1.0 1.0;
    ]
    connectivities = [
        [1, 2, 5], [1, 5, 4],
        [2, 3, 6], [2, 6, 5],
        [4, 5, 8], [4, 8, 7],
        [5, 6, 9], [5, 9, 8],
    ]
    cellshapes = fill(:tri3, length(connectivities))
    return Mesh(coordinates, connectivities, cellshapes)
end


function perturb_center(mesh; dx, dy)
    center = only(filter(node ->
        isapprox(node.coord.x, 0.5; atol = 1e-8) &&
        isapprox(node.coord.y, 0.5; atol = 1e-8),
        mesh.nodes,
    ))
    center.coord = center.coord + [dx, dy, 0.0]

    for cell in mesh.elems
        cell.quality = cell_quality(cell)
    end
    mesh.elem_fields["quality"] = Float64[cell.quality for cell in mesh.elems]

    return center
end


println("\nDefault Laplacian:")
mesh = tri_mesh()
perturb_center(mesh, dx=0.18, dy=-0.12)
initial_q = mean(mesh.elem_fields["quality"])
result = smooth(mesh; quiet = false, maxit = 10)
final_q = mean(mesh.elem_fields["quality"])
@test final_q > initial_q

println("\nSmart Laplacian smoothing preserves minimum quality:")
mesh = tri_mesh()
perturb_center(mesh, dx=0.18, dy=-0.12)
initial_qmin = minimum(mesh.elem_fields["quality"])
smooth(mesh; quiet = false, maxit = 10, smart = true)
final_qmin = minimum(mesh.elem_fields["quality"])
@test final_qmin >= initial_qmin

println("\nWeighted Laplacian smoothing updates quality data:")
mesh = tri_mesh()
perturb_center(mesh, dx=-0.18, dy=0.14)
smooth(mesh; quiet = false, maxit = 10, weighted = true)
@test haskey(mesh.elem_fields, "quality")
@test all(isfinite, mesh.elem_fields["quality"])
@test minimum(mesh.elem_fields["quality"]) > 0.0

println("\nBoundary nodes remain fixed when requested:")
mesh = tri_mesh()
perturb_center(mesh, dx=0.2, dy=-0.15)
surf_cells = get_outer_facets(mesh.elems)
surf_nodes, _ = get_patches(surf_cells)
boundary_ids = sort([node.id for node in surf_nodes])
boundary_before = copy(get_coords([mesh.nodes[id] for id in boundary_ids], mesh.ctx.ndim))
smooth(mesh; quiet = false, maxit = 10, fixed_boundary = true)
boundary_after = get_coords([mesh.nodes[id] for id in boundary_ids], mesh.ctx.ndim)
@test boundary_after ≈ boundary_before
