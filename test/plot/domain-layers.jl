using Serendip
using Test
using LinearAlgebra

function make_square_case()
    geo = GeoModel(quiet=true)
    add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0, nx=3, ny=3, shape=:quad4, tag="solids")
    mesh = Mesh(geo, quiet=true)

    mapper = RegionMapper()
    add_mapping(mapper, "solids", MechSolid, LinearElastic, E=30e6, nu=0.25)

    model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=1.0, quiet=true)
    ana = MechAnalysis(model)
    stage = add_stage(ana)
    add_bc(stage, :node, x == 0.0, ux=0.0)
    add_bc(stage, :node, (x == 0.0, y == 0.0), uy=0.0)
    add_bc(stage, :node, (x==1, y==1), uy=-0.05)
    run(ana, quiet=true)

    return (mesh=mesh, model=model)
end


function make_cube_case()
    geo = GeoModel(quiet=true)
    add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 1.0, nx=3, ny=3, nz=3, tag="solids")
    mesh = Mesh(geo, quiet=true)

    mapper = RegionMapper()
    add_mapping(mapper, "solids", MechSolid, LinearElastic, E=20e6, nu=0.25)

    model = FEModel(mesh, mapper, quiet=true)
    ana = MechAnalysis(model)
    stage = add_stage(ana)
    add_bc(stage, :node, x == 0.0, ux=0.0)
    add_bc(stage, :node, y == 0.0, uy=0.0)
    add_bc(stage, :node, z == 0.0, uz=0.0)
    add_bc(stage, :node, (x==1, y==1, z==1), uz=-0.05)
    run(ana, quiet=true)

    return (mesh=mesh, model=model)
end


function max_layer_offset(plot)
    coords1 = [node.coord for node in sort(plot.layers[1].nodes, by=node -> node.id)]
    coords2 = [node.coord for node in sort(plot.layers[2].nodes, by=node -> node.id)]
    return maximum(norm(c2 - c1) for (c1, c2) in zip(coords1, coords2))
end


case2d = make_square_case()
@test haskey(case2d.model.node_fields, "U")
@test maximum(abs, case2d.model.node_fields["U"]) > 0.0
case2d.mesh.node_fields["temp"] = collect(1.0:length(case2d.mesh.nodes))

plot2d = DomainPlot(size=(5cm, 4cm))
add_plot(plot2d, case2d.mesh, face_color=:aliceblue, view_mode=:surface_with_edges, field="temp", field_kind=:node, colorbar=:left)
add_plot(plot2d, case2d.model, warp=5.0, view_mode=:outline, line_color=:red, field="ux", field_kind=:node, vector_field="U", arrow_color=:red, colorbar=:right)

Serendip.configure!(plot2d)
@test length(plot2d.layers) == 2
@test length(plot2d.colorbars) == 2
@test length(plot2d.left_items) == 1
@test length(plot2d.right_items) == 1
@test minimum(getfield.(plot2d.layers[1].nodes, :id)) == 1
@test minimum(getfield.(plot2d.layers[2].nodes, :id)) == 1
@test max_layer_offset(plot2d) > 0.0
@test plot2d.layers[2].vector_field == "U"
@test maximum(norm.(eachrow(plot2d.layers[2].vector_values))) > 0.0

save(plot2d, "output/domain-layers-2d.pdf")
@test isfile("output/domain-layers-2d.pdf")


case3d = make_cube_case()
@test haskey(case3d.model.node_fields, "U")
@test maximum(abs, case3d.model.node_fields["U"]) > 0.0
case3d.mesh.node_fields["temp"] = collect(1.0:length(case3d.mesh.nodes))

plot3d = DomainPlot(size=(5cm, 4cm))
add_plot(plot3d, case3d.mesh, face_color=:aliceblue, view_mode=:surface_with_edges, field="temp", field_kind=:node, colorbar=:top)
add_plot(plot3d, case3d.model, warp=5.0, view_mode=:outline, line_color=:red, field="ux", field_kind=:node, vector_field="U", arrow_color=:red, colorbar=:bottom)

Serendip.configure!(plot3d)
@test length(plot3d.layers) == 2
@test length(plot3d.colorbars) == 2
@test length(plot3d.top_items) == 1
@test length(plot3d.bottom_items) == 1
@test max_layer_offset(plot3d) > 0.0
@test plot3d.layers[2].vector_field == "U"
@test maximum(norm.(eachrow(plot3d.layers[2].vector_values))) > 0.0

save(plot3d, "output/domain-layers-3d.pdf")
@test isfile("output/domain-layers-3d.pdf")
