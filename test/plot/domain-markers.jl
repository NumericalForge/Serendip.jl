using Serendip
using Test

isdir("output") || mkdir("output")

function marker_model_2d()
    geo = GeoModel(quiet=true)
    add_block(geo, [0, 0, 0], 2, 2, 0, nx=2, ny=2, shape=:quad4, tag="solids")
    mesh = Mesh(geo, quiet=true)

    mapper = RegionMapper()
    add_mapping(mapper, "solids", MechSolid, LinearElastic, E=1e4, nu=0.25)
    model = FEModel(mesh, mapper, stress_state=:plane_strain, quiet=true)
    ana = MechAnalysis(model; outdir="output", outkey="domain-markers")
    stage = add_stage(ana)

    add_dof(select(model.nodes, [1.0, 1.0, 0.0])[1], :rz, :mz)

    add_bc(stage, :node, [0.0, 0.0, 0.0], ux=0.0, uy=0.0)
    add_bc(stage, :node, [1.0, 0.0, 0.0], ux=0.0)
    add_bc(stage, :node, [2.0, 0.0, 0.0], ux=1.0)
    add_bc(stage, :node, [0.0, 1.0, 0.0], fy=1.0)
    add_bc(stage, :node, [1.0, 1.0, 0.0], rz=0.0)
    add_bc(stage, :face, y==2.0, ty=-1.0)

    return model, stage
end

function node_id(model, X)
    return select(model.nodes, X)[1].id
end

model, stage = marker_model_2d()

@test DomainPlot(model, quiet=true).mark == :none
@test DomainPlot(model, mark=:auto, mark_size=3.0, quiet=true).mark_size == 3.0
@test_throws ArgumentError DomainPlot(model, mark=:unknown, quiet=true)
@test_throws ArgumentError DomainPlot(model, mark=:circle, quiet=true)
@test_throws ArgumentError DomainPlot(model, mark=:square, quiet=true)
@test_throws ArgumentError DomainPlot(model, mark_size=0.0, quiet=true)
@test_throws ArgumentError DomainPlot(model, mark=:support, quiet=true)

kinds = Serendip._domain_resolve_node_marker_kinds(stage, model.ctx.ndim)
@test kinds[node_id(model, [0.0, 0.0, 0.0])] == :full_translation
@test kinds[node_id(model, [1.0, 0.0, 0.0])] == :zero_translation
@test kinds[node_id(model, [2.0, 0.0, 0.0])] == :essential
@test kinds[node_id(model, [0.0, 1.0, 0.0])] == :natural
@test kinds[node_id(model, [1.0, 1.0, 0.0])] == :clamp
@test kinds[node_id(model, [1.0, 2.0, 0.0])] == :natural

regular = DomainPlot(model, mark=:auto, quiet=true)
auto = DomainPlot(model, mark=:auto, stage=stage, quiet=true)
support = DomainPlot(model, mark=:support, stage=stage, quiet=true)

save(regular, "output/domain-markers-regular.pdf")
save(auto, "output/domain-markers-auto.pdf")
save(support, "output/domain-markers-support.pdf")

@test isfile("output/domain-markers-regular.pdf")
@test isfile("output/domain-markers-auto.pdf")
@test isfile("output/domain-markers-support.pdf")

geo3d = GeoModel(quiet=true)
add_block(geo3d, [0, 0, 0], 1, 1, 1, nx=1, ny=1, nz=1, tag="solids")
mesh3d = Mesh(geo3d, quiet=true)
mapper3d = RegionMapper()
add_mapping(mapper3d, "solids", MechSolid, LinearElastic, E=1e4, nu=0.25)
model3d = FEModel(mesh3d, mapper3d, quiet=true)
ana3d = MechAnalysis(model3d; outdir="output", outkey="domain-markers-3d")
stage3d = add_stage(ana3d)
add_bc(stage3d, :face, z==0.0, uz=0.0)

plot3d = DomainPlot(model3d, mark=:auto, stage=stage3d, quiet=true)
save(plot3d, "output/domain-markers-3d.pdf")
@test isfile("output/domain-markers-3d.pdf")
