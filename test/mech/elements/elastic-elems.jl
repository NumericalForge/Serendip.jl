using Serendip
using Test

dis = [ -0.012, -0.095 ]

top_node = nothing

for shape in (TRI3, TRI6, QUAD4, QUAD8, QUAD9)
    printstyled(shape.name, color=:cyan); println()
    geo = GeoModel()
    add_block(geo, [0, 0], 1, 1, 0, nx=2, ny=2, shape=shape)
    mesh = Mesh(geo)
    select(mesh, :face, y==0, tag="bottom") # bottom face
    select(mesh, :face, y==1, tag="top") # top face

    mapper = RegionModel(MechBulk, LinearElastic, E=100.0, nu=0.2)

    model = FEModel(mesh, mapper)

    ana = MechAnalysis(model)
    stage = add_stage(ana)
    add_bc(stage, :edge, "bottom", ux=0, uy=0)
    add_bc(stage, :edge, "top", ty=-10.)
    run(ana).successful

    global top_node = select(model, :node, y==1)[1]
    ux = top_node.dofs[:ux].vals[:ux]
    uy = top_node.dofs[:uy].vals[:uy]
    @test [ux, uy] ≈ dis atol=4e-2

    println( get_values(top_node) )

end

#for shape in (TET4, TET10, HEX8, HEX20, HEX27)
for shape in (TET4, TET10, HEX8, HEX20, HEX27)
    printstyled(shape.name, color=:cyan); println()
    geo = GeoModel()
    add_block(geo, [0, 0, 0], 1, 1, 1, nx=2, ny=2, nz=2, shape=shape)
    mesh = Mesh(geo)

    select(mesh, :face, z==0, tag="bottom") # bottom face
    select(mesh, :face, z==1, tag="top") # top face
    select(mesh, :face, x==0, tag="sides") # lateral face
    select(mesh, :face, x==1, tag="sides") # lateral face

    mapper = RegionModel(MechBulk, LinearElastic, E=100.0, nu=0.2)
    model = FEModel(mesh, mapper)

    ana = MechAnalysis(model)
    stage = add_stage(ana)
    add_bc(stage, :face, "bottom", ux=0, uy=0, uz=0)
    add_bc(stage, :face, "sides", ux=0)
    add_bc(stage, :face, "top", tz=-10.)
    run(ana).successful

    # top_node = model.nodes[:(z==1)][1]
    global top_node = select(model, :node, z==1)[1]

    uy = top_node.dofs[:uy].vals[:uy]
    uz = top_node.dofs[:uz].vals[:uz]

    # println( get_values(top_node) )
    @test [uy, uz] ≈ dis atol=1e-2
end
