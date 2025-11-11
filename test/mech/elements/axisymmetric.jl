using Serendip
using Test

for shape in (TRI3, TRI6, QUAD4, QUAD8)
    # Axisymmetric
    printstyled("$(shape.name)\n", color=:cyan)
    printstyled("axisymmetric\n", color=:cyan)

    geo = GeoModel()
    add_block(geo, [0, 0], 1, 1, 0, nx=4, ny=4, shape=shape, tag="solids")
    mesh = Mesh(geo)

    mapper = RegionModel(MechBulk, LinearElastic, E=100.0, nu=0.2)

    model = FEModel(mesh, mapper, stress_state=:axisymmetric)
    ana   = MechAnalysis(model)

    stage = add_stage(ana)
    add_bc(stage, :face, x==0, ux=0)
    add_bc(stage, :face, y==0, uy=0)
    add_bc(stage, :face, y==1, ty=-10)
    run(ana).successful

    sample_node = select(model, :node, (x==1, y==1))[1]
    uxr = get_dof(sample_node, :ux).vals[:ux]
    uyr = get_dof(sample_node, :uy).vals[:uy]
    # println( get_values(sample_node) )

    # 3D
    printstyled("3d version", color=:cyan); println()

    mesh = revolve(mesh, base=[0,0,0], axis=[0,1,0], n=12)

    model = FEModel(mesh, mapper)
    ana   = MechAnalysis(model)

    stage = add_stage(ana)
    add_bc(stage, :node, (x==0,y==0, z==0), uz=0)
    add_bc(stage, :node, (x==0,y==0), ux=0, uy=0)
    add_bc(stage, :face, y==0, uy=0)
    add_bc(stage, :face, y==1, ty=-10)
    run(ana).successful

    sample_node = select(model, :node, (x==1, y==1))[1]
    ux = get_dof(sample_node, :ux).vals[:ux]
    uy = get_dof(sample_node, :uy).vals[:uy]
    # println( get_values(sample_node) )

    # Verification
    @test [uxr, uyr] ≈ [ux, uy] atol=1e-3

end
