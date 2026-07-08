using Serendip
using Test

dis = [ -0.012, -0.095 ]

top_node = nothing

@announced_testset "2d elastic elements" begin
    for shape in (:tri3, :tri6, :quad4, :quad8, :quad9)
        @announced_testset "$shape" begin
            geo = GeoModel()
            add_block(geo, [0, 0], 1, 1, 0, nx=2, ny=2, shape=shape)
            mesh = Mesh(geo, sort=false)
            select(mesh, :face, y==0, tag="bottom")
            select(mesh, :face, y==1, tag="top")

            mapper = RegionModel(MechSolid, LinearElastic, E=100.0, nu=0.2)
            model = FEModel(mesh, mapper)
            ana = MechAnalysis(model)
            stage = add_stage(ana)
            add_bc(stage, :edge, "bottom", ux=0, uy=0)
            add_bc(stage, :edge, "top", ty=-10.)
            @test run(ana).successful

            global top_node = select(model, :node, y==1)[1]
            ux = top_node.dofs[:ux].vals[:ux]
            uy = top_node.dofs[:uy].vals[:uy]
            @test [ux, uy] ≈ dis atol=4e-2
        end
    end
end

@announced_testset "3d elastic elements" begin
    for shape in (:tet4, :tet10, :hex8, :hex20, :hex27)
        @announced_testset "$shape" begin
            geo = GeoModel()
            add_block(geo, [0, 0, 0], 1, 1, 1, nx=2, ny=2, nz=2, shape=shape)
            mesh = Mesh(geo, sort=false)

            select(mesh, :face, z==0, tag="bottom")
            select(mesh, :face, z==1, tag="top")
            select(mesh, :face, x==0, tag="sides")
            select(mesh, :face, x==1, tag="sides")

            mapper = RegionModel(MechSolid, LinearElastic, E=100.0, nu=0.2)
            model = FEModel(mesh, mapper)
            ana = MechAnalysis(model)
            stage = add_stage(ana)
            add_bc(stage, :face, "bottom", ux=0, uy=0, uz=0)
            add_bc(stage, :face, "sides", ux=0)
            add_bc(stage, :face, "top", tz=-10.)
            @test run(ana).successful

            global top_node = select(model, :node, z==1)[1]
            uy = top_node.dofs[:uy].vals[:uy]
            uz = top_node.dofs[:uz].vals[:uz]
            @test [uy, uz] ≈ dis atol=1e-2
        end
    end
end
