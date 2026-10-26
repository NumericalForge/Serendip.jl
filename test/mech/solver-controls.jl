using Serendip, Test


@testset "mechanical dof scaling" begin
    geo = GeoModel()
    add_block(geo, [0.0, 0.0, 0.0], 2.0, 0.0, 0.0, n=1, shape=:lin2, tag="beam")
    mesh = Mesh(geo, ndim=2)

    mapper = RegionMapper()
    add_map(mapper, "beam", MechBeam, LinearElastic, E=1.0e4, nu=0.0, b=0.1, h=0.2)
    model = FEModel(mesh, mapper)
    ana = MechAnalysis(model)
    stage = add_stage(ana)
    add_bc(stage, :node, x==0.0, ux=0.0, uy=0.0, rz=0.0)

    dofs, _ = Serendip.configure_dofs(model, stage.bcs)
    force_scale = Serendip._mechanical_force_scale(model.elems, dofs)

    for dof in dofs
        if dof.name == :rz
            @test force_scale[dof.eq_id] ≈ 0.5
        else
            @test force_scale[dof.eq_id] ≈ 1.0
        end
    end
end


@testset "force reference" begin
    external_force = [10.0, 20.0]
    internal_force = [9.0, -30.0]
    force_scale = [1.0, 0.5]

    @test Serendip._force_reference(external_force, internal_force, force_scale, 2:2) ≈ 25.0
end


@testset "iteration step control" begin
    @test Serendip._iteration_step_factor(1, 4) ≈ 1.5
    @test Serendip._iteration_step_factor(4, 4) ≈ 1.0
    @test Serendip._iteration_step_factor(16, 4) ≈ 0.5
end
