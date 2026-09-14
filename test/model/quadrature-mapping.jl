using Serendip
using Test
using LinearAlgebra

@announced_testset "Mapping-level quadrature for generic elements" begin
    geo_line = GeoModel()
    add_block(geo_line, [0.0, 0.0], 1.0, 0.0, 0.0; nx=1, shape=:lin3, tag="bar")
    mesh_line = Mesh(geo_line, ndim=2, quiet=true)

    mapper_line = RegionMapper()
    add_mapping(mapper_line, "bar", MechBar, LinearElastic; quadrature=(3,), E=1.0, nu=0.25, A=0.1)
    model_line = FEModel(mesh_line, mapper_line, quiet=true)
    @test length(model_line.elems[1].ips) == 3
    change_quadrature(model_line.elems, (2,))
    @test length(model_line.elems[1].ips) == 2

    geo_quad = GeoModel()
    add_block(geo_quad, [0.0, 0.0], 1.0, 1.0, 0.0; nx=1, ny=1, shape=:quad4, tag="solid")
    mesh_quad = Mesh(geo_quad, quiet=true)

    mapper_quad = RegionMapper()
    add_mapping(mapper_quad, "solid", MechSolid, LinearElastic; quadrature=(2, 2), E=1.0, nu=0.25)
    model_quad = FEModel(mesh_quad, mapper_quad, stress_state=:plane_strain, quiet=true)
    @test length(model_quad.elems[1].ips) == 4
    change_quadrature(model_quad.elems, (3, 3))
    @test length(model_quad.elems[1].ips) == 9

    geo_hex = GeoModel()
    add_block(geo_hex, [0.0, 0.0, 0.0], 1.0, 1.0, 1.0; nx=1, ny=1, nz=1, shape=:hex8, tag="solid")
    mesh_hex = Mesh(geo_hex, quiet=true)

    mapper_hex = RegionMapper()
    add_mapping(mapper_hex, "solid", MechSolid, LinearElastic; quadrature=(2, 2, 2), E=1.0, nu=0.25)
    model_hex = FEModel(mesh_hex, mapper_hex, quiet=true)
    @test length(model_hex.elems[1].ips) == 8
end

@announced_testset "Mapping-level quadrature for MechBeam" begin
    geo = GeoModel()
    add_block(geo, [0.0, 0.0], 1.0, 0.0, 0.0; nx=1, shape=:lin3, tag="beam")
    mesh = Mesh(geo, ndim=2, quiet=true)

    mapper_scalar = RegionMapper()
    add_mapping(mapper_scalar, "beam", MechBeam, LinearElastic; quadrature=2, E=1.0, nu=0.25, b=0.1, h=0.2)
    model_scalar = FEModel(mesh, mapper_scalar, quiet=true)
    @test length(model_scalar.elems[1].ips) == 4

    mapper_single = RegionMapper()
    add_mapping(mapper_single, "beam", MechBeam, LinearElastic; quadrature=(3,), E=1.0, nu=0.25, b=0.1, h=0.2)
    model_single = FEModel(mesh, mapper_single, quiet=true)
    @test length(model_single.elems[1].ips) == 6

    mapper_tuple = RegionMapper()
    add_mapping(mapper_tuple, "beam", MechBeam, LinearElastic; quadrature=(3, 4), E=1.0, nu=0.25, b=0.1, h=0.2)
    model_tuple = FEModel(mesh, mapper_tuple, quiet=true)
    @test length(model_tuple.elems[1].ips) == 12
    change_quadrature(model_tuple.elems, (2, 3))
    @test length(model_tuple.elems[1].ips) == 6
end

@announced_testset "Mapping-level quadrature for 3D MechBeam" begin
    geo = GeoModel()
    add_block(geo, [0.0, 0.0, 0.0], 1.0, 0.0, 0.0; nx=1, shape=:lin3, tag="beam")
    mesh = Mesh(geo, ndim=3, quiet=true)

    mapper = RegionMapper()
    add_mapping(mapper, "beam", MechBeam, LinearElastic; quadrature=(2, 2, 2), E=1.0, nu=0.25, b=0.1, h=0.2)
    model = FEModel(mesh, mapper, quiet=true)
    @test length(model.elems[1].ips) == 8
end

@announced_testset "Mapping-level quadrature for MechShell" begin
    geo_quad = GeoModel()
    add_block(geo_quad, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0; nx=1, ny=1, shape=:quad4, tag="shell")
    mesh_quad = Mesh(geo_quad, ndim=3, quiet=true)

    mapper_quad = RegionMapper()
    add_mapping(mapper_quad, "shell", MechShell, LinearElastic; quadrature=(2, 2, 2), E=1.0, nu=0.25, thickness=0.1)
    model_quad = FEModel(mesh_quad, mapper_quad, quiet=true)
    @test length(model_quad.elems[1].ips) == 8
    change_quadrature(model_quad.elems, (3, 3, 2))
    @test length(model_quad.elems[1].ips) == 18

    geo_tri = GeoModel()
    add_block(geo_tri, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0; nx=1, ny=1, shape=:tri3, tag="shell")
    mesh_tri = Mesh(geo_tri, ndim=3, quiet=true)

    mapper_tri = RegionMapper()
    add_mapping(mapper_tri, "shell", MechShell, LinearElastic; quadrature=(3, 2), E=1.0, nu=0.25, thickness=0.1)
    model_tri = FEModel(mesh_tri, mapper_tri, quiet=true)
    @test length(model_tri.elems[1].ips) == 6
end

@announced_testset "Shell drilling stiffness is independent of thickness quadrature" begin
    for (shape, surface_rule) in ((:quad4, (2, 2)), (:quad8, (3, 3)), (:tri3, 3), (:tri6, 6))
        geo = GeoModel()
        add_block(geo, [0.0, 0.0, 0.0], 2.0, 1.0, 0.0; nx=1, ny=1, shape=shape, tag="shell")
        mesh = Mesh(geo, ndim=3, quiet=true)
        reference_drilling = nothing

        for nth in (2, 3, 4)
            quadrature = surface_rule isa Tuple ? (surface_rule..., nth) : (surface_rule, nth)
            mapper = RegionMapper()
            add_mapping(mapper, "shell", MechShell, LinearElastic;
                quadrature=quadrature, E=100.0, nu=0.25, thickness=0.1, kappa=0.1)
            model = FEModel(mesh, mapper, quiet=true)
            elem = model.elems[1]
            K, _, _ = Serendip.elem_stiffness(elem)
            elem.etype.κ = 0.0
            K_material, _, _ = Serendip.elem_stiffness(elem)
            K_drilling = K - K_material

            # On a flat shell, the drilling-rotation block is κGt ∫ N Nᵀ dA.
            surface_ips = Serendip.get_ip_coords(elem.shape, surface_rule)
            coords = Serendip.get_coords(elem)
            expected = zeros(length(elem.nodes), length(elem.nodes))
            for qp in surface_ips
                N = elem.shape.func(qp.coord)
                J = coords' * elem.shape.deriv(qp.coord)
                area_scale = sqrt(det(J' * J))
                expected += 0.1 * (100.0 / (2 * 1.25)) * 0.1 * qp.w * area_scale * N * N'
            end
            @test K_drilling[6:6:end, 6:6:end] ≈ expected rtol=1e-11 atol=1e-12
            if reference_drilling === nothing
                reference_drilling = K_drilling
            else
                @test K_drilling ≈ reference_drilling rtol=1e-10 atol=1e-12
            end
        end
    end
end

@announced_testset "Invalid mapping-level quadrature" begin
    mapper = RegionMapper()
    @test_throws ErrorException add_mapping(mapper, :all, MechSolid, LinearElastic; quadrature=-1, E=1.0, nu=0.25)
    # @test_throws ErrorException add_mapping(mapper, :all, MechSolid, LinearElastic; quadrature=(2, 0), E=1.0, nu=0.25)
    # add_mapping(mapper, :all, MechSolid, LinearElastic; quadrature=(2, 0), E=1.0, nu=0.25)
    @test_throws ErrorException add_mapping(mapper, :all, MechSolid, LinearElastic; quadrature=(2, 2, 2, 2), E=1.0, nu=0.25)

    geo_tri = GeoModel()
    add_block(geo_tri, [0.0, 0.0], 1.0, 1.0, 0.0; nx=1, ny=1, shape=:tri3, tag="solid")
    mesh_tri = Mesh(geo_tri, quiet=true)
    bad_tri_mapper = RegionMapper()
    add_mapping(bad_tri_mapper, "solid", MechSolid, LinearElastic; quadrature=(2, 2), E=1.0, nu=0.25)
    # @test_throws ErrorException FEModel(mesh_tri, bad_tri_mapper, stress_state=:plane_strain, quiet=true)

    geo_beam3d = GeoModel()
    add_block(geo_beam3d, [0.0, 0.0, 0.0], 1.0, 0.0, 0.0; nx=1, shape=:lin3, tag="beam")
    mesh_beam3d = Mesh(geo_beam3d, ndim=3, quiet=true)
    bad_beam_mapper = RegionMapper()
    add_mapping(bad_beam_mapper, "beam", MechBeam, LinearElastic; quadrature=(2, 2), E=1.0, nu=0.25, b=0.1, h=0.2)
    # @test_throws ErrorException FEModel(mesh_beam3d, bad_beam_mapper, quiet=true)

    geo_shell = GeoModel()
    add_block(geo_shell, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0; nx=1, ny=1, shape=:quad4, tag="shell")
    mesh_shell = Mesh(geo_shell, ndim=3, quiet=true)
    bad_shell_mapper = RegionMapper()
    add_mapping(bad_shell_mapper, "shell", MechShell, LinearElastic; quadrature=(2, 2), E=1.0, nu=0.25, thickness=0.1)
    # @test_throws ErrorException FEModel(mesh_shell, bad_shell_mapper, quiet=true)
end
