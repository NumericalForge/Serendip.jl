using Serendip
using Test

@testset "Mapping-level quadrature for generic elements" begin
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

@testset "Mapping-level quadrature for MechBeam" begin
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

@testset "Mapping-level quadrature for 3D MechBeam" begin
    geo = GeoModel()
    add_block(geo, [0.0, 0.0, 0.0], 1.0, 0.0, 0.0; nx=1, shape=:lin3, tag="beam")
    mesh = Mesh(geo, ndim=3, quiet=true)

    mapper = RegionMapper()
    add_mapping(mapper, "beam", MechBeam, LinearElastic; quadrature=(2, 2, 2), E=1.0, nu=0.25, b=0.1, h=0.2)
    model = FEModel(mesh, mapper, quiet=true)
    @test length(model.elems[1].ips) == 8
end

@testset "Mapping-level quadrature for MechShell" begin
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

@testset "Invalid mapping-level quadrature" begin
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
