using Serendip
using Test

@testset "Path modes" begin
    @testset "Invalid path mode" begin
        geo = GeoModel(quiet=true)
        add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0, nx=1, ny=1, shape=:quad4)

        p1 = add_point(geo, [0.0, 0.0, 0.0])
        p2 = add_point(geo, [1.0, 0.0, 0.0])
        l1 = add_line(geo, p1, p2)

        @test_throws Exception add_path(geo, [l1]; mode=:bad)
        @test_throws MethodError add_path(geo, [l1]; embedded=true)
    end

    @testset "Interface mode creates interface elements" begin
        geo = GeoModel(quiet=true)
        add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0, nx=2, ny=2, shape=:quad4)

        p1 = add_point(geo, [0.0, 0.0, 0.0])
        p2 = add_point(geo, [1.0, 1.0, 0.0])
        l1 = add_line(geo, p1, p2)
        add_path(geo, [l1]; mode=:interface, tag="bars", interface_tag="joints")

        mesh = Mesh(geo, quiet=true)
        lines = select(mesh.elems, :line, "bars")
        interfaces = select(mesh.elems, :line_interface, "joints")

        @test length(lines) == 2
        @test length(interfaces) == 2
        @test all(length(joint.couplings) == 2 for joint in interfaces)
        @test all(!line.embedded for line in lines)
    end

    @testset "Embedded mode couples line elements to hosts" begin
        geo = GeoModel(quiet=true)
        add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0, nx=2, ny=2, shape=:quad4)

        p1 = add_point(geo, [0.0, 0.0, 0.0])
        p2 = add_point(geo, [1.0, 1.0, 0.0])
        l1 = add_line(geo, p1, p2)
        add_path(geo, [l1]; mode=:embedded, tag="bars")

        mesh = Mesh(geo, quiet=true)
        lines = select(mesh.elems, :line, "bars")

        @test length(lines) == 2
        @test all(line.embedded for line in lines)
        @test all(length(line.couplings) == 1 for line in lines)
        @test isempty(select(mesh.elems, :line_interface))
    end

    @testset "Conforming mode reuses linear mesh edges" begin
        geo = GeoModel(quiet=true)
        add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0, nx=2, ny=2, shape=:quad4)

        p1 = add_point(geo, [0.0, 0.0, 0.0])
        p2 = add_point(geo, [1.0, 0.0, 0.0])
        l1 = add_line(geo, p1, p2)
        add_path(geo, [l1]; mode=:conforming, tag="bars")

        mesh = Mesh(geo, quiet=true)
        lines = select(mesh.elems, :line, "bars")

        @test length(lines) == 2
        @test all(line.shape.kind == :lin2 for line in lines)
        @test all(isempty(line.couplings) for line in lines)
        @test isempty(select(mesh.elems, :line_interface))
    end

    @testset "Conforming mode reuses quadratic mesh edges" begin
        geo = GeoModel(quiet=true)
        add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0, nx=1, ny=1, shape=:quad8)

        p1 = add_point(geo, [0.0, 0.0, 0.0])
        p2 = add_point(geo, [1.0, 0.0, 0.0])
        l1 = add_line(geo, p1, p2)
        add_path(geo, [l1]; mode=:conforming, tag="bars")

        mesh = Mesh(geo, quiet=true)
        lines = select(mesh.elems, :line, "bars")

        @test length(mesh.nodes) == 8
        @test length(lines) == 1
        @test lines[1].shape.kind == :lin3
        @test lines[1].nodes[3] in mesh.nodes
    end

    @testset "Conforming mode requires existing endpoints" begin
        geo = GeoModel(quiet=true)
        add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0, nx=1, ny=1, shape=:quad4)

        p1 = add_point(geo, [0.25, 0.0, 0.0])
        p2 = add_point(geo, [1.0, 0.0, 0.0])
        l1 = add_line(geo, p1, p2)
        add_path(geo, [l1]; mode=:conforming, tag="bars")

        @test_throws Exception Mesh(geo, quiet=true)
    end

    @testset "Free mode creates standalone path elements" begin
        geo = GeoModel(quiet=true)
        add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0, nx=1, ny=1, shape=:quad4)

        p1 = add_point(geo, [0.0, 0.0, 0.0])
        p2 = add_point(geo, [1.0, 1.0, 0.0])
        l1 = add_line(geo, p1, p2)
        add_path(geo, [l1]; mode=:free, tag="bars", shape=:lin3)

        mesh = Mesh(geo, quiet=true)
        lines = select(mesh.elems, :line, "bars")

        @test length(mesh.nodes) == 5
        @test length(lines) == 1
        @test lines[1].shape.kind == :lin3
        @test isempty(lines[1].couplings)
        @test !lines[1].embedded
        @test isempty(select(mesh.elems, :line_interface))
        @test isempty(select(mesh.elems, :tip))
    end
end
