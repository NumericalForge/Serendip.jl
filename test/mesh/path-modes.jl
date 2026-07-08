using Serendip
using Test
using LinearAlgebra: norm

@testset "Path modes" begin
    @testset "Path closure policy is consistent" begin
        path_coords_auto = Path([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]; closed=:auto)
        path_coords_open = Path([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]; closed=false)
        path_coords_forced = Path([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]; closed=true)

        p1 = Point([0.0, 0.0, 0.0])
        p2 = Point([1.0, 0.0, 0.0])
        p3 = Point([0.0, 0.0, 0.0])
        edge1 = Edge(-1, "Line", [p1, p2])
        edge2 = Edge(-1, "Line", [p2, p3])
        path_edges_auto = Path([edge1, edge2]; closed=:auto)
        path_edges_open = Path([edge1, edge2]; closed=false)

        @test path_coords_auto.closed
        @test !path_coords_open.closed
        @test path_coords_forced.closed
        @test length(path_coords_forced.cmds) == 4
        @test path_coords_forced.cmds[end].key == :L
        @test path_coords_forced.cmds[end].idxs == [3, 1]
        @test path_edges_auto.closed
        @test !path_edges_open.closed
    end

    @testset "Invalid path mode" begin
        geo = GeoModel(quiet=true)
        add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0, nx=1, ny=1, shape=:quad4)

        p1 = add_point(geo, [0.0, 0.0, 0.0])
        p2 = add_point(geo, [1.0, 0.0, 0.0])
        l1 = add_line(geo, p1, p2)

        @test_throws Exception add_path(geo, [l1]; mode=:bad)
        @test_throws MethodError add_path(geo, [l1]; embedded=true)
        @test_throws Exception add_path(geo, [l1]; mode=:free, n=0)
        @test_throws Exception add_path(geo, [l1]; mode=:embedded, n=2)
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
        add_path(geo, [l1]; mode=:free, tag="bars")

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

    @testset "Free mode subdivides each command with scalar n" begin
        geo = GeoModel(quiet=true)
        add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0, nx=1, ny=1, shape=:quad4)

        add_path(geo, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]; mode=:free, tag="bars", quadratic=false, n=3)

        mesh = Mesh(geo, quiet=true)
        lines = select(mesh.elems, :line, "bars")

        @test length(lines) == 6
        @test all(line.shape.kind == :lin2 for line in lines)
    end

    @testset "Free mode adds implicit closing segment when forced closed" begin
        geo = GeoModel(quiet=true)
        add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0, nx=1, ny=1, shape=:quad4)

        add_path(geo, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0];
            mode=:free, tag="bars", quadratic=false, closed=true)

        mesh = Mesh(geo, quiet=true)
        lines = select(mesh.elems, :line, "bars")

        @test length(lines) == 3
        @test all(line.shape.kind == :lin2 for line in lines)
    end

    @testset "Free mode subdivides higher-order line shapes" begin
        geo = GeoModel(quiet=true)
        add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0, nx=1, ny=1, shape=:quad4)

        p1 = add_point(geo, [0.0, 0.0, 0.0])
        p2 = add_point(geo, [1.0, 0.0, 0.0])
        l1 = add_line(geo, p1, p2)
        add_path(geo, [l1]; mode=:free, tag="bars", n=2)

        mesh = Mesh(geo, quiet=true)
        lines = select(mesh.elems, :line, "bars")
        endpoints1 = Set(lines[1].nodes[1:2])
        endpoints2 = Set(lines[2].nodes[1:2])

        @test length(lines) == 2
        @test all(line.shape.kind == :lin3 for line in lines)
        @test length(intersect(endpoints1, endpoints2)) == 1
    end

    @testset "Free mode uses arc-length subdivision on curved paths" begin
        geo = GeoModel(quiet=true)
        add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0, nx=1, ny=1, shape=:quad4)

        pc = add_point(geo, [0.0, 0.0, 0.0])
        p1 = add_point(geo, [1.0, 0.0, 0.0])
        p2 = add_point(geo, [0.0, 1.0, 0.0])
        arc = add_circle_arc(geo, p1, pc, p2)
        add_path(geo, [arc]; mode=:free, tag="bars", quadratic=false, n=4)

        mesh = Mesh(geo, quiet=true)
        lines = select(mesh.elems, :line, "bars")
        lengths = [norm(line.nodes[2].coord - line.nodes[1].coord) for line in lines]

        @test length(lines) == 4
        @test maximum(lengths) - minimum(lengths) < 1e-3
    end
end
