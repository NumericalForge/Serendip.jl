using Serendip
using Test

function geo_split_quad_mesh(; tag)
    coords = [
        0.0 0.0
        1.0 0.0
        1.0 0.5
        0.0 0.5
        1.0 1.0
        0.0 1.0
    ]
    return Mesh(coords, [[1, 2, 3, 4], [4, 3, 5, 6]], [:quad4, :quad4], tag=tag, quiet=true)
end


function geo_split_quad8_mesh(; tag)
    coords = [
        0.0 0.0
        1.0 0.0
        1.0 0.5
        0.0 0.5
        1.0 1.0
        0.0 1.0
        0.5 0.0
        1.0 0.25
        0.5 0.5
        0.0 0.25
        1.0 0.75
        0.5 1.0
        0.0 0.75
    ]
    conns = [
        [1, 2, 3, 4, 7, 8, 9, 10],
        [4, 3, 5, 6, 9, 11, 12, 13],
    ]
    return Mesh(coords, conns, [:quad8, :quad8], tag=tag, quiet=true)
end


@testset "Geometry stored meshes" begin
    @testset "Adjacent OCC region conforms to stored mesh boundary" begin
        old = geo_split_quad_mesh(tag="old")
        geo = GeoModel(quiet=true)
        add_mesh(geo, old)
        add_rectangle(geo, [1.0, 0.0, 0.0], 1.0, 1.0, tag="new")

        mesh = Mesh(geo, quiet=true)

        @test length(select(mesh.elems, "old")) == 2
        @test length(select(mesh.elems, "new")) > 0
        @test any(node -> node.coord.x == 1.0 && node.coord.y == 0.5, mesh.nodes)

        interface_nodes = select(mesh, :node, :(x == 1.0))
        interface_keys = Set((node.coord.x, node.coord.y, node.coord.z) for node in interface_nodes)
        @test length(interface_nodes) == length(interface_keys)
    end

    @testset "Quadratic stored mesh constrains with vertex nodes" begin
        old = geo_split_quad8_mesh(tag="old")
        geo = GeoModel(quiet=true)
        add_mesh(geo, old)
        add_rectangle(geo, [1.0, 0.0, 0.0], 1.0, 1.0, tag="new")

        mesh = Mesh(geo, quadratic=true, quiet=true)

        @test length(select(mesh.elems, "old")) == 2
        @test length(select(mesh.elems, "new")) > 0

        interface_nodes = select(mesh, :node, :(x == 1.0))
        interface_y = sort([node.coord.y for node in interface_nodes])
        @test interface_y == [0.0, 0.25, 0.5, 0.75, 1.0]
    end

    @testset "Stored mesh combines with structured block" begin
        geo = GeoModel(quiet=true)
        old = geo_split_quad_mesh(tag="old")
        add_mesh(geo, old)
        add_block(geo, [1.0, 0.0, 0.0], 1.0, 1.0, 0.0, nx=1, ny=2, shape=:quad4, tag="block")

        mesh = Mesh(geo, quiet=true)

        @test length(select(mesh.elems, "old")) == 2
        @test length(select(mesh.elems, "block")) == 2
        @test length(old.nodes) == 6
        @test length(old.elems) == 2
    end

    @testset "Structured block constrains adjacent OCC region" begin
        geo = GeoModel(quiet=true)
        add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0, nx=1, ny=1, shape=:quad4, tag="block")
        add_rectangle(geo, [1.0, 0.0, 0.0], 1.0, 1.0, tag="new")

        mesh = Mesh(geo, quiet=true)

        @test length(select(mesh.elems, "block")) == 1
        @test length(select(mesh.elems, "new")) > 0

        interface_nodes = select(mesh, :node, :(x == 1.0))
        interface_keys = Set((node.coord.x, node.coord.y, node.coord.z) for node in interface_nodes)
        @test length(interface_nodes) == length(interface_keys)
        @test sort([node.coord.y for node in interface_nodes]) == [0.0, 1.0]
    end
end
