using Serendip
using Test

function quad_mesh(x1, x2; tag)
    coords = [
        x1 0.0
        x2 0.0
        x2 1.0
        x1 1.0
    ]
    return Mesh(coords, [[1, 2, 3, 4]], [:quad4], tag=tag, quiet=true)
end

function split_quad_mesh(x1, x2; tag)
    coords = [
        x1 0.0
        x2 0.0
        x2 0.5
        x1 0.5
        x2 1.0
        x1 1.0
    ]
    return Mesh(coords, [[1, 2, 3, 4], [4, 3, 5, 6]], [:quad4, :quad4], tag=tag, quiet=true)
end

function hex_mesh(x1, x2; tag)
    coords = [
        x1 0.0 0.0
        x2 0.0 0.0
        x2 1.0 0.0
        x1 1.0 0.0
        x1 0.0 1.0
        x2 0.0 1.0
        x2 1.0 1.0
        x1 1.0 1.0
    ]
    return Mesh(coords, [[1, 2, 3, 4, 5, 6, 7, 8]], [:hex8], tag=tag, quiet=true)
end

@testset "Mesh joining" begin
    @testset "add_mesh mutates one mesh" begin
        left = quad_mesh(0.0, 1.0, tag="left")
        right = quad_mesh(1.0, 2.0, tag="right")

        @test_throws MethodError add_mesh(left, right, quad_mesh(2.0, 3.0, tag="extra"))

        result = add_mesh(left, right)

        @test result === left
        @test length(left.nodes) == 6
        @test length(left.elems) == 2
        @test sort(getfield.(left.elems, :tag)) == ["left", "right"]
        @test length(right.nodes) == 4
        @test length(right.elems) == 1
    end

    @testset "join_meshes is non-mutating and variadic" begin
        m1 = quad_mesh(0.0, 1.0, tag="one")
        m2 = quad_mesh(1.0, 2.0, tag="two")
        m3 = quad_mesh(2.0, 3.0, tag="three")

        joined = join_meshes(m1, m2, m3)

        @test length(joined.nodes) == 8
        @test length(joined.elems) == 3
        @test sort(getfield.(joined.elems, :tag)) == ["one", "three", "two"]
        @test length(m1.nodes) == 4
        @test length(m2.nodes) == 4
        @test length(m3.nodes) == 4
        @test_throws Exception join_meshes(m1)
    end

    @testset "3D conforming interface" begin
        left = hex_mesh(0.0, 1.0, tag="left")
        right = hex_mesh(1.0, 2.0, tag="right")

        add_mesh(left, right)

        @test length(left.nodes) == 12
        @test length(left.elems) == 2
        @test sort(getfield.(left.elems, :tag)) == ["left", "right"]
    end

    @testset "Invalid joins" begin
        @test_throws Exception add_mesh(quad_mesh(0.0, 1.0, tag="a"), quad_mesh(0.0, 1.0, tag="b"))
        @test_throws Exception add_mesh(quad_mesh(0.0, 1.0, tag="coarse"), split_quad_mesh(1.0, 2.0, tag="fine"))
    end

end
