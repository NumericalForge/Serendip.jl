using Serendip
using Test

@announced_testset "Revolve line meshes" begin
    shapes = (:lin2, :lin3)
    data = ((48, 49), (48, 145))
    for (shape, expected) in zip(shapes, data)
        @announced_testset "$shape" begin
            geo = GeoModel()
            add_block(geo, [0, 0, 0], 1, 1, 0; n=4, shape=shape, tag="solids")
            mesh = Mesh(geo)
            mesh = revolve(mesh, base=[0, 0, 0], axis=[0, 1, 0], n=12)
            @test (length(mesh.elems), length(mesh.nodes)) == expected
        end
    end
end

@announced_testset "Revolve surface meshes" begin
    shapes = (:tri3, :tri6, :quad4, :quad8)
    data = ((768, 437), (384, 1113), (192, 245), (192, 921))
    for (shape, expected) in zip(shapes, data)
        @announced_testset "$shape" begin
            geo = GeoModel()
            add_block(geo, [0, 0, 0], 1, 1, 0; nx=4, ny=4, shape=shape, tag="solids")
            mesh = Mesh(geo)
            mesh = revolve(mesh, base=[0, 0, 0], axis=[0, 1, 0], n=12)
            @test (length(mesh.elems), length(mesh.nodes)) == expected
        end
    end
end

@announced_testset "Revolve node and cord" begin
    data = ((4, 9), (64, 193))

    mesh = revolve(Node([1, 0, 0]), minangle=0, maxangle=90, base=[0, 0, 0], axis=[0, 1, 0], n=4)
    @test (length(mesh.elems), length(mesh.nodes)) == data[1]

    mesh = revolve(mesh, base=[0, 0, 0], axis=[0, 0, 1], n=16)
    @test (length(mesh.elems), length(mesh.nodes)) == data[2]
end
