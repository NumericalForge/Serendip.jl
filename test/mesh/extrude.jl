using Serendip
using Test

@announced_testset "Extrude mesh" begin
    geo = GeoModel()
    add_block(geo, [0, 0, 0], 1, 1, 0; nx=3, ny=3, shape=:quad4)
    mesh = Mesh(geo)
    mesh = extrude(mesh, axis=[0, 0, 1], length=4, n=10, quiet=true)
    @test length(mesh.elems) == 90
end

@announced_testset "Extrude normal to mesh" begin
    geo = GeoModel()
    add_block(geo, [0, 0, 0], 1, 1, 0; n=6, shape=:lin3)
    mesh = Mesh(geo)
    mesh = revolve(mesh, minangle=45, maxangle=135, n=6, axis=[0, -1, 0], base=[0, 0, 0])
    mesh = extrude(mesh, length=4.25, n=6, quiet=true)
    @test length(mesh.elems) == 216
end
