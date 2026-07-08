using Serendip, Test

mktempdir() do dir
    filename = joinpath(dir, "box.step")

    source = GeoModel(size=0.5, quiet=true)
    add_box(source, [0.0, 0.0, 0.0], 1.0, 1.0, 1.0; tag="box")
    save(source, filename, true)

    geo = GeoModel(filename; size=0.5, quiet=true)
    mesh = Mesh(geo; ndim=3, quiet=true)

    @test mesh.ctx.ndim == 3
    @test length(mesh.elems) > 0

    unsupported = joinpath(dir, "box.stl")
    touch(unsupported)
    @test_throws Exception GeoModel(unsupported; quiet=true)
end
