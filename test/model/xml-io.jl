using Serendip
using Test


@announced_testset "FEModel XML generation" begin
    geometry = GeoModel()
    add_block(
        geometry,
        [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0],
        nx=1,
        ny=1,
        shape=:quad4,
        tag="solids",
    )
    mesh = Mesh(geometry, quiet=true)

    mapper = RegionMapper()
    add_mapping(mapper, "solids", MechSolid, LinearElastic, E=100.0, nu=0.2)
    model = FEModel(mesh, mapper, quiet=true)

    mktempdir() do dir
        filename = joinpath(dir, "model.xml")
        save(model, filename, quiet=true)
        @test isfile(filename)

        document = XmlDocument(filename)
        @test document.root.name == "FEModel"
        @test document.root("Materials") isa XmlElement
        @test document.root("Nodes") isa XmlElement
        @test document.root("Elements").children[1].name == "MechSolid"
    end
end
