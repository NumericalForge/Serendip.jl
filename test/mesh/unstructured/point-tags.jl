using Serendip, Test


function tagged_adjacent_occ_region(geo::GeoModel, tag::String)
    points = [
        add_point(geo, [1.0, 0.0, 0.0], tag=tag),
        add_point(geo, [2.0, 0.0, 0.0]),
        add_point(geo, [2.0, 1.0, 0.0]),
        add_point(geo, [1.0, 1.0, 0.0]),
    ]
    return add_polygon(geo, points, tag="occ")
end


@announced_testset "Geometric point tags" begin
    @announced_testset "Boundary and embedded point tags" begin
        geo = GeoModel(quiet=true)
        points = [
            add_point(geo, [0.0, 0.0, 0.0], tag="origin"),
            add_point(geo, [1.0, 0.0, 0.0]),
            add_point(geo, [1.0, 1.0, 0.0]),
            add_point(geo, [0.0, 1.0, 0.0]),
        ]
        add_polygon(geo, points, tag="body")
        add_point(geo, [0.5, 0.5, 0.0], embedded=true, tag="sensor")

        mesh = Mesh(geo, quiet=true)

        @test length(select(mesh, :node, "origin")) == 1
        @test length(select(mesh, :node, "sensor")) == 1
        @test !isempty(select(mesh.elems, "body"))
    end

    @announced_testset "Orphan tagged point does not create a node" begin
        geo = GeoModel(quiet=true)
        add_rectangle(geo, [0.0, 0.0, 0.0], 1.0, 1.0, tag="body")
        add_point(geo, [2.0, 2.0, 0.0], tag="orphan")

        mesh = Mesh(geo, quiet=true)

        @test isempty(select(mesh, :node, "orphan"))
        @test !any(node -> node.coord.x == 2.0 && node.coord.y == 2.0, mesh.nodes)
    end

    @announced_testset "Tag survives welding with a structured block" begin
        geo = GeoModel(quiet=true)
        add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0,
            nx=1, ny=1, shape=:quad4, tag="block")
        tagged_adjacent_occ_region(geo, "block-occ-interface")

        mesh = Mesh(geo, quiet=true)

        tagged = select(mesh, :node, "block-occ-interface")
        @test length(tagged) == 1
        @test tagged[1].coord.x == 1.0
        @test tagged[1].coord.y == 0.0
    end

    @announced_testset "Tag survives welding with a stored mesh" begin
        stored = Mesh(
            [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0],
            [[1, 2, 3, 4]],
            [:quad4],
            tag="stored",
            quiet=true,
        )
        geo = GeoModel(quiet=true)
        add_mesh(geo, stored)
        tagged_adjacent_occ_region(geo, "stored-occ-interface")

        mesh = Mesh(geo, quiet=true)

        tagged = select(mesh, :node, "stored-occ-interface")
        @test length(tagged) == 1
        @test tagged[1].coord.x == 1.0
        @test tagged[1].coord.y == 0.0
        @test all(isempty(node.tag) for node in stored.nodes)
    end
end
