using Serendip, Test

@announced_testset "GeoModel select" begin
    geo = GeoModel(quiet=true)
    base = add_rectangle(geo, [0.0, 0.0, 0.0], 2.0, 1.0, tag="body")

    edges = select(geo, :edge)
    @test length(edges) == 4

    bottom = select(geo, :edge, y == 0)
    @test length(bottom) == 1
    @test all(edge.tag == "" for edge in bottom)

    tagged_bottom = select(geo, :edge, (y == 0, x >= 0); tag="bottom")
    @test length(tagged_bottom) == 1
    @test tagged_bottom[1].tag == "bottom"
    @test select(geo, :edge, "bottom")[1] === tagged_bottom[1]
    @test [edge.id for edge in select(geo, :edge, "bottom", x <= 2)] == [edge.id for edge in tagged_bottom]

    left = select(geo, :edge, x == 0)
    @test length(left) == 1
    @test [edge.id for edge in select(geo, :edge, x == 0, [0.0, 0.5, 0.0])] == [edge.id for edge in left]
    @test length(select(geo, :edge, :all)) == length(edges)
    @test isempty(select(geo, :edge, :none))
    @test length(select(geo, :edge, y == 0, invert=true)) == 3

    top = select(geo, :surface, [1.0, 0.5, 0.0])
    @test length(top) == 1
    @test top[1] === base
end

@announced_testset "GeoModel select uses bbox corners conservatively" begin
    geo = GeoModel(quiet=true)
    p1 = add_point(geo, [0.0, 0.0, 0.0])
    p2 = add_point(geo, [1.0, 1.0, 0.0])
    add_line(geo, p1, p2, "diag")

    @test isempty(select(geo, :edge, x == y))
    @test length(select(geo, :edge, "diag")) == 1
end

@announced_testset "GeoModel select isolates transfinite vertical lines" begin
    geo = GeoModel(quiet=true)
    p1 = add_point(geo, [0, 0, -0.5])
    p2 = add_point(geo, [0, 0, 0])
    p3 = add_point(geo, [1, 0, -0.5])
    p4 = add_point(geo, [2, 0, 0])
    p5 = add_point(geo, [2, 0, -0.5])

    l1 = add_line(geo, p1, p2)
    a1 = add_circle_arc(geo, p2, p3, p4)
    l2 = add_line(geo, p4, p5)
    extrude(geo, [l1, a1, l2], [0, 3, 0], recombine=true)

    vlines = select(geo, :edge, or(y == 0, y == 3), z <= 0)
    @test length(vlines) == 4
end
