using Serendip
using Test
using LinearAlgebra: norm

@testset "Path Tips Are Global Endpoints" begin
    Xstart = [0.1, 0.2, 0.0]
    Xmid   = [0.5, 0.5, 0.0]
    Xend   = [0.9, 0.8, 0.0]
    atol   = 1e-6

    build_mesh(tips) = begin
        geo = GeoModel(quiet=true)
        add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0, nx=8, ny=8, shape=QUAD4)

        p1 = add_point(geo, Xstart)
        p2 = add_point(geo, Xmid)
        p3 = add_point(geo, Xend)
        l1 = add_line(geo, p1, p2)
        l2 = add_line(geo, p2, p3)

        add_path(geo, [l1, l2], tag="bars", interface_tag="joints", tips=tips, tip_tag="tips")
        return Mesh(geo)
    end

    tip_coords(mesh) = [e.nodes[end].coord for e in mesh.elems if e.role == :tip]
    has_tip_at(coords, X) = any(norm(c - X) <= atol for c in coords)

    mesh_none = build_mesh(:none)
    coords_none = tip_coords(mesh_none)
    @test length(coords_none) == 0

    mesh_start = build_mesh(:start)
    coords_start = tip_coords(mesh_start)
    @test length(coords_start) == 1
    @test has_tip_at(coords_start, Xstart)
    @test !has_tip_at(coords_start, Xend)

    mesh_end = build_mesh(:end)
    coords_end = tip_coords(mesh_end)
    @test length(coords_end) == 1
    @test !has_tip_at(coords_end, Xstart)
    @test has_tip_at(coords_end, Xend)

    mesh_both = build_mesh(:both)
    coords_both = tip_coords(mesh_both)
    @test length(coords_both) == 2
    @test has_tip_at(coords_both, Xstart)
    @test has_tip_at(coords_both, Xend)
end
