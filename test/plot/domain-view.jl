using Serendip
using Test

geo = GeoModel(quiet=true)
add_block(geo, [0, 0, 0], 1, 1, 1, nx=1, ny=1, nz=1, tag="solids")
mesh = Mesh(geo, quiet=true)

mapper = RegionMapper()
add_mapping(mapper, "solids", MechSolid, LinearElastic, E=1e4, nu=0.25)
model = FEModel(mesh, mapper, quiet=true)

geo2d = GeoModel(quiet=true)
add_block(geo2d, [0, 0, 0], 1, 1, 0, nx=1, ny=1, shape=:quad4, tag="solids")
mesh2d = Mesh(geo2d, quiet=true)

function projected_xy(coords, up)
    nodes = [Node(coord; id=i) for (i, coord) in enumerate(coords)]
    Serendip.project_to_2d!(nodes, 0.0, 0.0, 0.0, up)
    return [[node.coord[1], node.coord[2]] for node in sort(nodes, by=node -> node.id)]
end

function configured_coords(plot)
    Serendip.configure!(plot)
    return [[node.coord[1], node.coord[2], node.coord[3]] for node in sort(plot.nodes, by=node -> node.id)]
end

@test DomainPlot(model, quiet=true).up == :z
@test DomainPlot(model, quiet=true).field_kind == :auto
@test DomainPlot(model, up=:x, quiet=true).up == :x
@test DomainPlot(model, up=:y, quiet=true).up == :y
@test DomainPlot(model, up=:z, quiet=true).up == :z
@test_throws ArgumentError DomainPlot(model, up=:w, quiet=true)

default_coords = configured_coords(DomainPlot(model, quiet=true))
explicit_z_coords = configured_coords(DomainPlot(model, up=:z, quiet=true))
for (coord_default, coord_explicit) in zip(default_coords, explicit_z_coords)
    @test coord_default ≈ coord_explicit
end

zup_coords = projected_xy([[0.0, -1.0, -1.0], [0.0, 1.0, -1.0], [0.0, -1.0, 1.0], [0.0, 1.0, 1.0]], :z)
for (coord, expected) in zip(zup_coords, [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    @test coord ≈ expected
end

xup_coords = projected_xy([[-1.0, 0.0, -1.0], [1.0, 0.0, -1.0], [-1.0, 0.0, 1.0], [1.0, 0.0, 1.0]], :x)
for (coord, expected) in zip(xup_coords, [[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [1.0, 1.0]])
    @test coord ≈ expected
end

yup_coords = projected_xy([[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [-1.0, 1.0, 0.0], [1.0, 1.0, 0.0]], :y)
for (coord, expected) in zip(yup_coords, [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    @test coord ≈ expected
end

axes_grid = ChartGrid(
    title="DomainPlot Transformations",
    size=(16cm, 12cm),
)

add_chart(axes_grid, DomainPlot(mesh2d, title="2D", axes=:bottom_left, quiet=true), (1, 1))
add_chart(axes_grid, DomainPlot(model, title="up = :z", axes=:top_right, up=:z, quiet=true), (1, 2))
add_chart(axes_grid, DomainPlot(model, title="up = :x", axes=:top_right, up=:x, quiet=true), (2, 1))
add_chart(axes_grid, DomainPlot(model, title="up = :y", axes=:top_right, up=:y, quiet=true), (2, 2))

save(axes_grid, "output/domainplot-up-axes.pdf")
save(axes_grid, "output/domainplot-up-axes.png")

@test isfile("output/domainplot-up-axes.pdf")
@test isfile("output/domainplot-up-axes.png")
