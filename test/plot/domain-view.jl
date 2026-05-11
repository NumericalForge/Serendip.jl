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

function view_plot(mesh; up=:z, title="", axes=:none)
    plot = DomainPlot(up=up, title=title, axes=axes, quiet=true)
    add_plot(plot, mesh)
    return plot
end

function surface_role_mesh(ndim)
    mesh = Mesh(ndim)
    mesh.nodes = Node[
        Node(0.0, 0.0, 0.0, id=1),
        Node(1.0, 0.0, 0.0, id=2),
        Node(1.0, 1.0, 0.0, id=3),
        Node(0.0, 1.0, 0.0, id=4),
    ]
    mesh.elems = Cell[Cell(get_shape(:quad4), :surface, mesh.nodes, id=1)]

    if ndim == 3
        append!(mesh.nodes, Node[
            Node(0.0, 0.0, 1.0, id=5),
            Node(1.0, 0.0, 1.0, id=6),
            Node(1.0, 1.0, 1.0, id=7),
            Node(0.0, 1.0, 1.0, id=8),
        ])
        push!(mesh.elems, Cell(get_shape(:quad4), :surface, mesh.nodes[5:8], id=2))
    end

    return mesh
end

@test DomainPlot(quiet=true).up == :z
@test view_plot(model).layers[1].field_kind == :auto
@test DomainPlot(up=:x, quiet=true).up == :x
@test DomainPlot(up=:y, quiet=true).up == :y
@test DomainPlot(up=:z, quiet=true).up == :z
@test_throws ArgumentError DomainPlot(up=:w, quiet=true)

surface2d_plot = view_plot(surface_role_mesh(2))
Serendip.configure!(surface2d_plot)
@test length(surface2d_plot.layers[1].elems) == 1
@test surface2d_plot.layers[1].elems[1].role == :surface
@test length(surface2d_plot.nodes) == 4

surface3d_plot = view_plot(surface_role_mesh(3))
Serendip.configure!(surface3d_plot)
@test length(surface3d_plot.layers[1].elems) == 2
@test all(elem.role == :surface for elem in surface3d_plot.layers[1].elems)
@test length(surface3d_plot.render_elems) == 2
render_keys = Serendip._domain_render_depth_key.(surface3d_plot.render_elems)
@test render_keys == sort(render_keys)

default_coords = configured_coords(view_plot(model))
explicit_z_coords = configured_coords(view_plot(model, up=:z))
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

add_chart(axes_grid, view_plot(mesh2d, title="2D", axes=:bottom_left), (1, 1))
add_chart(axes_grid, view_plot(model, title="up = :z", axes=:top_right, up=:z), (1, 2))
add_chart(axes_grid, view_plot(model, title="up = :x", axes=:top_right, up=:x), (2, 1))
add_chart(axes_grid, view_plot(model, title="up = :y", axes=:top_right, up=:y), (2, 2))

save(axes_grid, "output/domainplot-up-axes.pdf")
save(axes_grid, "output/domainplot-up-axes.png")

@test isfile("output/domainplot-up-axes.pdf")
@test isfile("output/domainplot-up-axes.png")
