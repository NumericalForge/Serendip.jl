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

function make_render_elem(layer, role::Symbol, zcoords::Vector{Float64}; layer_index=1, index_in_layer=1, xcoords=nothing, ycoords=nothing, node_ids=nothing)
    xcoords === nothing && (xcoords = role == :line ? [0.0, 1.0] : [0.0, 1.0, 1.0, 0.0])
    ycoords === nothing && (ycoords = role == :line ? [0.0, 0.0] : [0.0, 0.0, 1.0, 1.0])
    node_ids === nothing && (node_ids = role == :line ? [10*index_in_layer + 1, 10*index_in_layer + 2] : [10*index_in_layer + 1, 10*index_in_layer + 2, 10*index_in_layer + 3, 10*index_in_layer + 4])
    nodes = if role == :line
        Node[
            Node(xcoords[1], ycoords[1], zcoords[1], id=node_ids[1]),
            Node(xcoords[2], ycoords[2], zcoords[2], id=node_ids[2]),
        ]
    else
        Node[
            Node(xcoords[1], ycoords[1], zcoords[1], id=node_ids[1]),
            Node(xcoords[2], ycoords[2], zcoords[2], id=node_ids[2]),
            Node(xcoords[3], ycoords[3], zcoords[3], id=node_ids[3]),
            Node(xcoords[4], ycoords[4], zcoords[4], id=node_ids[4]),
        ]
    end
    shape = role == :line ? get_shape(:lin2) : get_shape(:quad4)
    elem = Cell(shape, role, nodes, id=index_in_layer)
    xmin = minimum(xcoords)
    xmax = maximum(xcoords)
    ymin = minimum(ycoords)
    ymax = maximum(ycoords)
    zmin = minimum(zcoords)
    zmax = maximum(zcoords)
    zavg = sum(zcoords) / length(zcoords)
    fallback_depth = 0.9 * zavg + 0.1 * zmin
    return Serendip.DomainRenderElem(layer, elem, 1.0, layer_index, index_in_layer, xmin, xmax, ymin, ymax, zmin, zmax, zavg, fallback_depth)
end

@test DomainPlot(quiet=true).up == :z
@test view_plot(model).layers[1].field_kind == :auto
@test DomainPlot(up=:x, quiet=true).up == :x
@test DomainPlot(up=:y, quiet=true).up == :y
@test DomainPlot(up=:z, quiet=true).up == :z
@test_throws ArgumentError DomainPlot(up=:w, quiet=true)
@test view_plot(mesh2d, title="line elem color").layers[1].line_elem_color == :auto

wireframe_outline_plot = DomainPlot(quiet=true)
add_plot(wireframe_outline_plot, mesh2d, view_mode=:wireframe, show_outline=true)
@test wireframe_outline_plot.layers[1].view_mode == :wireframe
@test wireframe_outline_plot.layers[1].show_outline

line_elem_plot = DomainPlot(quiet=true)
add_plot(line_elem_plot, mesh2d, line_elem_color=:green)
@test line_elem_plot.layers[1].line_elem_color == Color(:green)

custom_color_plot = DomainPlot(quiet=true)
add_plot(
    custom_color_plot,
    mesh2d,
    face_color=Color(:steelblue),
    line_color=(0.25, 0.25, 0.25),
    outline_color=Color(:black),
    line_elem_color=(0.1, 0.6, 0.2),
)
@test custom_color_plot.layers[1].face_color == Color(:steelblue)
@test custom_color_plot.layers[1].line_color == Color((0.25, 0.25, 0.25))
@test custom_color_plot.layers[1].outline_color == Color(:black)
@test custom_color_plot.layers[1].line_elem_color == Color((0.1, 0.6, 0.2))

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
render_tol = Serendip._domain_render_depth_tolerance(surface3d_plot)
sorted_render_elems = sort(copy(surface3d_plot.render_elems), by=Serendip._domain_render_depth_key)
Serendip._domain_correct_render_order_3d!(sorted_render_elems, render_tol)
Serendip._domain_correct_shared_edge_priority_3d!(sorted_render_elems, render_tol)
@test surface3d_plot.render_elems == sorted_render_elems

render_layer = surface3d_plot.layers[1]
depth_tol = 1e-6

back_surface = make_render_elem(render_layer, :surface, [0.4, 0.4, 0.4, 0.4], index_in_layer=1)
front_surface = make_render_elem(render_layer, :surface, [-1.0, -1.0, -1.0, -1.0], index_in_layer=2)
@test Serendip._domain_render_order_3d(back_surface, front_surface, depth_tol) < 0

front_line = make_render_elem(render_layer, :line, [0.0, 0.0], index_in_layer=3)
coplanar_surface = make_render_elem(render_layer, :surface, [0.0, 0.0, 0.0, 0.0], index_in_layer=4)
@test Serendip._domain_render_order_3d(coplanar_surface, front_line, depth_tol) < 0

coplanar_solid_face = make_render_elem(render_layer, :solid, [0.0, 0.0, 0.0, 0.0], index_in_layer=41)
@test Serendip._domain_render_order_3d(coplanar_solid_face, front_line, depth_tol) < 0

front_edge_line = make_render_elem(render_layer, :line, [0.0, 1.0], index_in_layer=42)
shared_front_surface = make_render_elem(render_layer, :surface, [0.2, 0.5, 0.6, 1.0], index_in_layer=43)
@test Serendip._domain_render_order_3d(shared_front_surface, front_edge_line, depth_tol) < 0

shared_edge_line = make_render_elem(render_layer, :line, [0.3, 0.3], index_in_layer=44, xcoords=[0.0, 1.0], ycoords=[0.0, 0.0], node_ids=[501, 502])
shared_edge_surface = make_render_elem(render_layer, :surface, [0.1, 0.4, 0.8, 0.9], index_in_layer=45, xcoords=[0.0, 1.0, 1.0, 0.0], ycoords=[0.0, 0.0, 1.0, 1.0], node_ids=[501, 502, 503, 504])
@test Serendip._domain_shared_node_count(shared_edge_line.elem, shared_edge_surface.elem) == 2
@test Serendip._domain_render_shared_edge_order_3d(shared_edge_surface, shared_edge_line, depth_tol) < 0
@test Serendip._domain_render_shared_edge_order_3d(shared_edge_line, shared_edge_surface, depth_tol) > 0

hidden_line = make_render_elem(render_layer, :line, [0.8, 0.8], index_in_layer=5)
front_cover_surface = make_render_elem(render_layer, :surface, [0.1, 0.1, 0.1, 0.1], index_in_layer=6)
@test Serendip._domain_render_order_3d(hidden_line, front_cover_surface, depth_tol) < 0

ambiguous_far_surface = make_render_elem(render_layer, :surface, [-0.2, 0.4, 0.1, 0.0], index_in_layer=7)
ambiguous_near_surface = make_render_elem(render_layer, :surface, [-0.1, 0.5, 0.2, 0.2], index_in_layer=8)
@test Serendip._domain_render_order_3d(ambiguous_far_surface, ambiguous_near_surface, depth_tol) > 0

shared_interval_back_surface = make_render_elem(render_layer, :surface, [0.1, 0.5, 0.6, 0.8], index_in_layer=9)
shared_interval_front_surface = make_render_elem(render_layer, :surface, [-0.1, 0.4, 0.7, 0.9], index_in_layer=10)
@test Serendip._domain_render_order_3d(shared_interval_back_surface, shared_interval_front_surface, depth_tol) < 0

disjoint_front_line = make_render_elem(render_layer, :line, [-1.0, -1.0], index_in_layer=11, xcoords=[3.0, 4.0], ycoords=[3.0, 3.0])
disjoint_back_surface = make_render_elem(render_layer, :surface, [0.0, 0.0, 0.0, 0.0], index_in_layer=12)
@test !Serendip._domain_render_overlap_2d(disjoint_front_line, disjoint_back_surface, depth_tol)
@test Serendip._domain_render_order_3d(disjoint_front_line, disjoint_back_surface, depth_tol) < 0

baseline_pair = sort([ambiguous_far_surface, ambiguous_near_surface], by=Serendip._domain_render_depth_key)
Serendip._domain_correct_render_order_3d!(baseline_pair, depth_tol)
@test baseline_pair == [ambiguous_near_surface, ambiguous_far_surface]

shared_edge_chain = sort([shared_edge_line, shared_edge_surface], by=Serendip._domain_render_depth_key)
Serendip._domain_correct_render_order_3d!(shared_edge_chain, depth_tol)
Serendip._domain_correct_shared_edge_priority_3d!(shared_edge_chain, depth_tol)
@test shared_edge_chain[end] === shared_edge_line

many_surfaces = [
    make_render_elem(render_layer, :surface, [0.0, 0.0, 0.0, 0.0], index_in_layer=200+i)
    for i in 1:10
]
late_line = make_render_elem(render_layer, :line, [0.0, 0.0], index_in_layer=199)
correction_chain = sort([late_line; many_surfaces], by=Serendip._domain_render_depth_key)
@test correction_chain[1] === late_line
Serendip._domain_correct_render_order_3d!(correction_chain, depth_tol)
@test correction_chain[end] === late_line

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
