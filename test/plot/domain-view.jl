using Serendip
using Test
using LinearAlgebra

geo = GeoModel(quiet=true)
add_block(geo, [0, 0, 0], 1, 1, 1, nx=1, ny=1, nz=1, tag="solids")
mesh = Mesh(geo, quiet=true)

mapper = RegionMapper()
add_mapping(mapper, "solids", MechSolid, LinearElastic, E=1e4, nu=0.25)
model = FEModel(mesh, mapper, quiet=true)

geo2d = GeoModel(quiet=true)
add_block(geo2d, [0, 0, 0], 1, 1, 0, nx=1, ny=1, shape=:quad4, tag="solids")
mesh2d = Mesh(geo2d, quiet=true)
mesh.node_fields["flow"] = hcat([node.coord[1] + 1.0 for node in mesh.nodes], [0.5 * node.coord[2] + 0.25 for node in mesh.nodes], ones(length(mesh.nodes)))
mesh.node_fields["zero_flow"] = zeros(length(mesh.nodes), 3)
mesh.node_fields["short_flow"] = ones(length(mesh.nodes), 2)
mesh.node_fields["bad_flow"] = collect(1.0:length(mesh.nodes))
mesh2d.node_fields["flow"] = hcat([node.coord[1] + 1.0 for node in mesh2d.nodes], [0.5 * node.coord[2] + 0.25 for node in mesh2d.nodes])
mesh2d.node_fields["edge_flow"] = hcat(fill(0.5, length(mesh2d.nodes)), zeros(length(mesh2d.nodes)))
mesh2d.node_fields["tiny_flow"] = hcat(fill(1e-4, length(mesh2d.nodes)), zeros(length(mesh2d.nodes)))
mesh2d.node_fields["zero_flow"] = zeros(length(mesh2d.nodes), 2)
mesh2d.node_fields["short_flow"] = ones(length(mesh2d.nodes), 1)
mesh2d.node_fields["bad_flow"] = collect(1.0:length(mesh2d.nodes))
mesh2d.node_fields["temp"] = collect(1.0:length(mesh2d.nodes))

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

function layer_plot(mesh; kwargs...)
    plot = DomainPlot(quiet=true)
    add_plot(plot, mesh; kwargs...)
    return plot.layers[1]
end

function two_hex_interface_mesh()
    geo = GeoModel(quiet=true)
    add_block(geo, [0, 0, 0], 1, 1, 1, nx=2, ny=1, nz=1, tag="solids")
    mesh = Mesh(geo, quiet=true)

    for cell in select(mesh.elems, :solid)
        xmid = sum(node.coord.x for node in cell.nodes) / length(cell.nodes)
        cell.tag = xmid < 0.5 ? "left" : "right"
    end

    return mesh
end

function configured_outline_geom_keys(mesh; warp=0.0)
    plot = DomainPlot(quiet=true)
    add_plot(plot, mesh; warp=warp)
    Serendip.configure!(plot)

    keys = Tuple[]
    for edge in values(plot.layers[1].outline_edges_d)
        push!(keys, Tuple(sort!(collect(Serendip.node_pos_key(node) for node in edge.nodes))))
    end
    return sort!(keys)
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
@test DomainPlot(mesh2d, vector_field="flow").vector_field == "flow"
@test DomainPlot(mesh2d, vector_field="flow").arrow_length == 12.0
@test DomainPlot(mesh2d, vector_field="flow").arrow_width == 0.5
@test DomainPlot(mesh2d, vector_field="flow").arrow_color == Color(:black)
@test view_plot(model).layers[1].field_kind == :auto
@test DomainPlot(up=:x, quiet=true).up == :x
@test DomainPlot(up=:y, quiet=true).up == :y
@test DomainPlot(up=:z, quiet=true).up == :z
@test_throws ArgumentError DomainPlot(up=:w, quiet=true)
@testset "Line element styling" begin
    default_line_elem_layer = layer_plot(mesh2d)
    @test default_line_elem_layer.line_elem_color == :auto
    @test default_line_elem_layer.line_elem_width == 0.6
    @test Serendip._domain_default_line_elem_rgb(default_line_elem_layer) == (0.8, 0.2, 0.1)

    wireframe_outline_layer = layer_plot(mesh2d, view_mode=:wireframe, show_outline=true)
    @test wireframe_outline_layer.view_mode == :wireframe
    @test wireframe_outline_layer.show_outline
    @test wireframe_outline_layer.line_elem_width == 0.6
    @test Serendip._domain_default_line_elem_rgb(wireframe_outline_layer) == (0.8, 0.2, 0.1)

    override_line_elem_layer = layer_plot(mesh2d, line_elem_color=:green)
    @test override_line_elem_layer.line_elem_color == Color(:green)

    custom_style_layer = layer_plot(
        mesh2d,
        face_color=Color(:steelblue),
        line_color=(0.25, 0.25, 0.25),
        outline_color=Color(:black),
        line_elem_color=(0.1, 0.6, 0.2),
    )
    @test custom_style_layer.face_color == Color(:steelblue)
    @test custom_style_layer.line_color == Color((0.25, 0.25, 0.25))
    @test custom_style_layer.outline_color == Color(:black)
    @test custom_style_layer.line_elem_color == Color((0.1, 0.6, 0.2))

    outline_default_layer = layer_plot(mesh2d, view_mode=:outline, outline_width=0.8, outline_color=Color(:black))
    @test outline_default_layer.line_elem_width == 0.8
    @test outline_default_layer.line_elem_color == :auto
    @test Serendip._domain_default_line_elem_rgb(outline_default_layer) == (0.0, 0.0, 0.0)
end

@testset "Vector overlay configuration" begin
    vector_layer = layer_plot(mesh2d, vector_field="flow", arrow_length=9.0, arrow_width=0.8, arrow_color=:red)
    @test vector_layer.vector_field == "flow"
    @test vector_layer.arrow_length == 9.0
    @test vector_layer.arrow_width == 0.8
    @test vector_layer.arrow_color == Color(:red)

    scalar_vector_plot = DomainPlot(mesh2d, field="temp", field_kind=:node, vector_field="flow", arrow_color=:green)
    Serendip.configure!(scalar_vector_plot)
    scalar_vector_layer = scalar_vector_plot.layers[1]
    @test scalar_vector_layer.field == "temp"
    @test scalar_vector_layer.vector_field == "flow"
    @test length(scalar_vector_layer.values) == length(mesh2d.nodes)
    @test size(scalar_vector_layer.vector_values) == (length(mesh2d.nodes), 2)
    @test maximum(norm.(eachrow(scalar_vector_layer.vector_values))) > 0.0

    zero_vector_plot = DomainPlot(mesh2d, vector_field="zero_flow")
    Serendip.configure!(zero_vector_plot)
    @test all(iszero, zero_vector_plot.layers[1].vector_values)

    save(zero_vector_plot, "output/domain-vectors-zero.pdf")
    @test isfile("output/domain-vectors-zero.pdf")

    plain_plot = DomainPlot(mesh2d)
    Serendip.configure!(plain_plot)
    endpoint_plot = DomainPlot(mesh2d, vector_field="tiny_flow")
    Serendip.configure!(endpoint_plot)
    @test endpoint_plot.canvas.limits[3] - plain_plot.canvas.limits[3] > 0.01

    @test_throws ArgumentError layer_plot(mesh2d, arrow_length=0.0)
    @test_throws ArgumentError layer_plot(mesh2d, arrow_width=0.0)
    @test_throws ArgumentError layer_plot(mesh2d, arrow_color=:notacolor)
    @test_throws ErrorException begin
        plot = DomainPlot(mesh2d, vector_field="missing_flow")
        Serendip.configure!(plot)
    end
    @test_throws ErrorException begin
        plot = DomainPlot(mesh2d, vector_field="bad_flow")
        Serendip.configure!(plot)
    end
    @test_throws ErrorException begin
        plot = DomainPlot(mesh2d, vector_field="short_flow")
        Serendip.configure!(plot)
    end
    @test_throws ErrorException begin
        plot = DomainPlot(mesh, vector_field="short_flow")
        Serendip.configure!(plot)
    end
end

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
Serendip._domain_correct_shared_edge_priority_3d!(sorted_render_elems, render_tol)
@test surface3d_plot.render_elems == sorted_render_elems

cube_vector_plot = DomainPlot(mesh, vector_field="flow")
Serendip.configure!(cube_vector_plot)
@test size(cube_vector_plot.layers[1].vector_values) == (length(mesh.nodes), 2)
@test maximum(norm.(eachrow(cube_vector_plot.layers[1].vector_values))) > 0.0
@test length(Serendip._domain_overlay_nodes(cube_vector_plot)[1]) < length(cube_vector_plot.layers[1].nodes)

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

shared_edge_chain = sort([shared_edge_line, shared_edge_surface], by=Serendip._domain_render_depth_key)
Serendip._domain_correct_shared_edge_priority_3d!(shared_edge_chain, depth_tol)
@test shared_edge_chain[end] === shared_edge_line

many_surfaces = [
    make_render_elem(render_layer, :surface, [0.0, 0.0, 0.0, 0.0], index_in_layer=200+i)
    for i in 1:10
]
late_line = make_render_elem(render_layer, :line, [0.0, 0.0], index_in_layer=199)
correction_chain = sort([late_line; many_surfaces], by=Serendip._domain_render_depth_key)
@test correction_chain[1] === late_line

@testset "Outline smoke test across cohesive/contact interfaces" begin
    base_mesh = two_hex_interface_mesh()
    base_outline = configured_outline_geom_keys(base_mesh)
    @test length(base_outline) == 16

    cohesive_mesh = two_hex_interface_mesh()
    add_cohesive_elements(cohesive_mesh, "left", tag="cohesive", implicit=false, quiet=true)
    cohesive_outline = configured_outline_geom_keys(cohesive_mesh)
    @test cohesive_outline == base_outline

    contact_mesh = two_hex_interface_mesh()
    add_contact_elements(contact_mesh, "left", "right", tag="contact", quiet=true)
    contact_outline = configured_outline_geom_keys(contact_mesh)
    @test contact_outline == base_outline
end

@testset "Tolerant geometric matching" begin
    quad_a = Cell(get_shape(:quad4), :solid, Node[
        Node(0.0, 0.0, 0.0, id=1),
        Node(1.0, 0.0, 0.0, id=2),
        Node(1.0, 1.0, 0.0, id=3),
        Node(0.0, 1.0, 0.0, id=4),
    ])
    quad_b = Cell(get_shape(:quad4), :solid, Node[
        Node(5e-6, 0.0, 0.0, id=5),
        Node(1.0 + 5e-6, 0.0, 0.0, id=6),
        Node(1.0 + 5e-6, 1.0, 0.0, id=7),
        Node(5e-6, 1.0, 0.0, id=8),
    ])

    @test length(get_outline_edges([quad_a, quad_b])) == 8
    @test isempty(get_outline_edges([quad_a, quad_b], tol=1e-5))
end

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

save(view_plot(mesh2d, title="2D", axes=:bottom_left), "output/domainplot-up-axes-2d.pdf")
save(view_plot(model, title="up = :z", axes=:top_right, up=:z), "output/domainplot-up-axes-z.pdf")
save(view_plot(model, title="up = :x", axes=:top_right, up=:x), "output/domainplot-up-axes-x.pdf")
save(view_plot(model, title="up = :y", axes=:top_right, up=:y), "output/domainplot-up-axes-y.pdf")
