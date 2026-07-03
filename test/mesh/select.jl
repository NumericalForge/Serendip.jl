using Serendip
using Test

function make_tagged_mesh()
    geo = GeoModel(quiet=true)
    add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0, nx=2, ny=1, shape=:quad4, tag="solids")
    mesh = Mesh(geo, quiet=true)

    for cell in select(mesh, :element, :solid)
        xmid = sum(node.coord.x for node in cell.nodes) / length(cell.nodes)
        cell.tag = xmid < 0.5 ? "left" : "right"
    end

    return mesh
end

mesh = make_tagged_mesh()
left = select(mesh, :element, "left")
right = select(mesh, :element, "right")

@test length(left) == 1
@test length(right) == 1
@test [cell.id for cell in select(mesh, :element, any_of("left", "right"))] == [cell.id for cell in select(mesh, :element, :solid)]
@test [cell.id for cell in select(mesh, :element, :solid, none_of("left"))] == [cell.id for cell in right]
@test [cell.id for cell in select(mesh, :element, any_of("left", x >= 0.5))] == [cell.id for cell in select(mesh, :element, :solid)]

redirected_nodes = select(mesh, :element, any_of("left", "right"), :node, x == 0.0)
@test !isempty(redirected_nodes)
@test all(node.coord.x == 0.0 for node in redirected_nodes)

@test_throws ArgumentError select(mesh, :element, any_of(:node, "left"))
@test_throws ArgumentError select(mesh, :element, none_of(:ip, "left"))

quad_nodes = Node[
    Node(0.0, 0.0, 0.0, id=1),
    Node(1.0, 0.0, 0.0, id=2),
    Node(1.0, 1.0, 0.0, id=3),
    Node(0.0, 1.0, 0.0, id=4),
]
solid_cell = Cell(get_shape(:quad4), :solid, quad_nodes; tag="a", id=1)
line_cell = Cell(get_shape(:lin2), :line, quad_nodes[1:2]; tag="a", id=2)
mixed = Cell[solid_cell, line_cell]

@test isempty(select(mixed, :solid, none_of("a")))
@test [cell.id for cell in select(mixed, :solid, "a"; invert=true)] == [2]
