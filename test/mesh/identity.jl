using Serendip
using Test


@announced_testset "Mesh identity and topology keys" begin
    n1 = Node(0.0, 0.0)
    n2 = Node(1.0, 0.0)
    n1_equivalent = Node(0.0, 0.0)

    @test n1 !== n1_equivalent
    @test !isequal(n1, n1_equivalent)

    cell1 = Cell(Serendip.LIN2, :line, [n1, n2])
    cell2 = Cell(Serendip.LIN2, :line, [n2, n1])
    @test !isequal(cell1, cell2)
    @test length(Set([cell1, cell2])) == 2

    repeated_key = Serendip._topology_key([n1, n2, n1])
    @test repeated_key == Serendip._topology_key([n1, n1, n2])
    @test repeated_key != Serendip._topology_key([n1, n2])

    n1.id = 1
    n2.id = 2
    numbered_copy = [Node(1.0, 0.0, id=2), Node(0.0, 0.0, id=1)]
    @test Serendip._topology_key([n1, n2]) == Serendip._topology_key(numbered_copy)

    geometric_copy = Cell(Serendip.LIN2, :line, [Node(1.0, 0.0), Node(0.0, 0.0)])
    @test Serendip._cell_key(cell1) == Serendip._cell_key(geometric_copy)

    nodes, patches = Serendip.get_patches([cell1, cell2])
    @test length(nodes) == 2
    @test sort(length.(patches)) == [2, 2]
end


@announced_testset "Explicit collection operations" begin
    block = Block(
        [0.0, 0.0],
        1.0,
        1.0,
        0.0;
        nx=1,
        ny=1,
        shape=Serendip.QUAD4,
    )
    blocks = [block]

    shallow = copy(blocks)
    copied = copy.(blocks)
    @test shallow !== blocks
    @test shallow[1] === block
    @test copied[1] !== block
    @test Serendip.get_coords(copied[1].points) == Serendip.get_coords(block.points)

    elements = Element{MechSolid}[]
    @test sort!(elements; by=element -> element.id) === elements
end


@announced_testset "Explicit edge and slice keys" begin
    geometry = GeoModel(quiet=true)
    add_block(
        geometry,
        [0.0, 0.0, 0.0],
        1.0,
        1.0,
        1.0;
        nx=1,
        ny=1,
        nz=1,
        shape=:hex8,
    )
    mesh = Mesh(geometry, quiet=true)

    surface = get_outer_facets(mesh.elems)
    @test length(Serendip.get_edges(surface)) == 12

    cut = slice(mesh, base=[0.5, 0.0, 0.0], axis=[1.0, 0.0, 0.0])
    @test length(cut.nodes) == 4
    @test length(cut.elems) == 1
end
