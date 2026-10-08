using Serendip
using Test


@announced_testset "Cell measures and diameters" begin
    line = Cell(
        Serendip.LIN2,
        :line,
        [Node(0.0, 0.0), Node(3.0, 4.0)],
    )
    @test cell_measure(line) ≈ 5.0
    @test cell_diameter(line) ≈ 5.0

    triangle = Cell(
        Serendip.TRI3,
        :surface,
        [Node(0.0, 0.0), Node(4.0, 0.0), Node(0.0, 3.0)],
    )
    @test cell_measure(triangle) ≈ 6.0
    @test cell_diameter(triangle) ≈ 5.0

    quadrilateral = Cell(
        Serendip.QUAD4,
        :surface,
        [
            Node(0.0, 0.0),
            Node(3.0, 0.0),
            Node(3.0, 4.0),
            Node(0.0, 4.0),
        ],
    )
    @test cell_measure(quadrilateral) ≈ 12.0
    @test cell_diameter(quadrilateral) ≈ 5.0
end
