using Serendip
using Test

printstyled("\nMesh generation on solids\n", color=:blue, bold=true)

nodes_count = [ 221, 441, 121, 341, 441 ]
shapes = [TRI3, TRI6, QUAD4, QUAD8, QUAD9]

for (n, shape) in zip(nodes_count, shapes)
    println("Generating mesh using ", shape.name)
    geo = GeoModel()
    add_block(geo, [0,0], [1,1], nx=10, ny=10, shape=shape)
    mesh = Mesh(geo)
    TR = @test length(mesh.nodes) == n
    println(TR)
end

nodes_count = [ 1331, 4961, 9261, 1331, 9261 ]
shapes = [HEX8, HEX20, HEX27, TET4, TET10]

for (n, shape) in zip(nodes_count, shapes)
    println("Generating mesh using ", shape.name)
    geo = GeoModel()
    add_block(geo, [0,0,0], [1,1,1], nx=10, ny=10, nz=10, shape=shape)
    mesh = Mesh(geo)
    TR = @test length(mesh.nodes) == n
    println(TR)
end



printstyled("\nMesh generation on trusses\n", color=:blue, bold=true)

println("\n2D")
coord = [ 0 0; 9 0; 18 0; 0 9; 9 9; 18 9.]
conn  = [ [1, 2], [1, 5], [2, 3], [2, 6], [2, 5], [2, 4], [3, 6], [3, 5], [4, 5], [5, 6] ]
mesh = Mesh(coord, conn)
TR = @test length(mesh.nodes) == 6
println(TR)

println("\n3D")
coord = [ 0.0 0.0 0.0; 0.0 1.0 0.0; 0.0 1.0 1.0]
conn  = [ [1, 3], [1, 2], [2, 3] ]
mesh = Mesh(coord, conn)
TR = @test length(mesh.nodes) == 3
println(TR)



printstyled("\nMesh generation with insets\n", color=:blue, bold=true)

shapes = [ HEX8, HEX20 ]

for shape in shapes
    println("\nGenerating mesh using ", shape.name)
    geo = GeoModel()
    add_block(geo, [0, 0, 0], [1, 1, 1], nx=8, ny=8, nz=8, shape=shape)

    p1 = add_point(geo, [0, 0, 0])
    p2 = add_point(geo, [1, 1, 1])
    l1 = add_line(geo, p1, p2)
    add_path(geo, [l1])

    mesh = Mesh(geo)
    TR = @test length(select(mesh.elems, :line)) == 8
    println(TR)
end



printstyled("\nMesh generation of joint cells\n", color=:blue, bold=true)

nodes_count = [ 108, 108, 36, 72, 648, 1620, 216, 540]
elems_count = [ 84, 39, 21, 21, 432, 432, 81, 81 ]
shapes = [ TRI3, TRI6, QUAD4, QUAD8, TET4, TET10, HEX8, HEX20 ]

for (n, e, shape) in zip(nodes_count, elems_count, shapes)
# for shape in shapes
    println("\nGenerating mesh using ", shape.name)
    geo = GeoModel()
    if shape.ndim==2
        add_block(geo, [0,0], [1,1], nx=3, ny=3, shape=shape)
    else
        add_block(geo, [0,0,0], [1,1,1], nx=3, ny=3, nz=3, shape=shape)
    end
    mesh = Mesh(geo)
    add_cohesive_elements(mesh)

    TR = @test length(mesh.nodes) == n && length(mesh.elems) == e
    println(TR)
end


printstyled("\nMesh generation of joint cells with insets\n", color=:blue, bold=true)

nodes_count = [ 121, 43, 661, 223 ]
elems_count = [ 96, 27, 444, 87 ]
shapes = [ TRI3, QUAD4, TET4, HEX8 ]

for (n, e, shape) in zip(nodes_count, elems_count, shapes)
    println("\nGenerating mesh using ", shape.name)
    geo = GeoModel()
    if shape.ndim==2
        add_block(geo, [0,0], [1,1], nx=3, ny=3, shape=shape)
        p1 = add_point(geo, [0, 0, 0])
        p2 = add_point(geo, [1, 1, 0])
    else
        add_block(geo, [0,0,0], [1,1,1], nx=3, ny=3, nz=3, shape=shape)
        p1 = add_point(geo, [0, 0, 0])
        p2 = add_point(geo, [1, 1, 1])
    end
    l1 = add_line(geo, p1, p2)
    add_path(geo, [l1])

    mesh = Mesh(geo)
    add_cohesive_elements(mesh)

    # @show length(mesh.nodes), length(mesh.elems)

    TR = @test length(mesh.nodes) == n && length(mesh.elems) == e
    println(TR)
end


# println("\nGenerating mesh using QUAD8 and 3 layers")
# bl   = Block( [0 0; 1 1], nx=4, ny=4, shape=QUAD8)
# bli  = BlockInset( [ 0 0; 1 1] )
# mesh = Mesh(bl, bli)
# mesh = insert_cohesive_elements!(mesh, layers=3)
# @test length(mesh.nodes) == 158
# TR = @test length(mesh.elems) == 48
# println(TR)
