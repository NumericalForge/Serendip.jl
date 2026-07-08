using Serendip
using Test

@announced_testset "Solid blocks in 2d" begin
    nodes_count = [221, 441, 121, 341, 441]
    shapes = [:tri3, :tri6, :quad4, :quad8, :quad9]

    for (n, shape) in zip(nodes_count, shapes)
        @announced_testset "$shape" begin
            geo = GeoModel()
            add_block(geo, [0, 0, 0], 1, 1, 0, nx=10, ny=10, shape=shape)
            mesh = Mesh(geo)
            @test length(mesh.nodes) == n
        end
    end
end


@announced_testset "Coordinate matrix blocks" begin
    coords4 = [
        0.0  0.0
        2.0  0.0
        2.2  1.0
        0.0  1.0
    ]
    geo4 = GeoModel(quiet=true)
    block4 = add_block(geo4, coords4, nx=2, ny=3, shape=:quad4, tag="quad4-block")
    mesh4 = Mesh(geo4, quiet=true)
    @test block4.blockshape.kind == :quad4
    @test length(mesh4.nodes) == 12
    @test length(mesh4.elems) == 6
    @test all(cell.tag == "quad4-block" for cell in mesh4.elems)

    coords8 = [
        0.0  0.0
        1.0  0.0
        1.0  1.0
        0.0  1.0
        0.5 -0.1
        1.1  0.5
        0.5  1.1
       -0.1  0.5
    ]
    geo8 = GeoModel(quiet=true)
    block8 = add_block(geo8, coords8, nx=2, ny=2, shape=:quad8, tag="quad8-block")
    mesh8 = Mesh(geo8, quiet=true)
    @test block8.blockshape.kind == :quad8
    @test block8.shape.kind == :quad8
    @test length(mesh8.nodes) == 21
    @test length(mesh8.elems) == 4
    @test all(cell.shape.kind == :quad8 for cell in mesh8.elems)
    @test all(cell.tag == "quad8-block" for cell in mesh8.elems)
end


@announced_testset "Mesh refinement" begin
    geo_h = GeoModel(quiet=true)
    add_block(geo_h, [0, 0, 0], 1, 1, 0, nx=1, ny=1, shape=:quad4, tag="href-quad4")
    mesh_h = Mesh(geo_h, quiet=true)
    refined_h = refine(mesh_h, mode=:h, n=2)

    @test length(refined_h.nodes) == 9
    @test length(refined_h.elems) == 4
    @test all(cell.shape.kind == :quad4 for cell in refined_h.elems)
    @test all(cell.role == :solid for cell in refined_h.elems)
    @test all(cell.tag == "href-quad4" for cell in refined_h.elems)

    geo_p = GeoModel(quiet=true)
    add_block(geo_p, [0, 0, 0], 2, 1, 0, nx=2, ny=1, shape=:quad4, tag="pref-quad4")
    mesh_p = Mesh(geo_p, quiet=true)
    refined_p = refine(mesh_p, mode=:p)

    @test length(refined_p.nodes) == 13
    @test length(refined_p.elems) == 2
    @test all(cell.shape.kind == :quad8 for cell in refined_p.elems)
    @test all(cell.role == :solid for cell in refined_p.elems)
    @test all(cell.tag == "pref-quad4" for cell in refined_p.elems)
end


@announced_testset "remove_elements compacts field data" begin
    geo_rm = GeoModel(quiet=true)
    add_block(geo_rm, [0, 0, 0], 2, 1, 0, nx=2, ny=1, shape=:quad4, tag="bulk")
    mesh_rm = Mesh(geo_rm, quiet=true)

    for elem in mesh_rm.elems
        xmid = sum(node.coord.x for node in elem.nodes) / length(elem.nodes)
        elem.tag = xmid < 1.0 ? "left" : "right"
    end

    mesh_rm.node_fields["scalar"] = collect(10:10:10*length(mesh_rm.nodes))
    mesh_rm.node_fields["vec"] = hcat(collect(1.0:length(mesh_rm.nodes)), collect(101.0:100.0+length(mesh_rm.nodes)))
    mesh_rm.elem_fields["scalar"] = collect(1.0:length(mesh_rm.elems))
    mesh_rm.elem_fields["vec"] = hcat(collect(11.0:10.0+length(mesh_rm.elems)), collect(21.0:20.0+length(mesh_rm.elems)))

    kept_elems = select(mesh_rm, :element, "left")
    kept_nodes = unique!(n -> n.id, [node for elem in kept_elems for node in elem.nodes])
    kept_node_ids = getfield.(kept_nodes, :id)
    kept_elem_ids = getfield.(kept_elems, :id)

    expected_node_scalar = mesh_rm.node_fields["scalar"][kept_node_ids]
    expected_node_vec = mesh_rm.node_fields["vec"][kept_node_ids, :]
    expected_elem_scalar = mesh_rm.elem_fields["scalar"][kept_elem_ids]
    expected_elem_vec = mesh_rm.elem_fields["vec"][kept_elem_ids, :]

    remove_elements(mesh_rm, "right")

    @test length(mesh_rm.elems) == length(kept_elems)
    @test length(mesh_rm.nodes) == length(kept_nodes)
    @test mesh_rm.node_fields["scalar"] == expected_node_scalar
    @test mesh_rm.node_fields["vec"] == expected_node_vec
    @test mesh_rm.elem_fields["scalar"] == expected_elem_scalar
    @test mesh_rm.elem_fields["vec"] == expected_elem_vec
    @test mesh_rm.elem_fields["quality"] == [cell.quality for cell in mesh_rm.elems]
end


@announced_testset "Solid blocks in 3d" begin
    nodes_count = [1331, 4961, 9261, 1331, 9261]
    shapes = [:hex8, :hex20, :hex27, :tet4, :tet10]

    for (n, shape) in zip(nodes_count, shapes)
        @announced_testset "$shape" begin
            geo = GeoModel()
            add_block(geo, [0, 0, 0], 1, 1, 1, nx=10, ny=10, nz=10, shape=shape)
            mesh = Mesh(geo)
            @test length(mesh.nodes) == n
        end
    end
end

@announced_testset "Truss meshes" begin
    @announced_testset "2d" begin
        coord = [0 0; 9 0; 18 0; 0 9; 9 9; 18 9.]
        conn = [[1, 2], [1, 5], [2, 3], [2, 6], [2, 5], [2, 4], [3, 6], [3, 5], [4, 5], [5, 6]]
        mesh = Mesh(coord, conn)
        @test length(mesh.nodes) == 6
    end

    @announced_testset "3d" begin
        coord = [0.0 0.0 0.0; 0.0 1.0 0.0; 0.0 1.0 1.0]
        conn = [[1, 3], [1, 2], [2, 3]]
        mesh = Mesh(coord, conn)
        @test length(mesh.nodes) == 3
    end
end

@announced_testset "Meshes with insets" begin
    for shape in [:hex8, :hex20]
        @announced_testset "$shape" begin
            geo = GeoModel()
            add_block(geo, [0, 0, 0], 1, 1, 1, nx=8, ny=8, nz=8, shape=shape)

            p1 = add_point(geo, [0, 0, 0])
            p2 = add_point(geo, [1, 1, 1])
            l1 = add_line(geo, p1, p2)
            add_path(geo, [l1])

            mesh = Mesh(geo)
            @test length(select(mesh.elems, :line)) == 8
        end
    end
end

@announced_testset "Cohesive element generation" begin
    nodes_count = [25, 49, 16, 40, 64, 343, 64, 208]
    elems_count = [84, 39, 21, 21, 432, 432, 81, 81]
    shapes = [:tri3, :tri6, :quad4, :quad8, :tet4, :tet10, :hex8, :hex20]

    for (n, e, shape) in zip(nodes_count, elems_count, shapes)
        @announced_testset "$shape" begin
            shapeobj = get_shape(shape)
            geo = GeoModel()
            if shapeobj.ndim == 2
                add_block(geo, [0, 0, 0], 1, 1, 0, nx=3, ny=3, shape=shape)
            else
                add_block(geo, [0, 0, 0], 1, 1, 1, nx=3, ny=3, nz=3, shape=shape)
            end
            mesh = Mesh(geo)
            add_cohesive_elements(mesh, implicit=true)

            @test length(mesh.nodes) == n
            @test length(mesh.elems) == e
        end
    end
end

@announced_testset "Joint cells with insets" begin
    nodes_count = [38, 23, 77, 71]
    elems_count = [96, 27, 444, 87]
    shapes = [:tri3, :quad4, :tet4, :hex8]

    for (n, e, shape) in zip(nodes_count, elems_count, shapes)
        @announced_testset "$shape" begin
            shapeobj = get_shape(shape)
            geo = GeoModel()
            if shapeobj.ndim == 2
                add_block(geo, [0, 0, 0], 1, 1, 0, nx=3, ny=3, shape=shape)
                p1 = add_point(geo, [0, 0, 0])
                p2 = add_point(geo, [1, 1, 0])
            else
                add_block(geo, [0, 0, 0], 1, 1, 1, nx=3, ny=3, nz=3, shape=shape)
                p1 = add_point(geo, [0, 0, 0])
                p2 = add_point(geo, [1, 1, 1])
            end
            l1 = add_line(geo, p1, p2)
            add_path(geo, [l1])

            mesh = Mesh(geo)
            add_cohesive_elements(mesh, implicit=true)

            @test length(mesh.nodes) == n
            @test length(mesh.elems) == e
        end
    end
end
