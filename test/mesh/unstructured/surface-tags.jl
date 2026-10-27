using Serendip, Test


function tagged_surface_patch_mesh(; quadratic=false)
    geo = GeoModel(size=0.2, quiet=true)
    volume = add_box(
        geo,
        [0.0, 0.0, 0.0],
        1.0,
        1.0,
        1.0,
        tag="concrete",
    )

    points = [
        add_point(geo, [0.2, 0.3, 0.0]),
        add_point(geo, [0.8, 0.3, 0.0]),
        add_point(geo, [0.8, 0.7, 0.0]),
        add_point(geo, [0.2, 0.7, 0.0]),
    ]
    patch = add_polygon(geo, points, tag="EBR-face")
    fragment(geo, volume, patch)

    return Mesh(geo, quadratic=quadratic, quiet=true)
end


function assert_surface_patch_tags(mesh)
    tagged_faces = select(mesh, :face, "EBR-face")
    @test !isempty(tagged_faces)
    @test all(face -> face.owner.tag == "concrete", tagged_faces)
    @test all(tagged_faces) do face
        all(face.nodes) do node
            x, y, z = node.coord
            abs(z) <= 1e-8 &&
                0.2 - 1e-8 <= x <= 0.8 + 1e-8 &&
                0.3 - 1e-8 <= y <= 0.7 + 1e-8
        end
    end

    bottom_faces = select(mesh, :face, :(abs(z) <= 1e-8))
    @test length(tagged_faces) < length(bottom_faces)
    @test all(elem -> elem.tag == "concrete", select(mesh, :element, :solid))

    return Set(Serendip._cell_key(face) for face in tagged_faces)
end


@announced_testset "Geometric surface tags on boundary faces" begin
    @announced_testset "Linear and quadratic meshes" begin
        for quadratic in (false, true)
            mesh = tagged_surface_patch_mesh(quadratic=quadratic)
            expected_keys = assert_surface_patch_tags(mesh)

            Serendip.synchronize(mesh)
            @test Set(
                Serendip._cell_key(face)
                for face in select(mesh, :face, "EBR-face")
            ) == expected_keys
        end
    end

    @announced_testset "Copy and join preserve boundary tags" begin
        mesh = tagged_surface_patch_mesh()
        expected_keys = assert_surface_patch_tags(mesh)

        copied = copy(mesh)
        @test Set(
            Serendip._cell_key(face)
            for face in select(copied, :face, "EBR-face")
        ) == expected_keys

        other = Mesh(
            [
                2.0 0.0 0.0
                3.0 0.0 0.0
                3.0 1.0 0.0
                2.0 1.0 0.0
                2.0 0.0 1.0
                3.0 0.0 1.0
                3.0 1.0 1.0
                2.0 1.0 1.0
            ],
            [[1, 2, 3, 4, 5, 6, 7, 8]],
            [:hex8],
            tag="other",
            quiet=true,
        )
        joined = join_meshes(mesh, other)
        @test Set(
            Serendip._cell_key(face)
            for face in select(joined, :face, "EBR-face")
        ) == expected_keys
    end

    @announced_testset "Tags survive cohesive insertion and select boundary shells" begin
        mesh = tagged_surface_patch_mesh()
        expected_keys = assert_surface_patch_tags(mesh)

        add_cohesive_elements(mesh, "concrete", tag="cohesive", quiet=true)
        @test Set(
            Serendip._cell_key(face)
            for face in select(mesh, :face, "EBR-face")
        ) == expected_keys

        add_boundary_shell_elements(
            mesh,
            "EBR-face",
            tag="EBR-shell",
            contact_tag="EBR-interface",
            quiet=true,
        )
        @test length(select(mesh, :element, :surface, "EBR-shell")) == length(expected_keys)
        @test length(select(mesh, :element, :contact, "EBR-interface")) == length(expected_keys)
    end

    @announced_testset "Conflicting imported surface tags" begin
        tags = Dict{Tuple,String}()
        key = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
        Serendip._set_boundary_face_tag!(tags, key, "first", 1)
        @test_throws ErrorException Serendip._set_boundary_face_tag!(tags, key, "second", 2)
    end
end
