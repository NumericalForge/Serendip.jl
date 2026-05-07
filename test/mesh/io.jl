using Serendip
using Test
using JSON

function io_test_mesh()
    geo = GeoModel()
    add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 1.0, nx=2, ny=2, nz=2, shape=:hex8)
    mesh = Mesh(geo, quiet=true)
    mesh.elems[1].tag = "outer-tag"
    mesh.elems[end].tag = "inner-tag-123456"
    return mesh
end

function test_round_trip(mesh::Mesh, filename::String)
    @test save(mesh, filename, quiet=true) === nothing

    saved = Mesh(filename)
    @test length(saved.nodes) == length(mesh.nodes)
    @test length(saved.elems) == length(mesh.elems)
    @test Set(keys(saved.node_fields)) == Set(keys(mesh.node_fields))
    @test Set(keys(saved.elem_fields)) == Set(keys(mesh.elem_fields))
    @test haskey(saved.elem_fields, "tag-data")
    @test saved.elems[1].tag == "outer-tag"
    @test saved.elems[end].tag == "inner-tag-123456"
end

@testset "Mesh IO" begin
    dir = "."
    # mktempdir() do dir
        @testset "VTK round-trip" begin
            test_round_trip(io_test_mesh(), joinpath(dir, "mesh.vtk"))
        end

        @testset "VTU round-trip" begin
            test_round_trip(io_test_mesh(), joinpath(dir, "mesh.vtu"))
        end

        @testset "JSON output" begin
            mesh = io_test_mesh()
            filename = joinpath(dir, "mesh.json")
            json = mesh_json(mesh)
            @test save(mesh, filename, quiet=true) === nothing

            data = JSON.parse(json)
            saved_data = JSON.parsefile(filename)

            @test saved_data == data
            @test !haskey(data, "mesh")
            @test length(data["nodes"]) == length(mesh.nodes)
            @test length(data["elements"]) == length(mesh.elems)
            @test length(data["faces"]) == length(mesh.faces)
            @test data["nodes"][1] == [0.0, 0.0, 0.0]
            @test data["elements"][1]["shape"] == "hex8"
            @test data["elements"][1]["conn"] == [node.id for node in mesh.elems[1].nodes]

            node_fields = data["nodal_fields"]
            elem_fields = data["element_fields"]

            @test node_fields["node-id"] == [node.id for node in mesh.nodes]
            @test elem_fields["elem-id"] == [elem.id for elem in mesh.elems]
            @test elem_fields["cell-type"] == [Int(cell.shape.vtk_type) for cell in mesh.elems]
            @test elem_fields["tag"][1] != elem_fields["tag"][end]
        end
    # end
end
