using Serendip
using Test

printstyled("\nWriting and loading vtk format\n", color=:blue, bold=true)

geo = GeoModel()
add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 1.0, nx=10, ny=10, nz=10, shape=:hex8)
m1 = Mesh(geo)
m1.elems[1].tag = "outer-tag"
m1.elems[end].tag = "inner-tag-123456"

save(m1, "out.vtk")
sleep(0.1)
m2 = Mesh("out.vtk")
node_keys_match = Set(collect(keys(m1.node_fields))) == Set(collect(keys(m2.node_fields)))
elem_keys_match = Set(collect(keys(m1.elem_fields))) == Set(collect(keys(m2.elem_fields)))
t = length(m1.nodes)==length(m2.nodes) &&
    length(m1.elems)==length(m2.elems) &&
    node_keys_match &&
    elem_keys_match &&
    haskey(m2.elem_fields, "tag-data") &&
    m2.elems[1].tag == "outer-tag" &&
    m2.elems[end].tag == "inner-tag-123456"

TR = @test t
println(TR)

save(m1, "out.vtu")
sleep(0.1)
m2 = Mesh("out.vtu")
node_keys_match = Set(collect(keys(m1.node_fields))) == Set(collect(keys(m2.node_fields)))
elem_keys_match = Set(collect(keys(m1.elem_fields))) == Set(collect(keys(m2.elem_fields)))
t = length(m1.nodes)==length(m2.nodes) &&
    length(m1.elems)==length(m2.elems) &&
    node_keys_match &&
    elem_keys_match &&
    haskey(m2.elem_fields, "tag-data") &&
    m2.elems[1].tag == "outer-tag" &&
    m2.elems[end].tag == "inner-tag-123456"

TR = @test t
println(TR)
