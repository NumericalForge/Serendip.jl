using Serendip
using Test

geo = GeoModel()

outer = add_rectangle(geo, [0.0, 0.0, 0.0], 1.0, 1.0, tag="outer")

p1 = add_point(geo, [0.5, 0.0, 0.0])
p2 = add_point(geo, [0.5, 1.0, 0.0])
e1 = add_line(geo, p1, p2)

p3 = add_point(geo, [0.0, 0.75, 0.0])
p4 = add_point(geo, [1.0, 0.75, 0.0])
e2 = add_line(geo, p3, p4)

out = fragment(geo, [outer], [e1, e2]; remove_object=true, remove_tool=true)
nsurfs = length([e for e in out if e isa Surface])

mesh = Mesh(geo)

println(@test nsurfs >= 3)
println(@test mesh.elems[1].tag == "outer")
