using Serendip, Test

geo = GeoModel()
base = add_rectangle(geo, [0.0, 0.0, 0.0], 1.0, 1.0; tag="body")
extrude(geo, base, [0.0, 0.0, 1.0])

mesh = Mesh(geo)
println(@test mesh.ctx.ndim == 3)
println(@test length(mesh.elems) > 0)
println(@test mesh.elems[1].tag == "body")
