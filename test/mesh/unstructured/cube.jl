using Serendip, Test

geo = GeoModel()
base = add_rectangle(geo, [0.0, 0.0, 0.0], 1.0, 1.0; tag="body")
extrude(geo, base, [0.0, 0.0, 1.0])

mesh = Mesh(geo)
@test mesh.ctx.ndim == 3
@test length(mesh.elems) > 0
@test mesh.elems[1].tag == "body"
