using Serendip, Test

# 2D unstructured mesh with a hole
geo = GeoModel()

outer = add_rectangle(geo, [0.0, 0.0, 0.0], 1.0, 1.0, tag="outer")
inner = add_rectangle(geo, [0.25, 0.25, 0.0], 0.5, 0.5, tag="inner")
cut(geo, outer, inner; remove_object=true, remove_tool=true, tag="outer")

mesh = Mesh(geo, quadratic=true)
@test mesh.ctx.ndim == 2
@test length(mesh.elems) > 0
@test mesh.elems[1].tag == "outer"
