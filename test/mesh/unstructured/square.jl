using Serendip, Test

# 2D unstructured mesh
geo = GeoModel()
add_rectangle(geo, [0.0, 0.0, 0.0], 1.0, 1.0; tag="solid")
mesh1 = Mesh(geo)
println(@test mesh1.ctx.ndim == 2)
println(@test length(mesh1.elems) > 0)
println(@test mesh1.elems[1].tag == "solid")
