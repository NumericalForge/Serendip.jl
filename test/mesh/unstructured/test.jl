using Serendip, Test

geo = GeoModel()

outer = add_rectangle(geo, [0.0, 0.0, 0.0], 1.0, 1.0)
inner = add_rectangle(geo, [0.25, 0.25, 0.0], 0.5, 0.5)
surfs = cut(geo, outer, inner; remove_object=true, remove_tool=true)

for surf in surfs
    extrude(geo, surf, [0.0, 0.5, 0.0])
end

mesh = Mesh(geo)
println(@test mesh.ctx.ndim == 3)
println(@test length(mesh.elems) > 0)
