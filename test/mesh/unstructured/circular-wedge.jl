using Serendip, Test

geo = GeoModel()

p0 = add_point(geo, [0.0, 0.0, 0.0])
p1 = add_point(geo, [1.0, 0.0, 0.0])
p2 = add_point(geo, [0.0, 1.0, 0.0])
p3 = add_point(geo, [-1.0, 0.0, 0.0])

a1 = add_circle_arc(geo, p1, p0, p2)
a2 = add_circle_arc(geo, p2, p0, p3)
l1 = add_line(geo, p3, p1)

loop = add_loop(geo, [a1, a2, l1])
surf = add_plane_surface(geo, loop, tag="body")
extrude(geo, surf, [0.0, 0.0, 0.5])

mesh = Mesh(geo)
println(@test mesh.ctx.ndim == 3)
println(@test length(mesh.elems) > 0)
println(@test mesh.elems[1].tag == "body")
