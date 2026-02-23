using Serendip, Test

geo = GeoModel()

p1 = add_point(geo, [0.0, 0.0, 0.0])
p2 = add_point(geo, [1.0, 0.0, 0.0])
p3 = add_point(geo, [1.0, 0.0, 0.75])
p4 = add_point(geo, [0.0, 0.0, 1.0])

l1 = add_line(geo, p1, p2)
l2 = add_line(geo, p2, p3)
l3 = add_line(geo, p3, p4)
l4 = add_line(geo, p4, p1)

loop = add_loop(geo, [l1, l2, l3, l4])
surf = add_plane_surface(geo, loop, tag="body")

revolve(geo, [surf], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2pi)
select(geo, :volume, [0.5, 0.0, 0.5], tag="body")

mesh = Mesh(geo)
println(@test mesh.ctx.ndim == 3)
println(@test mesh.elems[1].tag == "body")
