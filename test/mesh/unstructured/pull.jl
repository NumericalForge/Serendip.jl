using Serendip, Test

geo = GeoModel()

base = add_rectangle(geo, [0.0, 0.0, 0.0], 2.0, 2.0, tag="body")
extrude(geo, base, [0.0, 0.0, 1.0])
select(geo, :volume, [1.0, 1.0, 0.5], tag="body")

p1 = add_point(geo, [1.0, 0.0, 1.0])
p2 = add_point(geo, [1.0, 1.0, 1.0])
p3 = add_point(geo, [0.0, 1.0, 1.0])
add_line(geo, p1, p2)
add_line(geo, p2, p3)

top = select(geo, :surface, [0.5, 0.5, 1.0])
extrude(geo, top, [0.0, 0.0, 0.7])
select(geo, :volume, [0.5, 0.5, 1.2], tag="body")

mesh = Mesh(geo)
println(@test mesh.ctx.ndim == 3)
println(@test mesh.elems[1].tag == "body")
