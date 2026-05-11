using Serendip
geo = GeoModel(size=0.15)
p1 = add_point(geo, [0,0,-0.5])
p2 = add_point(geo, [0,0,0])
p3 = add_point(geo, [1,0,-0.5])
p4 = add_point(geo, [2,0,0])
p5 = add_point(geo, [2,0,-0.5])

l1 = add_line(geo, p1, p2)
a1 = add_circle_arc(geo, p2, p3, p4)
l2 = add_line(geo, p4, p5)

extrude(geo, [l1, a1, l2], [0,3,0], recombine=true)

lines = select(geo, :edge)
set_transfinite_curve(geo, lines, 12)

vlines = select(geo, :edge, or(y==0,y==3), z<=0)
set_transfinite_curve(geo, vlines, 5)

surfs = select(geo, :surface)
set_transfinite_surface(geo, surfs)
set_recombine(geo, surfs)


mesh = Mesh(geo, quadratic=true)
chart = DomainPlot(azimuth=60, elevation=30)
add_plot(chart, mesh, show_outline=true, outline_width=0.3, outline_color=:black, outline_angle=90)
save(chart, "shell-mesh-cover.pdf")