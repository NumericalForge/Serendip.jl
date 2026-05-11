using Serendip

geo = GeoModel(size=0.07)

r1 = 0.15
r2 = 0.15

# cylinder surfaces
c1 = add_circle(geo, [0,0,0], [1,0,0], r1)
s1 = extrude(geo, c1, [1,0,0])

c2 = add_circle(geo, [0.5,0,0], [0,0.5,0], r2)
s2 = extrude(geo, c2, [0,0.5,0])

fragment(geo, s1, s2)

mesh = Mesh(geo, quadratic=true)

chart = DomainPlot(elevation=35, azimuth=45)
add_plot(chart, mesh, show_outline=true, outline_width=0.3, outline_color=:black, outline_angle=100)
save(chart, "tubes.pdf")
