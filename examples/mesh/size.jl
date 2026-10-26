using Serendip

# ❱❱❱ Geometry model ❰❰❰ #

ℓ = 0.4   # x direction
b = 0.025 # y direction
h = 0.1   # z direction

geo = GeoModel()
p1 = add_point(geo, [0, 0, 0], size=0.01)
p2 = add_point(geo, [ℓ, 0, 0], size=0.01)
p3 = add_point(geo, [ℓ, h, 0], size=0.01)
p4 = add_point(geo, [0.25*ℓ, h, 0], size=0.005)
p5 = add_point(geo, [0, h, 0], size=0.005)
p6 = add_point(geo, [0, h/2, 0], size=0.005)

add_polygon(geo, [p1, p2, p3, p4, p5, p6], tag="solid")

# ❱❱❱ Mesh generation ❰❰❰ #

mesh = Mesh(geo)

# ❱❱❱ Plotting ❰❰❰ #

plot = DomainPlot()
add_plot(plot, mesh)
save(plot, "size.pdf")

