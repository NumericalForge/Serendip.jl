
using Serendip

# ❱❱❱ Geometry model ❰❰❰ #

ℓ = 0.4   # x direction
b = 0.025 # y direction
h = 0.1   # z direction

geo = GeoModel()
add_rectangle(geo, [0,0,0], ℓ, h; tag="solid")

p1 = add_point(geo, [0.00*ℓ, 0.95*h, 0])
p2 = add_point(geo, [0.95*ℓ, 0.95*h, 0])
bar = add_line(geo, p1, p2)
add_path(geo, [bar], tag="bar", interface_tag="bar-interface")

# Refinement
set_refinement(geo, [0, h, 0], ℓ/2, h, 0.1, 0.005, 0.01; transition=0.5, roundness=0.8)

# Set the size at a particular point
pts = select(geo, :point, [ℓ, h, 0])
set_size(geo, pts, 0.005)

# ❱❱❱ Mesh generation ❰❰❰ #

mesh = Mesh(geo)

# ❱❱❱ Plotting ❰❰❰ #

plot = DomainPlot()
add_plot(plot, mesh)
save(plot, "refinement.pdf")

