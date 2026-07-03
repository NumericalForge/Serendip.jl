using Test
using Serendip
using Cairo: read_from_png, width, height, image_surface_get_data


geo = GeoModel(quiet=true)
add_block(geo, [0,0,0], 1, 1, 1, nx=1, ny=1, nz=1, tag="solids")
mesh = Mesh(geo, quiet=true)

mapper = RegionMapper()
add_mapping(mapper, "solids", MechSolid, LinearElastic, E=1e4, nu=0.25)
model = FEModel(mesh, mapper, quiet=true)

plot = DomainPlot(size=(5cm, 4cm))
add_plot(plot, model)
@test isapprox(plot.width, 5cm)
@test isapprox(plot.height, 4cm)
Serendip.save(plot, "output/domainplot-cm.png")
plot_png = read_from_png("output/domainplot-cm.png")
@test Int(width(plot_png)) == round(Int, Serendip._png_raster_scale * plot.width)
@test Int(height(plot_png)) == round(Int, Serendip._png_raster_scale * plot.height)
@test unsafe_load(image_surface_get_data(plot_png)) == 0xffffffff
