using Test
using Serendip
using Cairo: read_from_png, width, height, image_surface_get_data

chart = Chart(size=(5cm, 4cm), quiet=true)
@test isapprox(chart.width, 5cm)
@test isapprox(chart.height, 4cm)

add_line(chart, 0:1, 0:1; label="line")
Serendip.save(chart, "output/chart-cm.pdf")
@test isfile("output/chart-cm.pdf")
Serendip.save(chart, "output/chart-cm.png")
chart_png = read_from_png("output/chart-cm.png")
@test Int(width(chart_png)) == round(Int, Serendip._png_raster_scale * chart.width)
@test Int(height(chart_png)) == round(Int, Serendip._png_raster_scale * chart.height)
@test unsafe_load(image_surface_get_data(chart_png)) == 0xffffffff

grid = ChartGrid(size=(12cm, 8cm), quiet=true)
@test isapprox(grid.width, 12cm)
@test isapprox(grid.height, 8cm)
add_chart(grid, chart, (1, 1))
Serendip.save(grid, "output/chart-grid-cm.png")
grid_png = read_from_png("output/chart-grid-cm.png")
@test Int(width(grid_png)) == round(Int, Serendip._png_raster_scale * grid.width)
@test Int(height(grid_png)) == round(Int, Serendip._png_raster_scale * grid.height)

geo = GeoModel(quiet=true)
add_block(geo, [0,0,0], 1, 1, 1, nx=1, ny=1, nz=1, tag="solids")
mesh = Mesh(geo, quiet=true)

mapper = RegionMapper()
add_mapping(mapper, "solids", MechSolid, LinearElastic, E=1e4, nu=0.25)
model = FEModel(mesh, mapper, quiet=true)

plot = DomainPlot(model, size=(5cm, 4cm), quiet=true)
@test isapprox(plot.width, 5cm)
@test isapprox(plot.height, 4cm)
Serendip.save(plot, "output/domainplot-cm.png")
plot_png = read_from_png("output/domainplot-cm.png")
@test Int(width(plot_png)) == round(Int, Serendip._png_raster_scale * plot.width)
@test Int(height(plot_png)) == round(Int, Serendip._png_raster_scale * plot.height)
@test unsafe_load(image_surface_get_data(plot_png)) == 0xffffffff
