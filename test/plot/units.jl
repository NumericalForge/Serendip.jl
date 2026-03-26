using Test
using Serendip

chart = Chart(size=(4cm, 3cm), quiet=true)
@test isapprox(chart.width, 4cm)
@test isapprox(chart.height, 3cm)

add_line(chart, 0:1, 0:1; label="line")
save(chart, "output/chart-cm.pdf")
@test isfile("output/chart-cm.pdf")

grid = ChartGrid(size=(12cm, 8cm), quiet=true)
@test isapprox(grid.width, 12cm)
@test isapprox(grid.height, 8cm)

geo = GeoModel(quiet=true)
add_block(geo, [0,0,0], 1, 1, 1, nx=1, ny=1, nz=1, tag="solids")
mesh = Mesh(geo, quiet=true)

mapper = RegionMapper()
add_mapping(mapper, "solids", MechSolid, LinearElastic, E=1e4, nu=0.25)
model = FEModel(mesh, mapper, quiet=true)

plot = DomainPlot(model, size=(5cm, 4cm), quiet=true)
@test isapprox(plot.width, 5cm)
@test isapprox(plot.height, 4cm)
