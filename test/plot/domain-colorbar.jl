using Serendip
using Test

geo = GeoModel(quiet=true)
add_block(geo, [0, 0, 0], 1, 1, 1, nx=1, ny=1, nz=1, tag="solids")
mesh = Mesh(geo, quiet=true)

mapper = RegionMapper()
add_mapping(mapper, "solids", MechSolid, LinearElastic, E=1e4, nu=0.25)
model = FEModel(mesh, mapper, quiet=true)
model.node_fields["temp"] = collect(1.0:length(model.nodes))

grid = ChartGrid(
    title="DomainPlot Colorbar Locations",
    size=(16cm, 12cm),
)

add_chart(grid, DomainPlot(model, field="temp", field_kind=:node, colorbar=:left), (1, 1))
add_chart(grid, DomainPlot(model, field="temp", field_kind=:node, colorbar=:right), (1, 2))
add_chart(grid, DomainPlot(model, field="temp", field_kind=:node, colorbar=:top), (2, 1))
add_chart(grid, DomainPlot(model, field="temp", field_kind=:node, colorbar=:bottom), (2, 2))

save(grid, "output/domainplot.pdf")
save(grid, "output/domainplot.png")

@test isfile("output/domainplot.pdf")
@test isfile("output/domainplot.png")
