using Serendip

X = collect(0:0.25:2π)

chart1 = Chart(
    title="Sine Curve",
    xlabel=t"$x$",
    ylabel=t"$sin(x)$",
    legend=:outer_right,
    quiet=true,
)
add_line(chart1, X, sin.(X); label=t"$sin(x)$", mark=:circle)

chart2 = Chart(
    title="Bar Values",
    xlabel="Category",
    ylabel="Value",
    legend=:top_left,
    quiet=true,
)
add_bar(chart2, 1:5, [1.0, 1.5, 0.7, 1.2, 1.8]; color=:steelblue, label="bars")

geo = GeoModel(quiet=true)
add_block(geo, [0,0,0], 1,1,1, nx=1, ny=1, nz=1, tag="solids")
mesh = Mesh(geo, quiet=true)

mapper = RegionMapper()
add_mapping(mapper, "solids", MechSolid, LinearElastic, E=1e4, nu=0.25)
model = FEModel(mesh, mapper, quiet=true)
model.node_fields["temp"] = collect(1.0:length(model.nodes))

plot_left = DomainPlot(model,
    field="temp",
    field_kind=:node,
    colorbar=:left,
    quiet=true,
)

plot_top = DomainPlot(model,
    field="temp",
    field_kind=:node,
    colorbar=:top,
    quiet=true,
)

grid = ChartGrid(
    title="Composed Figure With Charts",
    size=(520, 420),
    quiet=true,
)

add_chart(grid, chart1, (1, 1))
add_chart(grid, chart2, (1, 2))
add_chart(grid, plot_left, (2, 1))
add_chart(grid, plot_top, (2, 2))

save(grid, "output/chart-grid.pdf")
save(chart1, "output/chart-grid-child.pdf")
