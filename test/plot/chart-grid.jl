using Serendip
using Test

X = collect(0:0.25:2π)

chart1 = Chart(
    title="Sine Curve",
    background=:white,
    legend_background=:old_paper,
    xlabel="`x`",
    ylabel="`sin(x)`",
    legend=:outer_right,
    quiet=true,
)
add_line(chart1, X, sin.(X); label="`sin(x)`", mark=:circle)

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
    title="Composed Figure With `sigma_n` Charts",
    size=(18cm, 14cm),
    background=:old_paper,
    column_headers=["`sin(x)`", "Bar Plot", "Ignored"],
    row_headers=["`u_x`", "Temperature", "Ignored"],
    quiet=true,
)

add_chart(grid, chart1, (1, 1))
add_chart(grid, chart2, (1, 2))
add_chart(grid, plot_left, (2, 1))
add_chart(grid, plot_top, (2, 2))

Serendip.configure!(grid)
temp_header_width = Serendip.getsize("Temperature", grid.font_size)[1]
@test grid.row_header_boxes[2].frame.width < temp_header_width

save(grid, "output/chart-grid.pdf")
save(grid, "output/chart-grid.png")
save(chart1, "output/chart-grid-child.pdf")
@test chart1.background == Color(:white)
