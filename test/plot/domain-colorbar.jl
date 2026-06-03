using Serendip
using Test

geo = GeoModel(quiet=true)
add_block(geo, [0, 0, 0], 1, 1, 1, nx=1, ny=1, nz=1, tag="solids")
mesh = Mesh(geo, quiet=true)

mapper = RegionMapper()
add_mapping(mapper, "solids", MechSolid, LinearElastic, E=1e4, nu=0.25)
model = FEModel(mesh, mapper, quiet=true)
model.node_fields["temp"] = collect(1.0:length(model.nodes))

function colorbar_plot(model, location)
    plot = DomainPlot()
    add_plot(plot, model; field="temp", field_kind=:node, colorbar=location)
    return plot
end

for location in (:left, :right, :top, :bottom)
    plot = colorbar_plot(model, location)
    save(plot, "output/domainplot-colorbar-$(location).pdf")
    save(plot, "output/domainplot-colorbar-$(location).png")
    @test isfile("output/domainplot-colorbar-$(location).pdf")
    @test isfile("output/domainplot-colorbar-$(location).png")
end

multi_plot = DomainPlot()
add_plot(multi_plot, model; field="temp", field_kind=:node, colorbar=:bottom, label="a")
add_plot(multi_plot, model; field="temp", field_kind=:node, colorbar=:bottom, label="b")
Serendip.configure!(multi_plot)

cb1, cb2 = multi_plot.bottom_items
frame_gap = cb2.frame.x - (cb1.frame.x + cb1.frame.width)
@test isapprox(frame_gap, 0.05 * multi_plot.canvas.frame.width)
