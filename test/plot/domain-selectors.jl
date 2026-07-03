using Serendip
using Test

function make_domain_selector_mesh()
    geo = GeoModel(quiet=true)
    add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0, nx=2, ny=1, shape=:quad4, tag="solids")
    mesh = Mesh(geo, quiet=true)

    for cell in select(mesh, :element, :solid)
        xmid = sum(node.coord.x for node in cell.nodes) / length(cell.nodes)
        cell.tag = xmid < 0.5 ? "left" : "right"
    end

    mesh.node_fields["temp"] = collect(1.0:length(mesh.nodes))
    return mesh
end

solid_tags(layer) = Set(cell.tag for cell in layer.elems if cell.role == :solid)

mesh = make_domain_selector_mesh()

compat_plot = DomainPlot(mesh, "left"; field="temp", field_kind=:node, colorbar=:left)
Serendip.configure!(compat_plot)
@test solid_tags(compat_plot.layers[1]) == Set(["left"])

left_plot = DomainPlot()
add_plot(left_plot, mesh, "left"; field="temp", field_kind=:node, colorbar=:left)
Serendip.configure!(left_plot)
@test solid_tags(left_plot.layers[1]) == Set(["left"])
@test length(left_plot.colorbars) == 1
@test length(left_plot.left_items) == 1
@test !isempty(left_plot.layers[1].values)

right_plot = DomainPlot()
add_plot(right_plot, mesh, none_of("left"); field="temp", field_kind=:node, colorbar=:right)
Serendip.configure!(right_plot)
@test solid_tags(right_plot.layers[1]) == Set(["right"])
@test length(right_plot.right_items) == 1

multi_plot = DomainPlot()
add_plot(multi_plot, mesh, "left"; field="temp", field_kind=:node, colorbar=:left, label="left")
add_plot(multi_plot, mesh, "right"; field="temp", field_kind=:node, colorbar=:right, label="right")
Serendip.configure!(multi_plot)
@test solid_tags(multi_plot.layers[1]) == Set(["left"])
@test solid_tags(multi_plot.layers[2]) == Set(["right"])
@test length(multi_plot.colorbars) == 2
