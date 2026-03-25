using Serendip

geo = GeoModel(quiet=true)
add_block(geo, [0,0,0], 1,1,1, nx=1, ny=1, nz=1, tag="solids")
mesh = Mesh(geo, quiet=true)

mapper = RegionMapper()
add_mapping(mapper, "solids", MechSolid, LinearElastic, E=1e4, nu=0.25)
model = FEModel(mesh, mapper, quiet=true)
model.node_fields["temp"] = collect(1.0:length(model.nodes))

for loc in (:left, :right, :top, :bottom)
    fig = DomainPlot(model,
        field = "temp",
        field_kind = :node,
        colorbar = loc,
        quiet = true,
    )
    save(fig, "output/domainplot-colorbar-$(loc).pdf")
end
