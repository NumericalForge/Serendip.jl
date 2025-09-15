using Serendip

geo = GeoModel()
add_block(geo, [0,0,0], [1,1,1], nx=2, ny=2, nz=4, tag="solids")
mesh = Mesh(geo, quiet=true)

mapper = RegionMapper()
add_mapping(mapper, "solids", MechBulk, LinearElastic, E=1e4, nu=0.25)
model = FEModel(mesh, mapper)

ana = MechAnalysis(model)
stage = add_stage(ana)

add_bc(stage, :node, z==0, ux=0, uy=0, uz=0 )
add_bc(stage, :face, z==1, tz=-100.0)

run(ana)

cmap = Colormap(:coolwarm, rev=true)

# plotting
plot = DomainPlot(model,
    azimuth = 30,
    elevation = 30,
    field = "uz",
    colormap = cmap,
    # warp = 20,
    # label = L"u_z",
    # fontsize = 8,
    # font = "Times New Roman",
    # colorbarscale = 0.8
)

save(plot, "mesh.pdf")