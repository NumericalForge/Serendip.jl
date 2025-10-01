using Serendip

# ❱❱❱ Geometry and generation

geo = GeoModel()
add_box(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 1.0)
mesh = Mesh(geo)
select(mesh, :element, tag="solids")

# ❱❱❱ Finite element modeling

mapper= RegionMapper()
add_mapping(mapper, "solids", MechBulk, LinearElastic, E=2e3, nu=0.2)

model = FEModel(mesh, mapper)
ana = MechAnalysis(model)
add_logger(ana, :nodalreduce, (z==1), "top-face.dat")

stage = add_stage(ana, nincs=4, nouts=1)
add_bc(stage, :node, (z==0), ux=0, uy=0, uz=0)
add_bc(stage, :face, (z==1), tz=:(-10*x))   # triangular load

run(ana)

# ❱❱❱ Post-processing

plot = DomainPlot(model,
    field = "σzz",
    colormap = :coolwarm,
    label = L"\sigma_z",
    elevation = 30,
    azimuth = -60,
    warp = 50
)
save(plot, "3d-static.pdf")