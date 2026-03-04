using Serendip

geo = GeoModel()
add_block(geo, [0.0, 0.0, 0.0], 0.2, 2.0, 0.2, nx=1, ny=8, nz=1, tag="solids")
mesh = Mesh(geo)

mapper = RegionMapper()
add_mapping(mapper, "solids", MechBulk, LinearElastic, E=36e6, nu=0.2, rho=24.0)
model = FEModel(mesh, mapper)
ana = DynamicAnalysis(model)

# Track displacement at one load application point.
load_log = add_logger(ana, :node, (x==0.2, y==1, z==0.2), "load-node.dat")

stage = add_stage(ana, nincs=100, tspan=0.1)
add_bc(stage, :node, (y==0, z==0), ux=0, uy=0, uz=0)
add_bc(stage, :node, (y==2, z==0), uz=0)
add_bc(stage, :node, (y==1, z==0.2), fz=-50)

run(ana, dTmax=0.05)

chart = Chart(
    xlabel = "Time [s]",
    ylabel = "z displacement",
)
add_line(chart, load_log.table[:t], load_log.table[:uz], mark=:circle)
save(chart, "dynamic-uz-vs-time.pdf")
