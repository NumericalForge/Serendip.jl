using Serendip

L = 2.0
H = 0.2
W = 0.2

E = 30e6
I = W*H^3/12

# Mesh generation
geo = GeoModel()
add_block(geo, [0.0, 0.0, 0.0], L, W, H, nx=10, ny=1, nz=1, shape=HEX20, tag="solids")
mesh = Mesh(geo)

# Model definition
mapper = RegionMapper()
add_mapping(mapper, "solids", MechBulk, LinearElastic, E=30e6, nu=0.2, rho=2.4)

model = FEModel(mesh, mapper)

ana = DynamicAnalysis(model)
# ana = MechAnalysis(model)
log1 = add_logger(ana, :nodalreduce, (x==L), "end.dat")
add_monitor(ana, :nodalreduce, (x==L), :uz)

stage = add_stage(ana, nincs=100, tspan=0.1)
add_bc(stage, :node, (x==0), ux=0, uy=0, uz=0)
add_bc(stage, :face, (x==L), tz=-1)

# run(ana, alpha=4.2038, beta=174.2803e-6)
run(ana)

chart = Chart(
    xlabel = "Time [s]",
    ylabel = "z displacement",
)
add_line(chart, log1.table[:t], log1.table[:uz], mark=:circle)
save(chart, "dynamic-uz-vs-time.pdf")