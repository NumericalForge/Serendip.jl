using Serendip
using Test

L = 2.0
H = 0.2
W = 0.2

E = 30e6
I = W*H^3/12
P = -1.0
rho = 2.4
A = W*H

# Mesh generation
geo = GeoModel()
add_block(geo, [0.0, 0.0, 0.0], L, W, H, nx=10, ny=1, nz=1, shape=:hex20, tag="solids")
mesh = Mesh(geo)

# Model definition
mapper = RegionMapper()
add_mapping(mapper, "solids", MechSolid, LinearElastic, E=30e6, nu=0.2, rho=rho)

model = FEModel(mesh, mapper)

ana = DynamicAnalysis(model)
# ana = MechAnalysis(model)
log1 = add_logger(ana, :nodalreduce, (x==L), "end.dat")
add_monitor(ana, :nodalreduce, (x==L), :uz)

stage = add_stage(ana, nincs=100, tspan=0.1)
add_bc(stage, :node, (x==0), ux=0, uy=0, uz=0)
add_bc(stage, :face, (x==L), tz=P/(H*W))

# run(ana, alpha=4.2038, beta=174.2803e-6)
run(ana)

t = Float64.(log1.table[:t])
uz = Float64.(log1.table[:uz])

# Step load on an undamped linear oscillator gives u_max = 2*u_static.
u_static = P*L^3/(3*E*I)
u_max_expected = 2*u_static
u_min = minimum(uz)
@test isapprox(u_min, u_max_expected; rtol=0.10)

beta1 = 1.875104068711961
f_expected = beta1^2/(2*pi)*sqrt(E*I/(rho*A*L^4))
f_num = get_frequency(log1.table, :uz; extrema=:minima)
println("Frequency: ", f_num, " Hz")
@test isapprox(f_num, f_expected; rtol=0.15)

chart = Chart(
    xlabel = "Time [s]",
    ylabel = "z displacement",
)
add_line(chart, t, uz, mark=:circle)
save(chart, "cantilever-solid-uz-vs-time.pdf")
