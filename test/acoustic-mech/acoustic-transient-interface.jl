using Serendip
using Test

h    = 1.0
c    = 1500.0
tmax = 0.1
q0   = -1.0e-5

geo = GeoModel()
add_block(geo, [0.0, 0.0], 0.1, h, 0.0, nx=1, ny=10, shape=QUAD4, tag="wall")
add_block(geo, [0.1, 0.0], 10.0, h, 0.0, nx=20, ny=10, shape=QUAD4, tag="fluid")

mesh = Mesh(geo)
add_contact_elements(mesh, "wall", "fluid", tag="interface")

mapper = RegionMapper()
add_mapping(mapper, "wall", MechSolid, LinearElastic, E=211e6, nu=0.3, rho=7.8)
add_mapping(mapper, "fluid", AcousticFluid, LinearAcousticFluid, rho=1.0, c=c)
add_mapping(mapper, "interface", AcousticMechInterface, AcousticMechInterfaceCoupling)

model = FEModel(mesh, mapper, thickness=1.0, g=9.81)
ana   = AcousticMechAnalysis(model)

log_wall  = add_logger(ana, :node, [0.0, h / 2, 0.0])
log_fluid = add_logger(ana, :node, [5.1, h / 2, 0.0])

stage = add_stage(ana, nincs=100, nouts=100, tspan=tmax)
add_bc(stage, :face, (x <= 0.1, y == 0.0), ux=0.0, uy=0.0)
add_bc(stage, :face, (x >= 0.1, y == h), fs=true)
add_bc(stage, :face, (x == 10.1), tq=:(t <= $tmax / 10 ? $q0 : 0.0))

status = run(ana, autoinc=true, quiet=false, tol=1e-2, dT0=0.01)

@test status.successful

t_wall = Float64.(log_wall.table[:t])
ux_wall = Float64.(log_wall.table[:ux])
t_fluid = Float64.(log_fluid.table[:t])
up_fluid = Float64.(log_fluid.table[:up])

@test !isempty(t_wall)
@test !isempty(t_fluid)
@test isapprox(t_wall[end], tmax; rtol=1e-8, atol=1e-12)
@test isapprox(t_fluid[end], tmax; rtol=1e-8, atol=1e-12)
# @test maximum(abs.(ux_wall)) > 0.0
# @test maximum(abs.(up_fluid)) > 0.0
