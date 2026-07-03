using Serendip
using Test

# ❱❱❱ Geometry and mesh

ℓ = 0.1
h = 0.2

geo = GeoModel()
add_block(geo, [0, 0, 0], ℓ, ℓ, h, shape=:hex8, nx=1, ny=1, nz=1)
mesh = Mesh(geo, quadratic=true)
select(mesh, :element, :bulk, tag="solids")

# ❱❱❱ Elements and constitutive models

fc = -30.0e3
ft = 3.0e3

mapper = RegionMapper()
add_mapping( mapper, "solids", MechSolid, ESCP, E=30.0e6, nu= 0.25, beta=1.15, fc= fc, epsc=-0.002, ft=ft, wc=0.0001 )

# Uniaxial compression

printstyled("\nUniaxial compression:\n\n", color=:yellow, bold=true)
model = FEModel(mesh, mapper)
ana = MechAnalysis(model)
log_c = add_logger(ana, :ip, [ℓ/2, ℓ/2, h/2])
log_c_top = add_logger(ana, :nodalreduce, z==h)

stage = add_stage(ana)
add_bc(stage, :node, x == 0, ux=0)
add_bc(stage, :node, y == 0, uy=0)
add_bc(stage, :node, z == 0, uz=0)
add_bc(stage, :node, z == h, uz=-0.0015)

status = run(ana, autoinc=true, tol=0.1, rspan=0.02, dT0=0.002, quiet=false)
@test status.successful
@test minimum(log_c.table["σzz"]) < 0.5*fc
@test log_c.table["εcp"][end] > 0.0

# Biaxial compression

printstyled("\nBiaxial compression:\n\n", color=:yellow, bold=true)
model = FEModel(mesh, mapper)
ana = MechAnalysis(model)
log_b = add_logger(ana, :ip, [ℓ/2, ℓ/2, h/2])
log_b_xl = add_logger(ana, :nodalreduce, x==ℓ)
log_b_yl = add_logger(ana, :nodalreduce, y==ℓ)

stage = add_stage(ana)
ub = -0.0006
add_bc(stage, :node, x == 0, ux=0)
add_bc(stage, :node, y == 0, uy=0)
add_bc(stage, :node, z == 0, uz=0)
add_bc(stage, :node, x == ℓ, ux=ub)
add_bc(stage, :node, y == ℓ, uy=ub)

status = run(ana, autoinc=true, tol=0.1, rspan=0.02, dT0=0.002, quiet=false)
@test status.successful
@test minimum(log_b.table["σxx"]) < 0.5*fc
@test minimum(log_b.table["σyy"]) < 0.5*fc
@test log_b.table["εcp"][end] > 0.0

# Tension

printstyled("\nUniaxial tension:\n\n", color=:yellow, bold=true)
model = FEModel(mesh, mapper)
ana = MechAnalysis(model)
log_t = add_logger(ana, :ip, [ℓ/2, ℓ/2, h/2])
log_t_top = add_logger(ana, :nodalreduce, z==h)

stage = add_stage(ana, nouts=20)
add_bc(stage, :node, x == 0, ux=0)
add_bc(stage, :node, y == 0, uy=0)
add_bc(stage, :node, z == 0, uz=0)
add_bc(stage, :node, z == h, uz=0.00035)

status = run(ana, autoinc=true, tol=0.01, dT0=0.001, dTmax=0.004, quiet=false)
@test status.successful
@test maximum(log_t.table["σzz"]) > 0.5*ft
@test log_t.table["εtp"][end] > 0.0
