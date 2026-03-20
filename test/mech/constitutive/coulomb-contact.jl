using Serendip
using Test

E = 27.0e6
nu = 0.2
kn = 2.0e8
ks = 1.8e8
mu = 0.9

geo = GeoModel()
add_block(geo, [0.0, 0.0], 0.1, 0.1, 0.0; nx=1, ny=1, shape=:quad4, tag="bulk")
add_block(geo, [0.1, 0.0], 0.1, 0.1, 0.0; nx=1, ny=1, shape=:quad4, tag="bulk")
mesh = Mesh(geo)

# Split regions so contact is created between left and right blocks.
select(mesh, :element, x <= 0.1, tag="left")
select(mesh, :element, x >= 0.1, tag="right")
add_contact_elements(mesh, "left", "right", tag="interface")
select(mesh, :element, "left", :node, tag="left")
select(mesh, :element, "right", :node, tag="right")

mapper = RegionMapper()
add_mapping(mapper, :bulk, MechSolid, LinearElastic, E=E, nu=nu)
add_mapping(mapper, "interface", MechContact, CoulombContact, kn=kn, ks=ks, mu=mu)

model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=1.0)
ana = MechAnalysis(model)

select(model, :element, "interface", :ip, tag="jips")
log = add_logger(ana, :ip, "jips")
add_monitor(ana, :ip, "jips", (:w, :σn, :τ))

# # Opening (tension): pure frictional model should not sustain significant tensile traction.
# stage = add_stage(ana, nouts=10)
# add_bc(stage, :node, "left", ux=0, uy=0)
# add_bc(stage, :node, "right", ux=0.002)

# status = run(ana, autoinc=true, tol=0.02, rspan=0.02, quiet=false)

# @test status.successful
# w  = log.table["w"][end]
# σn = log.table["σn"][end]
# τ  = log.table["τ"][end]
# @test w > 0.0
# @test isapprox(σn, 0.0, atol=1e-8)
# @test isapprox(τ, 0.0, atol=1e-8)

# Closing (compression): normal traction should become compressive.

stage = add_stage(ana, nouts=10)
add_bc(stage, :node, "left", ux=0, uy=0)
add_bc(stage, :node, "right", ux=-0.001)

status = run(ana, autoinc=true, tol=0.01, rspan=0.02, quiet=true)

@test status.successful
w  = log.table["w"][end]
σn = log.table["σn"][end]
τ  = log.table["τ"][end]
@test w < 0.0
@test σn <= 1.0e-8

# # Shear

# stage = add_stage(ana, nouts=10)
# add_bc(stage, :node, "left", ux=0, uy=0)
# add_bc(stage, :node, "right", ux=0, uy=0.001)

# status = run(ana, autoinc=true, tol=0.02, rspan=0.02, quiet=true)

# @test status.successful
# w  = log.table["w"][end]
# σn = log.table["σn"][end]
# τ  = log.table["τ"][end]
# @test w < 0.0
# @test σn <= 1.0e-8
# @test τ >= 0.0

# Compression + Shear

stage = add_stage(ana, nouts=10)
add_bc(stage, :node, "left", ux=0, uy=0)
add_bc(stage, :node, "right", ux=-0.002, uy=0.005)

status = run(ana, autoinc=true, tol=0.02, rspan=0.02, quiet=true)

@test status.successful
w  = log.table["w"][end]
σn = log.table["σn"][end]
τ  = log.table["τ"][end]
@test w < 0.0
@test σn <= 1.0e-8
@test τ >= 0.0


chart = Chart(xlabel=t"$σ_n$", ylabel=t"$τ$")
add_line(chart, log.table["σn"], log.table["τ"], mark=:circle)
save(chart, "coulomb-contact.pdf")

# chart = Chart(xlabel=t"$w$", ylabel=t"$σ_n$")
# add_line(chart, log.table["w"]*1000, log.table["σn"], mark=:circle)
# save(chart, "coulomb-contact.pdf")
