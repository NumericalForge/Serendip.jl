using Serendip
using Test

# ❱❱❱ Geometry and mesh

geo = GeoModel()
bl1  = add_block(geo, [0, 0], [0.1, 0.1]; nx=1, ny=1, shape=QUAD4, tag="bulk")
bl2  = add_block(geo, [0.1, 0], [0.2, 0.1]; nx=1, ny=1, shape=QUAD4, tag="bulk")
mesh = Mesh(geo)
add_cohesive_elements(mesh, tag="joint")

left_elem = select(mesh, :element, :bulk, x<=0.1)
select(get_nodes(left_elem), tag="left")

right_elem = select(mesh, :element, :bulk, x>=0.1)
select(get_nodes(right_elem), tag="right")

E = 27.e6
nu = 0.2
ft = 2.4e3
fc = -24e3
zeta = 5.0
wc = 1.7e-4

# ❱❱❱ Finite element analysis using MohrCoulombCohesive

printstyled("\nMohrCoulombCohesive:\n\n", color=:yellow, bold=true)

mapper = RegionMapper()
add_mapping(mapper, "bulk", MechBulk, LinearElastic, E=E, nu=nu)
add_mapping(mapper, "joint", MechInterface, MohrCoulombCohesive, E=E, nu=nu, ft=ft, mu=1.4, zeta=zeta, wc=wc)

model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=1.0)
ana = MechAnalysis(model)

joint = select(model.elems, "joint")
select( get_ips(joint), tag="jips" )
log1 = add_logger(ana, :ip, "jips")
add_monitor(ana, :ip, "jips", (:σn, :τ))

stage = add_stage(ana, nincs=80, nouts=20)
add_bc(stage, :node, "left", ux=0, uy=0)
add_bc(stage, :node, "right", ux=1e-9, uy=8e-5)

run(ana, autoinc=true, maxits=3, tol=0.001, rspan=0.01, dTmax=0.1, tangent_scheme=:ralston)


# ❱❱❱ Finite element analysis using PowerYieldCohesive

printstyled("\nPowerYieldCohesive:\n\n", color=:yellow, bold=true)

mapper = RegionMapper()
add_mapping(mapper, "bulk", MechBulk, LinearElastic, E=E, nu=nu)
add_mapping(mapper, "joint", MechInterface, PowerYieldCohesive, E=E, nu=nu, fc=fc, ft=ft, zeta=zeta, wc=wc, alpha=1.5, gamma=0.05, theta=1.5)

model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=1.0)
ana = MechAnalysis(model)

joint = select(model.elems, "joint")
select( get_ips(joint), tag="jips" )
log1 = add_logger(ana, :ip, "jips")
add_monitor(ana, :ip, "jips", (:σn, :τ))

stage = add_stage(ana, nincs=80, nouts=20)
add_bc(stage, :node, "left", ux=0, uy=0)
add_bc(stage, :node, "right", ux=1e-9, uy=8e-5)

run(ana, autoinc=true, maxits=3, tol=0.001, rspan=0.01, dTmax=0.1, tangent_scheme=:ralston)


# ❱❱❱ Finite element analysis using AsinhYieldCohesive

printstyled("\nAsinhYieldCohesive:\n\n", color=:yellow, bold=true)

mapper = RegionMapper()
add_mapping(mapper, "bulk", MechBulk, LinearElastic, E=E, nu=nu)
add_mapping(mapper, "joint", MechInterface, AsinhYieldCohesive, E=E, nu=nu, fc=fc, ft=ft, zeta=zeta, wc=wc, alpha=0.6, gamma=0.05, theta=1.5)

model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=1.0)
ana = MechAnalysis(model)

joint = select(model.elems, "joint")
select( get_ips(joint), tag="jips" )
log1 = add_logger(ana, :ip, "jips")
add_monitor(ana, :ip, "jips", (:σn, :τ))

stage = add_stage(ana, nincs=80, nouts=20)
add_bc(stage, :node, "left", ux=0, uy=0)
add_bc(stage, :node, "right", ux=1e-9, uy=8e-5)

run(ana, autoinc=true, maxits=3, tol=0.001, rspan=0.01, dTmax=0.1, tangent_scheme=:ralston)