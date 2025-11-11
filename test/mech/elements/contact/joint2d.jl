using Serendip
using Test

ℓ = 0.4
h = 0.1

geo = GeoModel(size=0.005)
# geo = GeoModel(size=0.015)
# bl1  = add_block(geo, [0, 0], [0.1, h]; nx=1, ny=1, shape=QUAD4, tag="bulk")
# bl2  = add_block(geo, [0.1, 0], [ℓ, h]; nx=1, ny=1, shape=QUAD4, tag="bulk")
# add_block(geo, [0, 0], [ℓ, h]; nx=4, ny=4, shape=QUAD4, tag="bulk")
# add_block(geo, [0, 0], [ℓ, h]; nx=3, ny=3, shape=QUAD4, tag="bulk")
# add_block(geo, [0, 0], [ℓ, h]; nx=2, ny=1, shape=QUAD4, tag="bulk")
add_rectangle(geo, [0,0, 0], ℓ, h; tag="bulk")

p1 = add_point(geo, [0.00*ℓ, 0.95*h, 0])
p2 = add_point(geo, [0.95*ℓ, 0.95*h, 0])
bar = add_line(geo, p1, p2)
add_path(geo, [bar], tag="bar", interface_tag="bar-interface") 

mesh = Mesh(geo)
# save(model, "mesh.vtk")


select(mesh.elems, :bulk, tag="bulk")
Serendip.add_cohesive_elements(mesh, x<0.3, tag="cohesive")
# Serendip.add_cohesive_elements2(mesh, :(x>0.1), tag="cohesive")


# finite element analysis
E     = 37.e6
ft    = 2.4e3
fc    = -24e3
alpha = 1.5
GF    = 0.07

mapper = RegionMapper()
add_mapping(mapper, "bulk", MechBulk, LinearElastic, E=E, nu=0.2)
# add_mapping(mapper, "cohesive", MechCohesive, PowerYieldCohesive, E=E, nu=0.2, fc=fc, ft=ft, zeta=5.0, GF=GF, alpha=alpha, gamma=0.05, theta=1.5)
add_mapping(mapper, "cohesive", MechCohesive, MohrCoulombCohesive, E=E, nu=0.2, ft=ft, GF=GF, mu=1.5, zeta=10)
# # add_mapping(mapper, "cohesive", MechContact, LinearInterface, ks=1e8, kn=1e8)
add_mapping(mapper, "bar", MechBar, LinearElastic, E=2e5, A = 0.005)
add_mapping(mapper, "bar-interface", MechBondSlip, LinearBondSlip, ks=1e6, kn=1e9, p=0.01) 

model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=1.0)
ana   = MechAnalysis(model)

joint = select(model.elems, "cohesive")
select( get_ips(joint), tag="jips" )
log1 = add_logger(ana, :ip, "jips")

stage = add_stage(ana, nincs=10, nouts=50)
add_bc(stage, :node, x==0, ux=0, uy=0)
# add_bc(stage, :node, x==ℓ, ux=0.4*1e-4, uy=0)
# add_bc(stage, :node, x==ℓ, uy=-4*1e-4)
add_bc(stage, :node, x==ℓ, uy=-1.5*1e-4)

# run(ana, autoinc=true, maxits=3, tol=0.1, rspan=0.1, dTmax=0.1, quiet=false)
run(ana, autoinc=true, maxits=3, tol=0.5, rspan=0.01, dTmax=0.1, tangent_scheme=:ralston, quiet=false)

# save(model, "mesh.vtk")

#     # if makeplots
#     #     table = log1.table

#     #     chart = Chart(; xlabel=L"$u_p$", ylabel=L"\sigma")
#     #     addplot!(chart, LinePlot(table[:jup], table[:jσn], mark=:circle))
#     #     save(chart, "up-sn.pdf")

#     #     chart = Chart(; xlabel=L"\sigma_n", ylabel=L"\tau")
#     #     addplot!(chart, LinePlot(table[:jσn], table[:js2], mark=:circle))
#     #     save(chart, "sn-tau.pdf")

