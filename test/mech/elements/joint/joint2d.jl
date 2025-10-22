using Serendip
using Test

geo = GeoModel(size=0.01)
# bl1  = add_block(geo, [0, 0], [0.1, 0.1]; nx=1, ny=1, shape=QUAD4, tag="bulk")
# bl2  = add_block(geo, [0.1, 0], [0.2, 0.1]; nx=1, ny=1, shape=QUAD4, tag="bulk")
# add_block(geo, [0, 0], [0.2, 0.1]; nx=4, ny=4, shape=QUAD4, tag="bulk")
# add_block(geo, [0, 0], [0.2, 0.1]; nx=2, ny=1, shape=QUAD4, tag="bulk")
add_rectangle(geo, [0,0, 0], 0.2, 0.1; tag="bulk")

mesh = Mesh(geo)
# save(model, "mesh.vtk")


select(mesh.elems, tag="bulk")
Serendip.add_cohesive_elements2(mesh, tag="cohesive")
# Serendip.add_cohesive_elements2(mesh, :(x<0.05), tag="cohesive")

# left_elem = select(mesh, :element, :bulk, x<=0.1)
# select(get_nodes(left_elem), tag="left")

# right_elem = select(mesh, :element, :bulk, x>=0.1)
# select(get_nodes(right_elem), tag="right")

# finite element analysis
E     = 27.e6
ft    = 2.4e3
fc    = -24e3
alpha = 1.5
GF    = 0.1

mapper = RegionMapper()
add_mapping(mapper, "bulk", MechBulk, LinearElastic, E=E, nu=0.2)
# add_mapping(mapper, "cohesive", MechCohesive, PowerYieldCohesive, E=E, nu=0.2, fc=fc, ft=ft, zeta=5.0, wc=1.7e-4, alpha=alpha, gamma=0.05, theta=1.5)
add_mapping(mapper, "cohesive", MechCohesive, MohrCoulombCohesive, E=E, nu=0.2, ft=ft, GF=GF, mu=1.5)
# # add_mapping(mapper, "cohesive", MechContact, LinearInterface, ks=1e8, kn=1e8)


model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=1.0)
ana   = MechAnalysis(model)


# # for i in (2,3)
# # for i in 2

joint = select(model.elems, "cohesive")
select( get_ips(joint), tag="jips" )
log1 = add_logger(ana, :ip, "jips")

stage = add_stage(ana, nincs=10, nouts=20)
add_bc(stage, :node, x==0, ux=0, uy=0)
# add_bc(stage, :node, x==0.2, ux=0.4*1e-4, uy=0)
add_bc(stage, :node, x==0.2, uy=-0.8*1e-4)

# run(ana, autoinc=true, maxits=3, tol=0.1, rspan=0.1, dTmax=0.1, quiet=false)
run(ana, autoinc=true, maxits=3, tol=0.1, rspan=0.01, dTmax=0.1, tangent_scheme=:ralston, quiet=false)

# save(model, "mesh.vtk")

#     # if makeplots
#     #     table = log1.table

#     #     chart = Chart(; xlabel=L"$u_p$", ylabel=L"\sigma")
#     #     addplot!(chart, LinePlot(table[:jup], table[:jσn], mark=:circle))
#     #     save(chart, "up-sn.pdf")

#     #     chart = Chart(; xlabel=L"\sigma_n", ylabel=L"\tau")
#     #     addplot!(chart, LinePlot(table[:jσn], table[:js2], mark=:circle))
#     #     save(chart, "sn-tau.pdf")

#     # end


# # end
