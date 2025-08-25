using Serendip
using Test

geo = GeoModel()
bl1  = add_block(geo, [0, 0], [0.1, 0.1]; nx=1, ny=1, shape=QUAD4, tag="bulk")
bl2  = add_block(geo, [0.1, 0], [0.2, 0.1]; nx=1, ny=1, shape=QUAD4, tag="bulk")
mesh = Mesh(geo)
add_cohesive_elements(mesh, tag="joint")

left_elem = select(mesh, :element, :bulk, x<=0.1)
select(get_nodes(left_elem), tag="left")

right_elem = select(mesh, :element, :bulk, x>=0.1)
select(get_nodes(right_elem), tag="right")

# finite element analysis
E = 27.e6
ft = 2.4e3
fc = -24e3
alpha=1.5

mapper = RegionMapper()
add_mapping(mapper, "bulk", MechBulk, LinearElastic, E=E, nu=0.2)
# add_mapping(mapper, "joint", MechInterface, PowerYieldCrack, E=E, nu=0.2, fc=fc, ft=ft, zeta=5.0, wc=1.7e-4, alpha=alpha, gamma=0.05, theta=1.5)
add_mapping(mapper, "joint", MechInterface, LinearInterface, ks=1e8, kn=1e8)

# for i in (2,3)
# for i in 2

    ctx = Context(stress_state=:plane_stress, thickness=1.0)
    model = FEModel(mesh, mapper, ctx)
    ana = MechAnalysis(model)

    joint = select(model.elems, "joint")
    select( get_ips(joint), tag="jips" )
    log1 = add_logger(ana, :ip, "jips")

    stage = add_stage(ana, nincs=80, nouts=20)
    add_bc(stage, :node, "left", ux=0, uy=0)
    add_bc(stage, :node, "right", ux=-1e-9, uy=8e-5)

    run(ana, autoinc=true, maxits=3, tol=0.001, rspan=0.01, dTmax=0.1, scheme=:Ralston)

    # if makeplots
    #     table = log1.table

    #     chart = Chart(; xlabel=L"$u_p$", ylabel=L"\sigma")
    #     addplot!(chart, LinePlot(table[:jup], table[:jσn], marker=:circle))
    #     save(chart, "up-sn.pdf")

    #     chart = Chart(; xlabel=L"\sigma_n", ylabel=L"\tau")
    #     addplot!(chart, LinePlot(table[:jσn], table[:js2], marker=:circle))
    #     save(chart, "sn-tau.pdf")

    # end


# end
