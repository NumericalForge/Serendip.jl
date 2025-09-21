using Serendip

# ❱❱❱ Geometry and mesh

h = 0.2

geo = GeoModel()

add_block(geo, [0,0,0], [0.1,0.1,h], shape=HEX8, nx=1, ny=1, nz=1)
# box = add_box(geo, [0,0,0], 0.1, 0.1, h)
# set_size(geo, box, 0.05)


mesh= Mesh(geo, quadratic=true)
select(mesh, :element, :bulk, tag="solids")

# ❱❱❱ Elements and constitutive models

fc = -30.87e3
ft = 2.95e3

mapper = RegionMapper()
add_mapping(mapper, "solids", MechBulk, UCP, 
    E=30e6, nu=0.2, 
    alpha=0.666, beta=1.15, 
    fc=fc, epsc=-0.002, eta=2.2,
    ft=ft, 
    GF=0.1,
    )

# ❱❱❱ Analysis 1: compression

# model = FEModel(mesh, mapper)
# ana = MechAnalysis(model, outkey="compression")

# log = add_logger(ana, :ip, [0.05, 0.05, 0.05])
# add_monitor(ana, :ip, [0.05, 0.05, 0.05], :σzz)

# stage1 = add_stage(ana, nouts=10)
# add_bc(stage1, :node, x==0, ux=0)
# add_bc(stage1, :node, y==0, uy=0)
# add_bc(stage1, :node, z==0, uz=0)
# # add_bc(stage1, :node, z==0, ux=0, uy=0, uz=0)
# add_bc(stage1, :node, z==0.1, uz=-0.001)

# run(ana, autoinc=true, tol=0.01)
# # run(ana, autoinc=true, tangent_tangent_scheme=:ralston)

# # Post-process
# table = log.table
# chart = Chart(
#     xlabel=L"$ε_{zz} × 10^3$", xmult=1000,
#     ylabel=L"$σ_{zz}$ [MPa]", ymult=1e-3,
# )
# add_series(chart, -table["εzz"], -table["σzz"], mark=:circle)
# save(chart, "compression.pdf")


# ❱❱❱ Analysis 2: tension

model = FEModel(mesh, mapper)
ana = MechAnalysis(model, outkey="tension")

log = add_logger(ana, :ip, [0.05, 0.05, h/2])
add_monitor(ana, :ip, [0.05, 0.05, h/2], :σzz)

stage1 = add_stage(ana, nouts=10)
add_bc(stage1, :node, x==0, ux=0)
add_bc(stage1, :node, y==0, uy=0)
add_bc(stage1, :node, z==0, uz=0)
# add_bc(stage1, :node, z==0, ux=0, uy=0, uz=0)
add_bc(stage1, :node, z==h, uz=0.0002)

run(ana, autoinc=true, dTmax=0.01)

# Post-process
table = log.table
chart = Chart(
    xlabel=L"$ε_{zz} × 10^3$", xmult=1000,
    ylabel=L"$σ_{zz}$ [MPa]", ymult=1e-3,
)
add_series(chart, table["εzz"], table["σzz"], mark=:circle)
save(chart, "tension.pdf")

