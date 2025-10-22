using Serendip

# ❱❱❱ Geometry and mesh

ℓ = 0.1
h = 0.2

geo = GeoModel()

add_block(geo, [0,0,0], [ℓ,ℓ,h], shape=HEX8, nx=2, ny=2, nz=4)
# box = add_box(geo, [0,0,0], ℓ, ℓ, h)
# set_size(geo, box, 0.05)


mesh= Mesh(geo, quadratic=true)
select(mesh, :element, :bulk, tag="solids")

# ❱❱❱ Elements and constitutive models

fc = -30.00e3
ft = 3.0e3

mapper = RegionMapper()
add_mapping(mapper, "solids", MechBulk, UCP, 
    E=30.0e6, nu=0.25, alpha=0.66, beta=1.15, 
    fc=fc, epsc=-0.002, eta=2.2, ft=ft, wc=0.0001 )
# add_mapping(mapper, "solids", MechBulk, LinearElastic, E=30.0e6, nu=0.25)

# "solids" => MechSolid => CSCP => (E=30e6, nu=0.25, alpha=0.66, beta=1.15, fc=fc, ft=ft, epsc=-0.002, GF=0.1, wc=0.00005, n=2.0)


# ❱❱❱ Analysis 1: compression

# model = FEModel(mesh, mapper)
# ana = MechAnalysis(model, outkey="compression")

# log = add_logger(ana, :ip, [0.0, 0.0, 0.0])
# add_monitor(ana, :ip, [0.0, 0.0, 0.0], :σzz)

# stage1 = add_stage(ana)
# uu = -0.002
# add_bc(stage1, :node, x==0, ux=0)
# add_bc(stage1, :node, y==0, uy=0)
# add_bc(stage1, :node, z==0, uz=0)
# # add_bc(stage1, :node, x==ℓ, ux=uu)
# # add_bc(stage1, :node, y==ℓ, uy=uu)
# add_bc(stage1, :node, z==h, uz=uu)

# run(ana, autoinc=true, tol=0.006, rspan=0.02, dT0=0.002,)

# ❱❱❱ Analysis 2: tension

model = FEModel(mesh, mapper)
ana = MechAnalysis(model, outkey="tension")

log = add_logger(ana, :ip, [0.0, 0.0, 0.0], "ip.dat")
add_monitor(ana, :ip, [0.0, 0.0, 0.0], :σxx)

stage1 = add_stage(ana, nouts=10)
uu = 0.0002
add_bc(stage1, :node, x==0, ux=0)
add_bc(stage1, :node, y==0, uy=0)
add_bc(stage1, :node, z==0, uz=0)
# add_bc(stage1, :node, x==ℓ, ux=uu)
add_bc(stage1, :node, z==h, uz=uu)

run(ana, autoinc=true, tol=0.1)

# ❱❱❱ Post-process

table = log.table

chart = Chart(
    xlabel=L"$ε_{zz}$",
    ylabel=L"$σ_{zz}$",
)
add_series(chart, table["εzz"], table["σzz"], mark=:circle)
save(chart, "σ-ε.pdf")

chart = Chart(
    xlabel=L"$ξ$",
    ylabel=L"$ρ",
)
add_series(chart, table["ξ"], table["ρ"], mark=:circle)
save(chart, "ξ-ρ.pdf")

