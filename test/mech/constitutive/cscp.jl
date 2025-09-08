using Serendip
using Test

geo = GeoModel()

add_block(geo, [0,0,0], [0.1,0.1,0.1], shape=HEX8, nx=1, ny=1, nz=1, tag="solids")

mesh= Mesh(geo)

fc = -30.87e3
ft = 2.95e3


mapper = RegionMapper()
add_mapping(mapper, "solids", MechBulk, UCP, 
    E=30e6, nu=0.2, 
    alpha=0.666, beta=1.15, 
    ft=ft, 
    fc=fc, 
    epsc=-0.002, 
    GF=0.1,
    wc=0.00005
    )

model = FEModel(mesh, mapper)
ana = MechAnalysis(model)

add_logger(ana, :ip, [0.05, 0.05, 0.05], "cscp.dat")
# add_monitor(ana, :ip, [0.05, 0.05, 0.05], :(sxx, syy), stop=:( rho<0.3*rho_max))

stage = add_stage(ana)

add_bc(stage, :node, x==0, ux=0)
add_bc(stage, :node, y==0, uy=0)
add_bc(stage, :node, z==0, uz=0)
add_bc(stage, :node, z==0.1, ux=-0.0014)

run(ana, autoinc=true)
