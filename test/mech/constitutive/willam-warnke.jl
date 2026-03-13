using Serendip
# mesh

geo = GeoModel()
add_block(geo, [0, 0, 0], 0.1, 0.1, 0.2, nx=2, ny=2, nz=2, tag="solids")
mesh = Mesh(geo)

# elements and constitutive model

mapper = RegionMapper()
add_mapping(
    mapper,
    "solids",
    MechSolid,
    WillamWarnke,
    E=30e6,
    nu=0.25,
    fc=-30e3,
    epsc=-0.002,
    ft=3e3,
    GF=0.1,
    ft_law=:hordijk,
    fc_law=:popovics,
    beta=1.1,
)

model = FEModel(mesh, mapper)
ana = MechAnalysis(model)

log = add_logger(ana, :ip, [0.05, 0.05, 0.1], "cscp.dat")

# boundary conditions

stage = add_stage(ana, nincs=10, nouts=10)
add_bc(stage, :node, x==0, ux=0)
add_bc(stage, :node, y==0, uy=0)
add_bc(stage, :node, z==0, uz=0)
add_bc(stage, :node, z==0.2, uz=-0.001)
add_bc(stage, :node, x==0.1, ux=-0.0005)

run(ana, tol=1e-1, autoinc=true, quiet=false).successful

# Plotting

table = log.table

chart = Chart(
    xlabel="ezz x 1000",
    ylabel="szz",
    xmult=1e3,
)
add_series(chart, table["εzz"], table["σzz"])
save(chart, "chart.pdf")
