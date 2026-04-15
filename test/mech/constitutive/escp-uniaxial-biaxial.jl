using Serendip
using Test

# ❱❱❱ Geometry and mesh

ℓ = 0.1
h = 0.2

geo = GeoModel()
add_block(geo, [0, 0, 0], ℓ, ℓ, h, shape=:hex8, nx=4, ny=4, nz=4)
mesh = Mesh(geo, quadratic=true)
select(mesh, :element, :bulk, tag="solids")

# ❱❱❱ Elements and constitutive models

fc = -30.0e3
ft = 3.0e3

mapper = RegionMapper()
add_mapping(
    mapper,
    "solids",
    MechSolid,
    ESCP,
    E     = 30.0e6,
    nu    = 0.25,
    beta  = 1.15,
    fc    = fc,
    epsc  = -0.002,
    ft    = ft,
    wc    = 0.0001
)

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

status = run(ana, autoinc=true, tol=0.01, dT0=0.001, dTmax=0.005, quiet=false)
@test status.successful
@test maximum(log_t.table["σzz"]) > 0.5*ft
@test log_t.table["εtp"][end] > 0.0


# ❱❱❱ Post-process

table_c = log_c.table
table_b = log_b.table
table_t = log_t.table
table_c_top = log_c_top.table
table_b_xl = log_b_xl.table
table_b_yl = log_b_yl.table
table_t_top = log_t_top.table

chart11 = Chart(
    title="Uniaxial vs Biaxial Compression (IP)",
    xlabel="`ε`",
    ylabel="`σ`",
    legend=:top_left,
)
add_line(chart11, table_c["εzz"], table_c["σzz"], mark=:circle, label="uniaxial compression")
add_line(chart11, table_b["εxx"], table_b["σxx"], mark=:square, label="biaxial compression xx")
add_line(chart11, table_b["εyy"], table_b["σyy"], mark=:diamond, label="biaxial compression yy")

chart12 = Chart(
    title="Uniaxial vs Biaxial Compression (face)",
    xlabel="`ε`",
    ylabel="`σ`",
    legend=:top_left,
)
add_line(chart12, table_c_top["uz"]/h, table_c_top["fz"]/ℓ^2, mark=:circle, label="uniaxial compression")
add_line(chart12, table_b_xl["ux"]/ℓ, table_b_xl["fx"]/(h*ℓ), mark=:square, label="biaxial compression xx")
add_line(chart12, table_b_yl["uy"]/ℓ, table_b_yl["fy"]/(h*ℓ), mark=:diamond, label="biaxial compression yy")

chart21 = Chart(
    title="Uniaxial Tension (IP)",
    xlabel="`ε`",
    ylabel="`σ`",
    legend=:right,
)
add_line(chart21, table_t["εzz"], table_t["σzz"], mark=:utriangle, mark_size=1, label="uniaxial tension")

chart22 = Chart(
    title="Uniaxial Tension (face)",
    xlabel="`ε`",
    ylabel="`σ`",
    legend=:right,
)
add_line(chart22, table_t_top["uz"]/h, table_t_top["fz"]/ℓ^2, mark=:utriangle, mark_size=1, label="uniaxial tension")

chart31 = Chart(
    title="Trajectories in `ξ-ρ` space",
    xlabel="`ξ`",
    ylabel="`ρ`",
    legend=:bottom_left,
)
add_line(chart31, table_c["ξ"], table_c["ρ"], mark=:circle, label="uniaxial compression")
add_line(chart31, table_b["ξ"], table_b["ρ"], mark=:square, label="biaxial compression")
add_line(chart31, table_t["ξ"], table_t["ρ"], mark=:utriangle, label="uniaxial tension")

grid = ChartGrid(
    size = (14cm, 18cm),
    background =:old_paper, 
    quiet=true,
)

add_chart(grid, chart11, (1, 1))
add_chart(grid, chart12, (1, 2))
add_chart(grid, chart21, (2, 1))
add_chart(grid, chart22, (2, 2))
add_chart(grid, chart31, (3, 1))

save(grid, "escp-uniaxial-biaxial.pdf")
