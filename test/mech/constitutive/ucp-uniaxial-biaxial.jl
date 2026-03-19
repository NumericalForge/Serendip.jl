using Serendip
using Test

# ❱❱❱ Geometry and mesh

ℓ = 0.1
h = 0.2

geo = GeoModel()
add_block(geo, [0, 0, 0], ℓ, ℓ, h, shape=HEX8, nx=2, ny=2, nz=2)
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
    UCP,
    E     = 30.0e6,
    nu    = 0.25,
    alpha = 0.66,
    beta  = 1.15,
    fc    = fc,
    epsc  = -0.002,
    eta   = 2.2,
    ft    = ft,
    wc    = 0.0001
)

# Uniaxial compression

printstyled("\nUniaxial compression:\n\n", color=:yellow, bold=true)
model = FEModel(mesh, mapper)
ana = MechAnalysis(model)
log_c = add_logger(ana, :ip, [ℓ/2, ℓ/2, h/2])

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

stage = add_stage(ana, nouts=10)
add_bc(stage, :node, x == 0, ux=0)
add_bc(stage, :node, y == 0, uy=0)
add_bc(stage, :node, z == 0, uz=0)
add_bc(stage, :node, z == h, uz=0.00035)

status = run(ana, autoinc=true, tol=0.01, dT0=0.001, quiet=false)
@test status.successful
@test maximum(log_t.table["σzz"]) > 0.5*ft
@test log_t.table["εtp"][end] > 0.0


# ❱❱❱ Post-process

table_c = log_c.table
table_b = log_b.table
table_t = log_t.table

chart = Chart(
    xlabel=t"$ε$",
    ylabel=t"$σ$",
    legend=:top_left,
)
add_line(chart, table_c["εzz"], table_c["σzz"], mark=:circle, label="uniaxial compression")
add_line(chart, table_b["εxx"], table_b["σxx"], mark=:square, label="biaxial compression xx")
add_line(chart, table_b["εyy"], table_b["σyy"], mark=:diamond, label="biaxial compression yy")
add_line(chart, table_t["εzz"], table_t["σzz"], mark=:utriangle, mark_size=1, label="uniaxial tension")
save(chart, "ucp-σ-ε.pdf")

chart = Chart(
    xlabel=t"$ξ$",
    ylabel=t"$ρ$",
    legend=:top_left,
)
add_line(chart, table_c["ξ"], table_c["ρ"], mark=:circle, label="uniaxial compression")
add_line(chart, table_b["ξ"], table_b["ρ"], mark=:square, label="biaxial compression")
add_line(chart, table_t["ξ"], table_t["ρ"], mark=:utriangle, label="uniaxial tension")
save(chart, "ucp-ξ-ρ.pdf")
