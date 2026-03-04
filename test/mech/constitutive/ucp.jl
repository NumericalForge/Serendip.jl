using Serendip

# ❱❱❱ Geometry and mesh

ℓ = 0.1
h = 0.2

geo = GeoModel()

add_block(geo, [0,0,0], ℓ,ℓ,h, shape=HEX8, nx=2, ny=2, nz=2)


mesh= Mesh(geo, quadratic=true)
select(mesh, :element, :bulk, tag="solids")

# ❱❱❱ Elements and constitutive models

fc = -30.00e3
ft = 3.0e3

mapper = RegionMapper()
add_mapping(mapper, "solids", MechBulk, UCP, 
    E=30.0e6, nu=0.25, alpha=0.66, beta=1.15, 
    fc=fc, epsc=-0.002, eta=2.2, ft=ft, wc=0.0001 )

# ❱❱❱ Analysis 1: compression

# printstyled("\nCompression test\n", color=:yellow, bold=true)

model = FEModel(mesh, mapper)
ana = MechAnalysis(model, outkey="compression")

log_c = add_logger(ana, :ip, [ℓ/2, ℓ/2, h/2], "ip-compression.dat")
add_monitor(ana, :ip, [ℓ/2, ℓ/2, h/2], :σzz)

stage1 = add_stage(ana)
uu = -0.0015
add_bc(stage1, :node, x==0, ux=0)
add_bc(stage1, :node, y==0, uy=0)
add_bc(stage1, :node, z==0, uz=0)
# add_bc(stage1, :node, x==ℓ, ux=uu)
# add_bc(stage1, :node, y==ℓ, uy=uu)
add_bc(stage1, :node, z==h, uz=uu)

run(ana, autoinc=true, tol=0.1, rspan=0.02, dT0=0.002, quiet=false)

# ❱❱❱ Analysis 2: tension

printstyled("\nTension test\n", color=:yellow, bold=true)

model = FEModel(mesh, mapper)
ana = MechAnalysis(model, outkey="tension")

log_t = add_logger(ana, :ip, [ℓ/2, ℓ/2, h/2], "ip-tension.dat")
add_monitor(ana, :ip, [ℓ/2, ℓ/2, h/2], :σzz)

stage1 = add_stage(ana, nouts=10)
uu = 0.0002
add_bc(stage1, :node, x==0, ux=0)
add_bc(stage1, :node, y==0, uy=0)
add_bc(stage1, :node, z==0, uz=0)
# add_bc(stage1, :node, x==ℓ, ux=uu)
add_bc(stage1, :node, z==h, uz=uu)

run(ana, autoinc=true, tol=0.1, quiet=false)


# ❱❱❱ Post-process

table_c = log_c.table
table_t = log_t.table

chart = Chart(
    xlabel=t"$ε_(z z)$",
    ylabel=t"$σ_(z z)$",
)
add_line(chart, table_c["εzz"], table_c["σzz"], mark=:circle)
add_line(chart, table_t["εzz"], table_t["σzz"], mark=:utriangle, mark_size=1)
save(chart, "σ-ε.pdf")

chart = Chart(
    xlabel=t"$ξ$",
    ylabel=t"$ρ$",
)
add_line(chart, table_c["ξ"], table_c["ρ"], mark=:circle)
add_line(chart, table_t["ξ"], table_t["ρ"], mark=:utriangle)
save(chart, "ξ-ρ.pdf")
