#❱❱❱ Cracking in a 3D cantilever beam ❰❰❰#

using Serendip

ℓ = 0.4   # x direction
b = 0.025 # y direction
h = 0.1   # z direction

geo = GeoModel(size=0.01)
box = add_box(geo, [0,0, 0], ℓ, b, h; tag="bulk")

p1 = add_point(geo, [ℓ/2, 0, h])
p2 = add_point(geo, [ℓ/2, b, h])
l1 = add_line(geo, p1, p2)
fragment(geo, l1, box)

p1 = add_point(geo, [0.00*ℓ, 0.25*b, 0.95*h])
p2 = add_point(geo, [0.95*ℓ, 0.25*b, 0.95*h])
bar = add_line(geo, p1, p2)
p = add_path(geo, [bar], tag="bar", interface_tag="bar-interface") 
add_array(geo, p, ny=3, dy=0.25*b)

mesh = Mesh(geo)

select(mesh.elems, :bulk, tag="bulk")
add_cohesive_elements(mesh, tag="cohesive")

# finite element analysis
E     = 37.e6
ft    = 2.4e3
fc    = -24e3
alpha = 1.5
GF    = 0.07

mapper = RegionMapper()
add_mapping(mapper, "bulk", MechBulk, LinearElastic, E=E, nu=0.2)
add_mapping(mapper, "cohesive", MechCohesive, MohrCoulombCohesive, E=E, nu=0.2, ft=ft, GF=GF, mu=1.5)
add_mapping(mapper, "bar", MechBar, LinearElastic, E=2e5, A=0.005)
add_mapping(mapper, "bar-interface", MechBondSlip, LinearBondSlip, ks=1e10, kn=1e9, p=0.01) 

model = FEModel(mesh, mapper)
ana   = MechAnalysis(model, outkey="crack-3d", outdir="crack-3d")

log1 = add_logger(ana, :nodalreduce, x==ℓ)
add_monitor(ana, :node, (x==ℓ, y==0, z==0), :uz)

stage = add_stage(ana, nincs=10, nouts=50)
add_bc(stage, :node, x==0, ux=0, uy=0, uz=0)
add_bc(stage, :node, x==ℓ, uz=-0.0002)

run(ana, autoinc=true, maxits=3, tol=0.5, rspan=0.03, tangent_scheme=:ralston, quiet=false)

# ❱❱❱ Post-processing

plot = DomainPlot(model,
    field      = "σxx",
    field_mult = 1e-3,
    warp       = 200,
    azimuth    = -80,
    colormap   = :spectral,
    diverging  = true,
    edge_color = :gray,
    edge_width = 0.1,
    colorbar   = :bottom,
    label      = L"$\sigma_{xx}$ [MPa]",
)
save(plot, "beam-crack-3d.pdf")

chart = Chart(
    xlabel = L"$u_z$ [mm]",
    ylabel = L"$σ_{zz}$ [kN]",
)
add_line(chart, -log1.table[:uz]*1e3, -log1.table[:fz], mark=:circle)
save(chart, "beam-crack-3d-chart.pdf")