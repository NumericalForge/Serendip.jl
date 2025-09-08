using Serendip

# Mesh generation
geo = GeoModel()
add_box(geo, [0.0, 0.0, 0.0], 0.5, 6.0, 0.5)
p1 = add_point(geo, [0.05, 0.05, 0.05])
p2 = add_point(geo, [0.05, 5.95, 0.05])
l1 = add_line(geo, p1, p2)
path = add_path(geo, [l1], tag="bars", interface_tag="interface")
add_array(geo, path, nx=3, dx=0.2)

mesh = Mesh(geo)
select(mesh, :element, :bulk, tag="solids")

# Fem analysis
E  = 24e6 # kPa
Es = 210e6 # kPa
fy = 240e3 # kPa
φ  = 0.01
p  = π*φ
A  = π*φ^2/4
H  = 0.0
nu = 0.3

mapper = RegionMapper()
add_mapping(mapper, "solids", MechBulk, LinearElastic, E=E, nu=nu)
add_mapping(mapper, "bars", MechBar, VonMises, E=Es, A=A, fy=fy, H=H)
add_mapping(mapper, "interface", MechBondSlip, LinearBondSlip, ks=1e7, kn=1e9, p=p)

model = FEModel(mesh, mapper)
ana = MechAnalysis(model)

stage = add_stage(ana, nouts=5)
add_bc(stage, :node, (y==0, z==0), ux=0, uy=0, uz=0)
add_bc(stage, :node, (y==6, z==0), ux=0, uz=0)
add_bc(stage, :face, (z==0.5), tz=-0.001)  # uniform load

run(ana, autoinc=true)

plot = DomainPlot(model,
    field = "σx´",
    field_factor = 1e-3,
    colormap = :rainbow,
    view_mode = :wireframe,
    label = L"\sigma_{x'}",
    warp = 20,
)
save(plot, "embedded-bar.pdf")
