using Serendip

# Geometry and generation >>>

geo = GeoModel()
add_box(geo, [0.0, 0.0, 0.0], 0.3, 6.0, 0.5)
p1 = add_point(geo, [0.05, 0.05, 0.05])
p2 = add_point(geo, [0.05, 5.95, 0.05])
l1 = add_line(geo, p1, p2)
path1 = add_path(geo, [l1], tag="bars", interface_tag="interface")
add_array(geo, path1, nx=3, dx=0.1)

p3 = add_point(geo, [0.05, 0.05, 0.45])
p4 = add_point(geo, [0.05, 5.95, 0.45])
l2 = add_line(geo, p3, p4)
path2 = add_path(geo, [l2], tag="bars", interface_tag="interface")
add_array(geo, path2, nx=2, dx=0.2)

mesh = Mesh(geo)
select(mesh, :element, :bulk, tag="solids")

# Finite element modeling >>>

E  = 24e6 # kPa
Es = 210e6 # kPa
fy = 240e3 # kPa
φ  = 0.01 # m
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
add_bc(stage, :node, (y==0), ux=0, uy=0, uz=0) # clamp
add_bc(stage, :node, (y==6, z==0), ux=0, uz=0) # roller
add_bc(stage, :face, (z==0.5), tz=-0.001)  # uniform load

run(ana, autoinc=true)

# ❱❱❱ Post-processing

plot = DomainPlot(model,
    field = "σx´", # stress in the bars
    field_factor = 1e-3,
    colormap = :rainbow,
    view_mode = :wireframe,
    label = L"\sigma_{x'}",
    warp = 20,
)
save(plot, "reinforced-beam.pdf")
