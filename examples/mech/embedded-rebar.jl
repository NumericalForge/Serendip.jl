using Serendip

using Serendip

# Mesh generation
geo = GeoModel()
add_box(geo, [0.0, 0.0, 0.0], 0.5, 6.0, 0.5)
p1 = add_point(geo, [0.05, 0.05, 0.05])
p2 = add_point(geo, [0.05, 5.95, 0.05])
l1 = add_line(geo, p1, p2)
path = add_path(geo, [l1], tag="embedded", embedded=true)
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
add_mapping(mapper, "embedded", MechBar, VonMises, E=Es, A=A, fy=fy, H=H)

model = FEModel(mesh, mapper)
ana   = MechAnalysis(model)

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
save(plot, "bar-interface.pdf")

plot = DomainPlot(model,
    field = "σx´",
    field_factor = 1e-3,
    colormap = :rainbow,
    view_mode = :wireframe,
    label = L"\sigma_{x'}",
    warp = 20,
)
save(plot, "bar-interface.pdf")



# # Mesh generation
# bl  = Block( [0 0 0; 1.0 6.0 1.0], nx=1, ny=10, nz=3, tag="solids")
# bl1 = BlockInset( [0.2 0.2 0.2; 0.2 5.8 0.2], curvetype="polyline", embedded=true, tag="bars")
# bl2 = copy(bl1)
# move!(bl2, dx=0.6)
# bls = [ bl, bl1, bl2 ]

# msh = Mesh(bls)

# # FEM analysis
# mats = [
#         "solids" => MechBulk => LinearElastic => (E=1.e4, nu=0.25),
#         "bars"  => MechBar => VonMises => (E=1.e8, A=0.005, fy=500e3),
#        ]
# model = Model(msh, mats, ana)

# bcs = [
#     :(y==0 && z==0) => NodeBC(ux=0, uy=0, uz=0),
#     :(y==6 && z==0) => NodeBC(ux=0, uy=0, uz=0),
#     :(z==1) => SurfaceBC(tz=-1000),
# ]
# addstage!(model, bcs, nincs=20)

# solve!(model)
