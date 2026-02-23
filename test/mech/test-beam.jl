using Serendip

geo = GeoModel(size=0.1)
# add_rectangle(geo, [0,0,0], 7, 0.7)
add_block(geo, [0,0,0], 7, 0.7, 0.0, nx=60, ny=30, shape=TRI3)

mesh = Mesh(geo)

mapper = RegionMapper()
add_mapping(mapper, :bulk, MechBulk, LinearElastic, E=25e6, nu=0.25)

model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=1.0)
ana   = MechAnalysis(model)

stage = add_stage(ana)

# add_monitor(ana, :node, 

add_bc(stage, :node, (x==2.8, y==0.7), fy=-10000)
add_bc(stage, :node, (x==0), ux=0, uy=0)
add_bc(stage, :face, (y==0.7), ty=-2000)

run(ana)

save(model, "test-beam.vtu")