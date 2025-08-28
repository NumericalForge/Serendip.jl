using Serendip

geo = GeoModel()
add_block(geo, [0.0, 0.0, 0.0], [0.2, 2.0, 0.2], nx=2, ny=12, nz=2, shape=HEX8, tag="solids")
mesh = Mesh(geo)

mapper = RegionMapper()
add_mapping(mapper, "solids", MechBulk, LinearElastic, E=36e6, nu=0.2, rho=24.0)
model = FEModel(mesh, mapper)
ana = DynamicAnalysis(model)

add_logger(ana, :node, (x==0.2, y==2, z==0.2), "end-node.dat")

stage = add_stage(ana, nincs=10)
add_bc(stage, :node, (y==0, z==0), ux=0, uy=0, uz=0)
add_bc(stage, :node, (y==2, z==0), uz=0)
add_bc(stage, :node, (y==1, z==0.2), fz=-10)

run(ana)
