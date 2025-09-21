using Serendip

geo = GeoModel()
add_block(geo, [0.0, 0.0], [2.0, 1.0]; n=2, shape=LIN2, tag="frame")

mesh = Mesh(geo, ndim=2)

mapper = RegionMapper()
add_mapping(mapper, "frame", MechFrame, LinearElastic, E=10.0, A=1.0, I=1.0)

model = FEModel(mesh, mapper)
ana = MechAnalysis(model)

stage = add_stage(ana)
add_bc(stage, :node, (x==0.0), ux=0.0, uy=0.0)
add_bc(stage, :node, (x==2.0), ux=0.0, uy=0.0)
add_bc(stage, :node, (x==1.0), mz=1.0)
add_bc(stage, :node, (x==1.0), fy=-1.0)
add_bc(stage, :edge, (x>=1.0), qy=-12)

@test run(ana).successful