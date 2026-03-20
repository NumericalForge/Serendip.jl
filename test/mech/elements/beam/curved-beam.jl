using Serendip

geo = GeoModel()
p1  = add_point(geo, [0.0, 0.0, 0.0])
p2  = add_point(geo, [0.5, 0.0, 0.5])
p3  = add_point(geo, [1.0, 0.0, 0.0])
bz  = add_bezier(geo, [p1, p2, p3])
set_transfinite_curve(geo, bz, 3)

mesh = Mesh(geo, ndim=3, quadratic=true)

mapper = RegionModel(MechBeam, LinearElastic, E=200e6, nu=0.2, b=0.1, h=0.1)

model = FEModel(mesh, mapper)
ana   = MechAnalysis(model)

log1 = add_logger(ana, :node, x==1)
add_monitor(ana, :node, x==0, :(my, fx, fz))
add_monitor(ana, :node, x==1, :uz)

stage = add_stage(ana)
add_bc(stage, :node, x==0, rx=0, ry=0, rz=0, ux=0, uy=0, uz=0)
add_bc(stage, :node, x==1, fx=1, fz=2)

run(ana, quiet=false)

@test log1.table[:uz][end] ≈ 0.00052768 atol=1e-5