using Serendip

# # 2D beam

# coord = [ 0 0; 1 0; 0.5 0.5]
# n     = 2
# bl    = Block(coord, nx=n, cellshape = LIN3, tag="beam")
# mesh   = Mesh(bl, ndim=2)

# mats  = [ "beam" => MechBeam => ElasticBeam => (E=1e4, nu=0, thy=0.1, thz=0.1) ]

# ana = MechAnalysis()
# model = FEModel(mesh, mats, ana)

# monitors = [
#     :(x==0) => NodeMonitor(:(mz))
#     :(x==0) => NodeMonitor(:(fx))
#     :(x==0) => NodeMonitor(:(fy))
#     :(x==1) => NodeMonitor(:(uy))
# ]
# setmonitors!(model, monitors)

# bcs =
# [
#    :(x==0) => NodeBC(rz=0, ux=0, uy=0),
#    :(x==1) => NodeBC(fx=1, fy=2),
# ]
# addstage!(model, bcs)
# solve!(model)


# 3D beam

# coord = [ 0 0 0; 1 0 0; 0.5 0 0.5]
# n     = 2
# bl    = Block(coord, nx=n, cellshape = LIN3, tag="beam")
# mesh   = Mesh(bl, ndim=3)

geo = GeoModel()
p1  = add_point(geo, [0.0, 0.0, 0.0])
p2  = add_point(geo, [0.5, 0, 0.5])
p3  = add_point(geo, [1.0, 0.0, 0.0])
bz  = add_bezier(geo, [p1, p2, p3])
set_transfinite_curve(geo, bz, 3)

# add_arc(geo,
# add_block(geo, [0 0 0; 1 0 0], [0.5 0 0.5]; nx=2, shape=LIN3, tag="beam")
mesh = Mesh(geo, ndim=3, quadratic=true)

mapper = RegionModel(MechBeam, LinearElastic, E=1e4, nu=0.0, thy=0.1, thz=0.1)

# mats  = [ "beam" => MechBeam => LinearElastic => (E=1e4, nu=0, thy=0.1, thz=0.1) ]
# mats  = [ "beam" => MechBeam => ElasticBeam => (E=1e4, nu=0, A=0.01) ]

model = FEModel(mesh, mapper)
ana = MechAnalysis(model)

# changequadrature!(model.elems, 3)

add_monitor(ana, :node, x==0, :my)
add_monitor(ana, :node, x==0, :fx)
add_monitor(ana, :node, x==0, :fz)
add_monitor(ana, :node, x==1, :uz)

# monitors = [
#     :(x==0) => NodeMonitor(:(my))
#     :(x==0) => NodeMonitor(:(fx))
#     :(x==0) => NodeMonitor(:(fz))
#     :(x==1) => NodeMonitor(:(uz))
# ]
# setmonitors!(ana, monitors)

stage = add_stage(ana)
add_bc(stage, :node, x==0, rx=0,  ry=0, rz=0, ux=0, uy=0, uz=0)
add_bc(stage, :node, x==1, fx=1, fz=2)

run(ana, quiet=false)
# bcs =
# [
#    :(x==0) => NodeBC(rx=0, ry=0, rz=0, ux=0, uy=0, uz=0),
#    :(x==1) => NodeBC(fx=1, fz=2),
# ]
# addstage!(ana, bcs)
# solve!(ana, quiet=false)