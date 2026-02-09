using Serendip
using Test

coord = [ 0. 0.; 1. 0. ]
conn  = [ [1, 2] ]

mesh = Mesh(coord, conn)

mapper = RegionModel(MechBar, VonMises, E=210e6, fy=500e3, A=0.01)

model = FEModel(mesh, mapper)
ana   = MechAnalysis(model)

add_monitor(ana, :ip, [0,0,0], :σx´)
add_monitor(ana, :node, [1,0,0], :ux)
add_monitor(ana, :node, [1,0,0], :fx)
log = add_logger(ana, :node, [1,0,0])

stage = add_stage(ana)

add_bc(stage, :node, (x==0, y==0), ux=0, uy=0)
add_bc(stage, :node, (x==1, y==0), uy=0)
add_bc(stage, :node, (x==1, y==0), fx=500)

run(ana, autoinc=true).successful

@test log.table[:ux][end] ≈ 0.000238 atol=1e-6
