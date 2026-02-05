using Serendip
using Test

coord = [ 0. 0.; 1. 0. ]
conn  = [ [1, 2] ]

mesh = Mesh(coord, conn)

mapper = RegionModel(MechBar, VonMises, E=210e6, fy=500e3, A=0.01)

model = FEModel(mesh, mapper)
ana   = MechAnalysis(model)
stage = add_stage(ana)

add_bc(stage, :node, (x==0, y==0), ux=0, uy=0)
add_bc(stage, :node, (x==1, y==0), uy=0)
add_bc(stage, :node, (x==1, y==0), ux=0.003)
@test run(ana).successful
