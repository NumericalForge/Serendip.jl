using Serendip
using Test

coord = [ 0. 0.; 1. 0. ]
conn  = [ [1, 2] ]

mesh = Mesh(coord, conn)

mapper = RegionModel(MechBar, VonMises, E=210e6, fy=500e3, A=0.01)

# tag!(mesh.elems, "bars")

# mats = [
#         # "bars" => MechBar => PPBar => (E=210e6, fy=500e3, A=0.01),
#         # "bars" => MechBar => VonMises => (E=210e6, fy=500e3, A=0.01),
#         "bars" => MechBar => PPRod => (E=210e6, fy=500e3, A=0.01),
#        ]

ctx = Context()
model = FEModel(mesh, mapper, ctx)
ana = MechAnalysis(model)
stage = add_stage(ana)


# bcs = [
#     :(x==0 && y==0) => NodeBC(ux=0, uy=0),
#     :(x==1 && y==0) => NodeBC(uy=0),
#     :(x==1 && y==0) => NodeBC(ux=0.003),
#     ]
# addstage!(ana, bcs, nincs=10)


# @test solve!(ana).success

add_bc(stage, :node, (x==0, y==0), ux=0, uy=0)
add_bc(stage, :node, (x==1, y==0), uy=0)
add_bc(stage, :node, (x==1, y==0), ux=0.003)
@test run(ana).success
