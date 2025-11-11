using Serendip
using Test

geo = GeoModel()
add_block(geo, [0,0,0], 1,1,1, nx=2, ny=2, nz=2, tag="solids")
mesh = Mesh(geo)

mapper = RegionMapper()
add_mapping(mapper, "solids", MechBulk, LinearElastic, E=100.0, nu=0.2)

model = FEModel(mesh, mapper)
ana = MechAnalysis(model)

add_logger(ana, :node, (x==1, y==1, z==1), "node.dat")
add_logger(ana, :node, (x==1, y==0))
add_logger(ana, :nodegroup, (x==1, y==1), "nodes.dat")

# add_logger(ana, :face, z==1, "face.dat")
# add_logger(ana, :edge, z==1, "edge.dat")
add_logger(ana, :ip, (x>0.5, y>0.5, z>0.5), "ip.dat")

add_logger(ana, :ipgroup, (x>0.5, y>0.5), "ips.dat")
add_logger(ana, :ip, [0.5,0.5,0.0])

# add_logger(ana, :point, [0.5,0.5,0.0], "point.dat")
# add_logger(ana, :segment, [0.5 0.5 0; 1 1 0], "segment.dat")

stage = add_stage(ana, nincs=4, nouts=4)
add_bc(stage, :node, z==0, ux=0, uy=0, uz=0 )
add_bc(stage, :face, z==1, tz=-10.0)

@test run(ana).successful