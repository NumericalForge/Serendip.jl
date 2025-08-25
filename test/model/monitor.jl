using Serendip
using Test

geo = GeoModel()
add_block(geo, [0,0,0], [1,1,1], nx=2, ny=2, nz=2, tag="solids")
mesh = Mesh(geo)

mapper = RegionMapper()
add_mapping(mapper, "solids", MechBulk, LinearElastic, E=100.0, nu=0.2)

model = FEModel(mesh, mapper)
ana = MechAnalysis(model)

add_monitor(ana, :node, (x==1, y==1, z==1), :uy, "node.table")
add_monitor(ana, :node, (x==1, y==0), :uy)
add_monitor(ana, :nodalreduce, (x==1, y==1), :uy, "nodes.table")

add_monitor(ana, :ip, (x>0.5, y>0.5, z>0.5), :σyy, "ip.table")
add_monitor(ana, :ip, [0.5,0.5,0.0], :σyy)
add_monitor(ana, :ipgroup, (x>0.5, y>0.5), :σyy, "ips.table")


# bl = Block( [0 0 0; 1 1 1], nx=2, ny=2, nz=2, shape=HEX8, tag="solids")
# mesh = Mesh(bl)

# mats = [ "solids" => MechBulk => LinearElastic => (E=100.0, nu=0.2) ]

# ctx = Context()
# model = FEModel(mesh, mats, ctx)
# ana = MechAnalysis(model)

# # Monitors
# monitors = [
#     (x==1, y==1, z==1) => NodeMonitor(:uy, "node.table")
#     (x==1, y==0) => NodeMonitor(:uy)
#     y==1 => NodeSumMonitor(:uy, "nodessum.book")

#     (x>0.5, y>0.5, z>0.5) => IpMonitor(:σyy, "ip.table")
#     (x>0.5, y>0.5) => IpGroupMonitor(:σyy, "ips.book")
#     [0.5,0.5,0.0] => IpMonitor(:σyy)
# ]

# setloggers!(ana, loggers)

stage = add_stage(ana, nincs=4, nouts=4)
add_bc(stage, :node, z==0, ux=0, uy=0, uz=0 )
add_bc(stage, :face, z==1, tz=-10.0)

@test run(ana).success