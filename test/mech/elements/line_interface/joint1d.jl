using Serendip
using Test

geo = GeoModel()
add_block(geo, [0.0, 0.0, 0.0], 1.0, 6.0, 1.0, nx=4, ny=6, nz=2, tag="solids")
p1 = add_point(geo, [0.2, 0.2, 0.2])
p2 = add_point(geo, [0.2, 5.8, 0.2])
edge = add_line(geo, p1, p2)
path = add_path(geo, [edge]; tag="bars", interface_tag="interface")
add_array(geo, path; nx=3, dx=0.3)

# Mesh generation
# bl  = Block( [0 0 0; 1.0 6.0 1.0], nx=4, ny=6, nz=2, tag="solids")
# bl1 = BlockInset( [0.2 0.2 0.2; 0.2 5.8 0.2], curvetype="polyline", tag="bars", jointtag="interface")
# bl2 = move( copy(bl1), dx=0.6)
# bl3 = move( copy(bl1), dx=0.3)
# bls = [bl, bl1, bl2, bl3 ]

mesh = Mesh(geo)
save(mesh, "mesh.vtu")

# FEM analysis
mapper = RegionMapper()
add_mapping(mapper, "solids", MechBulk, LinearElastic, E=1.e4, nu=0.)
add_mapping(mapper, "bars", MechBar, LinearElastic, E=1.e8, A=0.005)
add_mapping(mapper, "interface", MechBondSlip, LinearBondSlip, ks=1.e5, kn=1.e5, p=0.25)

# mats = [
    # "solids" => MechBulk => LinearElastic => (E=1.e4, nu=0.),
    # "interface" => MechBondSlip => ElasticRSJoint => (ks=1.e5, kn=1.e5, p=0.25),
    # "bars"   => MechBar => LinearElastic => (E=1.e8, A=0.005),
# ]

model = FEModel(mesh, mapper)
ana = MechAnalysis(model)

add_logger(ana, :node, (x==0.5, y==3.0, z==0.5))

stage = add_stage(ana, nincs=5)
add_bc(stage, :node, (y==0.0, z==0.0), ux=0, uy=0, uz=0)
add_bc(stage, :node, (y==6.0, z==0.0), uz=0.0)
add_bc(stage, :face, z==1.0, tz=-1000)

@test run(ana).successful

# bcs = [
    #    :(y==0 && z==0) => NodeBC(uy=0, uz=0),
    #    :(y==6 && z==0) => NodeBC(uz=0),
    #    :(z==1) => SurfaceBC(tz=-1000 ),
    #   ]

# addlogger!(ana, :(x==0.5 && y==3.0 && z==0.5) => NodeLogger() )
# addstage!(ana, bcs, nincs=20)
# @test solve!(ana).successful