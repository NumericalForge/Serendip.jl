using Serendip
using Test

# Mesh:
geo = GeoModel()

add_block(geo, [0.0, 0.0, 0.0], [1.0, 6.0, 1.0], nx=1, ny=10, nz=1, tag="solids")
p1 = add_point(geo, [0.5, 3.0, 0.2], tag="bars")
p2 = add_point(geo, [0.5, 6.0, 0.2], tag="bars")
l1 = add_line(geo, p1, p2)

p1 = add_point(geo, [0.5, 3.0, 0.8], tag="bars")
p2 = add_point(geo, [0.5, 6.0, 0.8], tag="bars")
l2 = add_line(geo, p1, p2)

add_path(geo, [l1], tag="bars", interface_tag="joints")
add_path(geo, [l2], tag="bars", interface_tag="joints", tips=:both, tip_tag="tips")

mesh = Mesh(geo)
save(mesh, "mesh.vtk")

# Finite elements:
mapper = RegionMapper()
add_mapping(mapper, "solids", MechBulk, LinearElastic, E=24e3, nu=0.2)
add_mapping(mapper, "bars", MechBar, LinearElastic, E=200e6, A=0.00011)
add_mapping(mapper, "joints", MechBondSlip, CebBondSlip, taumax=12, taures=3, s1=0.001, s2=0.0011, s3=0.004, alpha=0.5, beta= 0.5, ks=(12/0.001)*5, kn=5000, p=0.25)
add_mapping(mapper, "tips", MechBondTip, LinearTipContact, k=1e8)

model = FEModel(mesh, mapper)
ana   = MechAnalysis(model)

stage = add_stage(ana, nincs=20)
add_bc(stage, :node, (y==0, z==0), ux=0, uy=0, uz=0)
add_bc(stage, :node, (y==6, z==0), uz=0)
add_bc(stage, :edge, (y==3, z==0), qz=-1.0)

run(ana, autoinc=true, scheme=:Ralston, tol=0.01, maxits=3).success
