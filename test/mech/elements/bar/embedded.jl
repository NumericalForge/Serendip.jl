using Serendip
using Test

# Mesh generation
b = 0.5
h = 0.6
ℓ = 4.0

geo = GeoModel()
add_block(geo, [0,0,0], b, ℓ, h; nx=1, ny=25, nz=3, tag="solid")

p1 = add_point(geo, [0.2, 0.1, 0.1])
p2 = add_point(geo, [0.2, 3.6, 0.1])
pc = add_point(geo, [0.2, 3.6, 0.4])
p3 = add_point(geo, [0.2, 3.9, 0.4])

l1 = add_line(geo, p1, p2)
arc = add_circle_arc(geo, p2, pc, p3)
add_path(geo, [l1, arc]; embedded=true, tag="embedded", shape=LIN3)

mesh = Mesh(geo)

# FEM analysis
mapper = RegionMapper()
add_mapping(mapper, "solid", MechBulk, LinearElastic, E=1.e4, nu=0.25)
add_mapping(mapper, "embedded", MechBar, VonMises, E=1.e8, fy=500e3, A=0.005)

model = FEModel(mesh, mapper)
ana   = MechAnalysis(model)

stage = add_stage(ana)
add_bc(stage, :node, (y==0, z==0), ux=0, uy=0, uz=0)
add_bc(stage, :node, (y==ℓ, z==0), ux=0, uy=0, uz=0)
add_bc(stage, :face, (z==h), tz=-1000)

run(ana, autoinc=true)

save(model, "embedded.vtu")
