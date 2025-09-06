using Serendip
using Test

geo = GeoModel()
add_block(geo, [0.0, 0.0, 0.0], [1.0, 1.0, 0.5], nx=2, ny=2, nz=2, tag="solids")
mesh = Mesh(geo)

mapper = RegionMapper()
add_mapping(mapper, "solids", MechBulk, DruckerPrager, E=100., nu=0.25, alpha=0.05, kappa=0.1)

model = FEModel(mesh, mapper)
ana   = MechAnalysis(model)
stage = add_stage(ana, nincs=10)

add_bc(stage, :node, (z==0.0), ux=0, uy=0, uz=0)
add_bc(stage, :node, (z==0.5), uz=-0.033)
add_bc(stage, :node, :(x==0 || x==1.0), ux=0, uy=0)
add_bc(stage, :node, :(y==0 || y==1.0), ux=0, uy=0)

stage = add_stage(ana, nincs=10)
add_bc(stage, :node, (z==0.0), ux=0, uy=0, uz=0)
add_bc(stage, :node, (z==0.5), uz=+0.008)
add_bc(stage, :node, :(x==0 || x==1.0), ux=0, uy=0)
add_bc(stage, :node, :(y==0 || y==1.0), ux=0, uy=0)

@test run(ana, autoinc=true, tol=1e-2).success
