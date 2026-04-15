# Modal analysis of a fluid block with a free surface.

using Serendip
using Test

h = 1.0
c = 1500.0

geo = GeoModel()
add_block(geo, [0.0, 0.0], 0.1, h, 0.0, nx=1, ny=10, shape=:quad4, tag="wall")
add_block(geo, [0.1, 0.0], 10.0, h, 0.0, nx=10, ny=10, shape=:quad4, tag="fluid")

mesh = Mesh(geo)
add_contact_elements(mesh, "wall", "fluid", tag="interface")

mapper = RegionMapper()
add_mapping(mapper, "wall", MechSolid, LinearElastic, E=211e6, nu=0.3, rho=7.8)
add_mapping(mapper, "fluid", AcousticFluid, LinearAcousticFluid, rho=1.0, c=c)
add_mapping(mapper, "interface", AcousticMechInterface, AcousticMechInterfaceCoupling)


model = FEModel(mesh, mapper, thickness=1.0, g=9.81)
ana = AcousticModalAnalysis(model)

stage = add_stage(ana)
add_bc(stage, :face, (x<=0.1, y==0), ux=0.0, uy=0.0)
add_bc(stage, :face, (x>=0.1, y==h), fs=true)

run(ana, nmodes=10, eigen_solver=:lapack)

ana.frequencies

# @test length(ana.frequencies) >= 3
# @test isapprox(ana.frequencies[1], 59; rtol=1e+1)
# @test isapprox(ana.frequencies[2], 704; rtol=1e+1)
# @test isapprox(ana.frequencies[3], 2361; rtol=1e+1)
