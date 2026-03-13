# Modal analysis of a fluid block with a free surface.

using Serendip
using Test

ℓ = 10.0
h = 10.0
c = 1500.0

geo = GeoModel()
add_block(geo, [0.0, 0.0], ℓ, h, 0.0, nx=10, ny=10, shape=QUAD4, tag="fluid")
mesh = Mesh(geo)

mapper = RegionMapper()
add_mapping(mapper, "fluid", AcousticFluid, LinearAcousticFluid, rho=1.0, c=c)

model = FEModel(mesh, mapper, thickness=1.0, g=9.81)
ana = AcousticModalAnalysis(model)

stage = add_stage(ana)
add_bc(stage, :face, (x==h), fs=true)

run(ana, nmodes=10, eig_method=:arpack)

@test length(ana.frequencies) >= 3
# @test isapprox(ana.frequencies[1], 1.75965313; rtol=1e-8)
# @test isapprox(ana.frequencies[2], 2.52412127; rtol=1e-8)
# @test isapprox(ana.frequencies[3], 3.15675566; rtol=1e-8)

ana.frequencies
