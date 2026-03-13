using Serendip
using Test

ℓ = 10.0
c = 1500.0

geo = GeoModel()
add_block(geo, [0.0, 0.0], ℓ, ℓ, 0.0, nx=10, ny=10, shape=QUAD4, tag="fluid")
mesh = Mesh(geo)

mapper = RegionMapper()
add_mapping(mapper, "fluid", AcousticFluid, LinearAcousticFluid, rho=1.0, c=c)

model = FEModel(mesh, mapper, thickness=1.0)
ana = AcousticModalAnalysis(model)

stage = add_stage(ana)
add_bc(stage, :face, (y==ℓ), up=0.0)

run(ana, nmodes=10, quiet=false)

@test length(ana.frequencies) >= 3
# @test isapprox(ana.frequencies[1], 235.861759691578; rtol=1e-8)
# @test isapprox(ana.frequencies[2], 528.7052455577332; rtol=1e-8)
# @test isapprox(ana.frequencies[3], 713.4155966354814; rtol=1e-8)

ana.frequencies