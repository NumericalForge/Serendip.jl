using Serendip
using Test

geo = GeoModel()
add_block(geo, [0.0, 0.0, 0.0], 1.0, 2.0, 0.0, nx=1, ny=4, shape=QUAD4, tag="solids")
mesh = Mesh(geo)

k = 100.62
rho = 1.6
cv = 486e3

mapper = RegionMapper()
add_mapping(mapper, "solids", ThermoSolid, ConstConductivity; k=k, rho=rho, cv=cv)

model = FEModel(mesh, mapper, T0=0.0)

ana = ThermoAnalysis(model)

stage = add_stage(ana; tspan=10_000.0, nincs=5, nouts=1)
add_bc(stage, :node, (y == 0), ut=10.0)
add_bc(stage, :node, (y == 2), ut=20.0)

status = run(ana; tol=0.1, quiet=true)

@test status.successful
@test minimum(model.node_fields["ut"]) ≈ 10.0 atol=1e-8
@test maximum(model.node_fields["ut"]) ≈ 20.0 atol=1e-8
