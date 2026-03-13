using Serendip
using Test

geo = GeoModel()
add_block(geo, [0.0, 0.0, 0.0], 1.0, 2.0, 0.0, nx=10, ny=10, shape=QUAD4, tag="solids")
mesh = Mesh(geo)

k = 0.0502
rho = 7.8
cv = 486.0
E = 200e6
nu = 0.3
alpha = 1.2e-5

mapper = RegionMapper()
add_mapping(mapper, "solids", TMSolid, LinearElasticThermo; E=E, nu=nu, k=k, alpha=alpha, rho=rho, cv=cv)

model = FEModel(mesh, mapper, T0=0.0, stress_state=:plane_strain)

mktempdir() do outdir
    ana = ThermoMechAnalysis(model; outdir=outdir)

    stage = add_stage(ana; tspan=3_000_000.0, nincs=10, nouts=2)
    add_bc(stage, :node, (x == 0), ut=10.0)
    add_bc(stage, :node, (y == 2), ut=20.0)
    add_bc(stage, :node, (y == 0), ux=0.0, uy=0.0)
    add_bc(stage, :node, (x == 1), fx=100.0)

    status = run(ana; tol=0.1, quiet=true)

    @test status.successful
    @test isfinite(maximum(abs, model.node_fields["ux"]))
    @test minimum(model.node_fields["ut"]) >= 10.0 - 1e-8
    @test maximum(model.node_fields["ut"]) <= 20.0 + 1e-8
end
