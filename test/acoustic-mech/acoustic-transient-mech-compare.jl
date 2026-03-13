using Serendip
using Test
using LinearAlgebra: norm

function build_cantilever_model()
    L = 2.0
    H = 0.2
    W = 0.2

    E = 30e6
    rho = 2.4

    geo = GeoModel()
    add_block(geo, [0.0, 0.0, 0.0], L, W, H, nx=10, ny=1, nz=1, shape=HEX20, tag="solids")
    mesh = Mesh(geo)

    mapper = RegionMapper()
    add_mapping(mapper, "solids", MechBulk, LinearElastic, E=E, nu=0.2, rho=rho)

    model = FEModel(mesh, mapper)
    return model, (; L, H, W, E, rho)
end

function run_cantilever(analysis_type)
    model, pars = build_cantilever_model()
    analysis_type === AcousticMechAnalysis && (model.ctx.stress_state = :d3)
    ana = analysis_type(model)
    log = add_logger(ana, :nodalreduce, (x == pars.L))

    stage = add_stage(ana, nincs=300, tspan=0.3)
    add_bc(stage, :node, (x == 0.0), ux=0.0, uy=0.0, uz=0.0)
    add_bc(stage, :face, (x == pars.L), tz=-1.0 / (pars.H * pars.W))

    status = run(ana, quiet=true)
    return status, log.table, pars
end

status_dyn, table_dyn, pars = run_cantilever(DynamicAnalysis)
status_acm, table_acm, _ = run_cantilever(AcousticMechAnalysis)

@test status_dyn.successful
@test status_acm.successful

t_dyn = Float64.(table_dyn[:t])
uz_dyn = Float64.(table_dyn[:uz])
t_acm = Float64.(table_acm[:t])
uz_acm = Float64.(table_acm[:uz])

@test length(t_dyn) == length(t_acm)
@test t_dyn ≈ t_acm

u_min_dyn = minimum(uz_dyn)
u_min_acm = minimum(uz_acm)
@test 0.70 <= u_min_acm / u_min_dyn <= 0.85
@test 1.00 <= norm(uz_acm) / norm(uz_dyn) <= 1.30

u_static = (-1.0) * pars.L^3 / (3 * pars.E * (pars.W * pars.H^3 / 12))
@test isapprox(u_min_dyn, 2 * u_static; rtol=0.15)
