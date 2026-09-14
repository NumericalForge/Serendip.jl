using Serendip
using Test
using LinearAlgebra

h  = 0.1
th = 0.05
L  = 1.0
E  = 210e6 # kPa
fy = 240e3 # kPa
H  = 0
nu = 0.3

geo = GeoModel()
add_block(geo, [0.0, 0.0], L, h, 0, nx=50, ny=2, shape=:quad8, tag="beam")
mesh = Mesh(geo)

mapper = RegionMapper()
add_mapping(mapper, "beam", MechSolid, VonMises, E=E, nu=nu, fy=fy, H=H)

model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=th)

ana   = MechAnalysis(model)
log = add_logger(ana, :node, (y==h/2, x==1))
mon = add_monitor(ana, :node, (y==h/2, x==1), :fy)

stage = add_stage(ana, nincs=30, nouts=1)
add_bc(stage, :node, (x==0), ux=0)
add_bc(stage, :node, (x==0, y==h/2), uy=0)
add_bc(stage, :node, (x==1, y==h/2), uy = -0.03)

run(ana, autoinc=true)

@test log.table["fy"][end] ≈ -30 atol=0.6

makeplots = false
if @isdefined(makeplots) && makeplots
    tab = log.table
    chart = Chart(;
        xlabel = "Displacement uy [m]",
        ylabel = "Force fy [kN]",
    )
    add_series(chart, -tab["uy"], -tab["fy"], mark=:circle)
    save(chart, "vm-2d.pdf")
end


@testset "normalized plane-stress hardening" begin
    ctx = Context(ndim=2, stress_state=:plane_stress)
    mat = VonMises(E=200_000.0, nu=0.3, fy=250.0, H=1_000.0)

    cstate = Serendip.VonMisesState(ctx)
    state  = copy(cstate)
    Δε     = [0.01, -0.002, 0.0, 0.0, 0.0, 0.003]

    _, status = update_state(mat, state, cstate, Δε)
    @test status.successful
    @test state.ε ≈ Δε
    @test state.εpa - cstate.εpa ≈ state.Δλ
    @test abs(√(3*Serendip.J2(state.σ)) - mat.σy - mat.H*state.εpa) <= 1e-6*mat.σy

    De = Serendip.calcDe(mat.E, mat.ν, :plane_stress)
    s  = Serendip.dev(state.σ)
    n  = √1.5*s/norm(s)
    D  = De - De*n*n'*De/(dot(n, De*n) + mat.H)
    @test Serendip.calcD(mat, state) ≈ D

    cstate = copy(state)
    state  = copy(cstate)
    _, status = update_state(mat, state, cstate, 0.5*Δε)
    @test status.successful
    @test state.εpa >= cstate.εpa
    @test state.εpa - cstate.εpa ≈ state.Δλ
    @test abs(√(3*Serendip.J2(state.σ)) - mat.σy - mat.H*state.εpa) <= 1e-6*mat.σy
end
