using Serendip
using Test
using LinearAlgebra


@testset "normalized von Mises solid hardening update" begin
    E, nu = 200_000.0, 0.3
    fy, H = 250.0, 1_000.0
    mat = VonMises(E=E, nu=nu, fy=fy, H=H)
    ctx = Context(ndim=3)
    De = Matrix(Serendip.calcDe(E, nu))
    G = E/(2*(1 + nu))

    increments = (
        [0.004, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.004],
        [0.003, -0.001, 0.0005, 0.001, 0.0007, -0.0008],
    )

    for Δε in increments
        cstate = Serendip.VonMisesState(ctx)
        state = copy(cstate)
        σtr = Serendip.Vec6(De*Δε)
        ftr = √(3*Serendip.J2(σtr)) - fy
        expected_Δλ = ftr/(3*G + H)
        @test ftr > 0.0

        _, status = update_state(mat, state, cstate, Δε)
        @test status.successful
        @test state.Δλ ≈ expected_Δλ
        @test state.εpa - cstate.εpa ≈ state.Δλ

        Δεp = De \ (σtr - state.σ)
        @test √(2/3)*norm(Δεp) ≈ state.Δλ
        @test abs(√(3*Serendip.J2(state.σ)) - fy - H*state.εpa) <= 1e-6*fy

        s = Serendip.dev(state.σ)
        g = √1.5*s/norm(s)
        De_g = De*g
        D_cont = De - (De_g*De_g')/(dot(g, De_g) + H)
        @test Matrix(Serendip.calcD(mat, state)) ≈ D_cont
    end

    cstate = Serendip.VonMisesState(ctx)
    for Δε in fill([0.002, 0.0, 0.0, 0.0, 0.0, 0.0], 3)
        state = copy(cstate)
        σtr = Serendip.Vec6(cstate.σ + De*Δε)
        ftr = √(3*Serendip.J2(σtr)) - fy - H*cstate.εpa
        @test ftr > 0.0

        _, status = update_state(mat, state, cstate, Δε)
        @test status.successful
        @test state.Δλ ≈ ftr/(3*G + H)
        @test state.εpa - cstate.εpa ≈ state.Δλ

        Δεp = De \ (σtr - state.σ)
        @test √(2/3)*norm(Δεp) ≈ state.Δλ
        @test abs(√(3*Serendip.J2(state.σ)) - fy - H*state.εpa) <= 1e-6*fy
        cstate = copy(state)
    end
end


th = 0.05
E  = 210e6 # kPa
fy = 240e3 # kPa
H  = 0
nu = 0.3

geo = GeoModel()
add_block(geo, [0, 0, -th], th, 1.0, 2*th, nx=1, ny=15, nz=2, shape=:hex20, tag="beam")
mesh = Mesh(geo)


mapper = RegionMapper()
add_mapping(mapper, "beam", MechSolid, VonMises, E=E, nu=nu, fy=fy, H=H)
model = FEModel(mesh, mapper)

ana = MechAnalysis(model)
log = add_logger(ana, :node, (x==th/2, y==1, z==0))
add_monitor(ana, :node, (x==th/2, y==1, z==0), :fz)

stage = add_stage(ana, nincs=20, nouts=1)
add_bc(stage, :node, (y==0), uy=0)
add_bc(stage, :node, (y==0, z==0), uz=0)
add_bc(stage, :node, (x==th/2, y==0, z==0), ux=0)
add_bc(stage, :node, (x==th/2, y==1, z==0), uz=-0.03)

run(ana, autoinc=true, tol=0.01)
@test log.table.fz[end]≈-31.0 atol=0.1

makeplots = false
if @isdefined(makeplots) && makeplots
    tab = log.table
    chart = Chart(;
        xlabel = "Displacement uz [m]",
        ylabel = "Force fz [kN]",
    )
    add_series(chart, -tab["uz"], -tab["fz"], mark=:circle)
    save(chart, "vm-3d.pdf")
end
