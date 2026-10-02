using Serendip
using Test
using LinearAlgebra


@testset "ConcreteCompression constructor" begin
    @test_throws Exception ConcreteCompression(E=30_000.0, fc=-40.0, epsref=-0.002, Epost=20_000.0)
    @test_throws Exception ConcreteCompression(E=30_000.0, fc=40.0, epsref=-0.002, Epost=3_000.0)
end


@testset "ConcreteCompression uniaxial target curve" begin
    E = 32_400.0
    nu = 0.2
    fc = -47.9
    epsref = -0.00230
    Epost = 0.1E
    mat = ConcreteCompression(E=E, nu=nu, fc=fc, epsref=epsref, Epost=Epost)
    ctx = Context(ndim=3)

    fe = 0.4fc
    εe = fe/E
    εq = (fc-Epost*epsref)/(E-Epost)
    fq = E*εq
    ε̅cp = fc/E-epsref

    εcurve(t) = (1-t)^2*εe+2*(1-t)*t*εq+t^2*epsref
    σcurve(t) = (1-t)^2*fe+2*(1-t)*t*fq+t^2*fc
    Etcurve(t) = ((1-t)*(fq-fe)+t*(fc-fq))/((1-t)*(εq-εe)+t*(epsref-εq))

    cstate = Serendip.ConcreteCompressionState(ctx)
    for t in range(0.0, 1.0, length=11)
        εa = εcurve(t)
        σa = σcurve(t)
        εtarget = Serendip.Vec6(-nu*σa/E, -nu*σa/E, εa, 0.0, 0.0, 0.0)
        state = copy(cstate)
        _, status = update_state(mat, state, cstate, Vector(εtarget-cstate.ε))

        @test status.successful
        @test state.σ ≈ Serendip.Vec6(0.0, 0.0, σa, 0.0, 0.0, 0.0) atol=2e-9
        @test state.ε̅c ≈ ε̅cp*t^2 atol=2e-12
        @test -Serendip.tr(state.εp) ≈ state.ε̅c atol=2e-12

        if 0.0 < t < 1.0
            D = Matrix(Serendip.calcD(mat, state))
            Et = D[3,3]-dot(D[3,1:2], D[1:2,1:2]\D[1:2,3])
            @test Et ≈ Etcurve(t) rtol=2e-10
        end
        cstate = copy(state)
    end

    εa = epsref-0.0005
    σa = fc+Epost*(εa-epsref)
    εtarget = Serendip.Vec6(-nu*σa/E, -nu*σa/E, εa, 0.0, 0.0, 0.0)
    state = copy(cstate)
    _, status = update_state(mat, state, cstate, Vector(εtarget-cstate.ε))
    @test status.successful
    @test state.σ[3] ≈ σa atol=2e-9
    @test state.ε̅c ≈ σa/E-εa atol=2e-12

    D = Matrix(Serendip.calcD(mat, state))
    Et = D[3,3]-dot(D[3,1:2], D[1:2,1:2]\D[1:2,3])
    @test Et ≈ Epost rtol=2e-10

    cstate = copy(state)
    unload_stress = 5.0
    Δε = Matrix(Serendip.calcDe(E, nu))\Serendip.Vec6(0.0, 0.0, unload_stress, 0.0, 0.0, 0.0)
    state = copy(cstate)
    _, status = update_state(mat, state, cstate, Vector(Δε))
    @test status.successful
    @test state.active == 0
    @test state.ε̅c == cstate.ε̅c
    @test Serendip.calcD(mat, state) == Serendip.calcDe(E, nu)
end


@testset "ConcreteCompression unit scaling" begin
    nu = 0.2
    epsref = -0.00230
    t = 0.65

    for scale in (1.0, 1e3, 1e6)
        E = 32_400.0*scale
        fc = -47.9*scale
        Epost = 0.1E
        mat = ConcreteCompression(E=E, nu=nu, fc=fc, epsref=epsref, Epost=Epost)
        fe = 0.4fc
        εe = fe/E
        εq = (fc-Epost*epsref)/(E-Epost)
        fq = E*εq
        εa = (1-t)^2*εe+2*(1-t)*t*εq+t^2*epsref
        σa = (1-t)^2*fe+2*(1-t)*t*fq+t^2*fc
        εtarget = Serendip.Vec6(-nu*σa/E, -nu*σa/E, εa, 0.0, 0.0, 0.0)

        cstate = Serendip.ConcreteCompressionState(Context(ndim=3))
        state = copy(cstate)
        _, status = update_state(mat, state, cstate, Vector(εtarget))
        @test status.successful
        @test state.active == 1
        @test state.σ[3] ≈ σa rtol=2e-12
    end
end


@testset "ConcreteCompression active sets" begin
    E = 30_000.0
    nu = 0.2
    mat = ConcreteCompression(E=E, nu=nu, fc=-40.0, epsref=-0.0023, Epost=3_000.0)
    ctx = Context(ndim=3)
    De = Matrix(Serendip.calcDe(E, nu))

    trial_stresses = (
        Serendip.Vec6(-20.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        Serendip.Vec6(-20.0, -20.0, 0.0, 0.0, 0.0, 0.0),
        Serendip.Vec6(-20.0, -20.0, -20.0, 0.0, 0.0, 0.0),
    )

    for (expected_active, σtr) in enumerate(trial_stresses)
        cstate = Serendip.ConcreteCompressionState(ctx)
        state = copy(cstate)
        Δε = De\σtr
        _, status = update_state(mat, state, cstate, Vector(Δε))

        @test status.successful
        @test state.active == expected_active
        @test state.ε̅c > 0.0
        @test -Serendip.tr(state.εp) ≈ state.ε̅c atol=2e-12

        principal = sort(collect(Serendip.eigvals(state.σ)))
        current_fc = Serendip.compression_strength(mat, state.ε̅c)
        @test all(abs(principal[a]-current_fc) <= 2e-9 for a in 1:expected_active)
        @test all(principal[a] >= current_fc-2e-9 for a in expected_active+1:3)
        D = Matrix(Serendip.calcD(mat, state))
        @test D ≈ D'

        direction = -sum(state.projectors[:,1:expected_active], dims=2)[:]
        h = 1e-8
        next_state = copy(state)
        Δσ, next_status = update_state(mat, next_state, state, h*direction)
        @test next_status.successful
        @test next_state.active == expected_active
        @test Δσ/h ≈ D*direction rtol=2e-4
    end

    cstate = Serendip.ConcreteCompressionState(ctx)
    state = copy(cstate)
    tensile_trial = Serendip.Vec6(10.0, 4.0, 2.0, 0.0, 0.0, 0.0)
    _, status = update_state(mat, state, cstate, Vector(De\tensile_trial))
    @test status.successful
    @test state.active == 0
    @test state.ε̅c == 0.0
    @test state.σ ≈ tensile_trial
end


@testset "ConcreteCompression rejects plane stress" begin
    mat = ConcreteCompression(E=30_000.0, nu=0.2, fc=-40.0, epsref=-0.0023, Epost=3_000.0)
    ctx = Context(ndim=2, stress_state=:plane_stress)
    cstate = Serendip.ConcreteCompressionState(ctx)
    state = copy(cstate)
    @test_throws ErrorException update_state(mat, state, cstate, zeros(6))
    @test_throws ErrorException Serendip.calcD(mat, state)
end


@testset "ConcreteCompression solid mapping" begin
    geo = GeoModel()
    add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 1.0, nx=1, ny=1, nz=1, shape=:hex8, tag="solid")
    mesh = Mesh(geo, quiet=true)
    mapper = RegionMapper()
    add_mapping(
        mapper,
        "solid",
        MechSolid,
        ConcreteCompression,
        E=32.4e6,
        nu=0.2,
        fc=-47.9e3,
        epsref=-0.0023,
        Epost=3.24e6,
        rho=2.4,
    )
    model = FEModel(mesh, mapper, quiet=true)

    @test model.elems[1].cmodel isa ConcreteCompression
    @test model.elems[1].etype.ρ == 2.4
    @test model.elems[1].ips[1].state isa Serendip.ConcreteCompressionState
end
