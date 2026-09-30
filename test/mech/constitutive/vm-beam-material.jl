using Serendip
using Test
using LinearAlgebra


@testset "normalized von Mises beam update" begin
    E, nu = 200_000.0, 0.3
    mat = VonMises(E=E, nu=nu, fy=250.0, H=1_000.0)
    ctx = Context(ndim=3)

    G  = E/(2*(1 + nu))
    De = [E, 2*G, 2*G]
    Pd = [2/3, 1.0, 1.0]

    increments = (
        [0.002, 0.0, 0.0],
        [0.0, 0.002, 0.0],
        [0.002, 0.001, 0.0005],
    )

    for Δε in increments
        cstate = Serendip.VonMisesBeamState(ctx)
        state  = copy(cstate)

        σtr    = De.*Δε
        _, status = update_state(mat, state, cstate, Δε)
        @test status.successful
        @test state.εpa - cstate.εpa ≈ state.Δλ

        σvm = √(state.σ[1]^2 + 1.5*(state.σ[2]^2 + state.σ[3]^2))
        Δγ  = 3*state.Δλ/(2*σvm)
        A   = 1.0 .+ Δγ.*De.*Pd

        @test state.σ ≈ σtr./A atol=1e-4 rtol=1e-7
        @test (σtr-state.σ)./De ≈ Δγ.*Pd.*state.σ atol=1e-9
        @test σvm - mat.σy - mat.H*state.εpa ≈ 0.0 atol=1e-6*mat.σy
    end
end


@testset "nonproportional von Mises beam path" begin
    E, nu = 200_000.0, 0.3
    mat = VonMises(E=E, nu=nu, fy=250.0, H=1_000.0)
    ctx = Context(ndim=3)
    path = (
        [0.0030, 0.0, 0.0],
        [0.0, 0.0, 0.0040/√2],
        [-0.0045, 0.0, 0.0],
        [0.0, 0.0, -0.0060/√2],
    )

    function integrate_path(nsub)
        cstate = Serendip.VonMisesBeamState(ctx)
        cstate.αs = 5/6

        for segment in path
            Δε = segment/nsub
            for i in 1:nsub
                state = copy(cstate)
                _, status = update_state(mat, state, cstate, Δε)
                status.successful || error("von Mises beam material update failed")
                cstate = state
            end
        end
        return cstate
    end

    reference = integrate_path(2_000)
    coarse    = integrate_path(10)
    refined   = integrate_path(100)
    scale     = max(norm(reference.σ), mat.σy)

    @test norm(refined.σ-reference.σ)/scale < norm(coarse.σ-reference.σ)/scale
    @test norm(refined.σ-reference.σ)/scale < 5e-3
    @test reference.εpa > 0.0
end


@testset "von Mises beam consistent tangent" begin
    E, nu = 200_000.0, 0.3
    mat = VonMises(E=E, nu=nu, fy=250.0, H=1_000.0)
    ctx = Context(ndim=3)

    cstate = Serendip.VonMisesBeamState(ctx)
    state  = copy(cstate)
    _, status = update_state(mat, state, cstate, [0.002, 0.001, 0.0005])
    @test status.successful
    @test state.Δλ > 0.0

    Dalg = Matrix(Serendip.calcD(mat, state))
    Dfd  = zeros(3, 3)
    h    = 1e-8

    for j in 1:3
        Δεp = [0.002, 0.001, 0.0005]
        Δεm = copy(Δεp)
        Δεp[j] += h
        Δεm[j] -= h

        statep = copy(cstate)
        statem = copy(cstate)
        _, statusp = update_state(mat, statep, cstate, Δεp)
        _, statusm = update_state(mat, statem, cstate, Δεm)
        @test statusp.successful && statusm.successful
        @test statep.Δλ > 0.0 && statem.Δλ > 0.0
        Dfd[:,j] = (statep.σ-statem.σ)/(2h)
    end

    @test Dalg ≈ Dfd rtol=1e-7
    @test Dalg ≈ Dalg' rtol=1e-12
end
