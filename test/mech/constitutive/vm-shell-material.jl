using Serendip
using Test
using LinearAlgebra


@testset "von Mises shell consistent tangent" begin
    E, nu = 200_000.0, 0.3
    fy    = 250.0
    αs    = 5/6
    ctx   = Context(ndim=3, stress_state=:plane_stress)
    active = [1, 2, 4, 5, 6]
    Δε = [0.002, -0.0003, 0.0, 0.0005, 0.0002, 0.001]

    for H in (0.0, 1_000.0)
        mat = VonMises(E=E, nu=nu, fy=fy, H=H)
        cstate = Serendip.VonMisesState(ctx)
        cstate.αs = αs
        state = copy(cstate)

        _, status = update_state(mat, state, cstate, Δε)
        @test status.successful
        @test state.Δλ > 0.0

        Dalg = Matrix(Serendip.calcD(mat, state))
        Dfd  = zeros(6, 6)
        h    = 1e-8

        for j in active
            Δεp = copy(Δε)
            Δεm = copy(Δε)
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

        relative_error = norm(Dalg[active,active]-Dfd[active,active])/norm(Dfd[active,active])
        relative_asymmetry = norm(Dalg-Dalg')/norm(Dalg)
        @test relative_error < 1e-6
        @test relative_asymmetry < 1e-12
        @test iszero(Dalg[3,:])
        @test iszero(Dalg[:,3])
    end
end
