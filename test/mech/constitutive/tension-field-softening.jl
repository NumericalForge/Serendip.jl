using Serendip
using Test
using LinearAlgebra
using Random

@testset "TensionFieldSoftening parameters and initialization" begin
    @test_throws ArgumentError TensionFieldSoftening(E=-1, fy=2)
    @test_throws ArgumentError TensionFieldSoftening(E=200, fy=2, nu=0.5)
    @test_throws ArgumentError TensionFieldSoftening(E=200, fy=2, GF=0)
    @test_throws ArgumentError TensionFieldSoftening(E=200, fy=2, beta=0.5)
    mat = TensionFieldSoftening(E=200, nu=0.3, fy=2, GF=0.2)
    state = Serendip.TensionFieldSofteningState(Context(ndim=3))
    @test norm(Serendip.calcD(mat, state)) > 0
    _, status = update_state(mat, state, copy(state), [0.01, 0., 0.])
    @test !status.successful
    @test state.branch == :initial
    @test_throws DimensionMismatch update_state(mat,state,copy(state),[0.0,0.0])
    _,status=update_state(mat,state,copy(state),[NaN,0.0,0.0])
    @test !status.successful && state.branch == :initial
end

@testset "TensionFieldSoftening states and exact tangents" begin
    rng = MersenneTwister(92)
    counts = Dict{Symbol,Int}()
    for beta in (0., 1., 2.)
        mat = beta == 0 ? TensionFieldSoftening(E=200, nu=0.3, fy=2) :
            TensionFieldSoftening(E=200, nu=0.3, fy=2, GF=0.4/(beta+1), beta=beta)
        cs = Serendip.TensionFieldSofteningState(Context(ndim=3))
        cs.h = 1.0
        cs.εpa = 0.015
        # Include a nonzero, noncoaxial committed plastic tensor.
        cs.εp = Serendip.Vec3(0.003, -0.001, 0.0007)
        for i in 1:100
            angle = 6rand(rng)-3
            c, s = cos(angle), sin(angle)
            V = [c -s; s c]
            e = V*Diagonal(0.095rand(rng,2).-0.03)*V'
            strain = [e[1,1], e[2,2], √2*e[1,2]] + cs.εp
            state = copy(cs)
            _, status = update_state(mat, state, cs, strain)
            @test status.successful
            status.successful || continue
            counts[state.branch] = get(counts, state.branch, 0)+1
            st = [state.σ[1] state.σ[3]/√2; state.σ[3]/√2 state.σ[2]]
            wt = [state.ω[1] state.ω[3]/√2; state.ω[3]/√2 state.ω[2]]
            @test minimum(eigvals(Symmetric(st))) >= -1e-9
            @test minimum(eigvals(Symmetric(wt))) >= -1e-9
            @test abs(dot(state.σ,state.ω)) < 1e-8
            @test state.σ ≈ Serendip.membrane_elastic_matrix(200.,0.3)*(state.ε-state.εp+state.ω) atol=1e-8
            @test state.εpa >= cs.εpa
            Dfd = zeros(3,3)
            dh = 1e-7
            for j in 1:3
                step = zeros(3); step[j] = dh
                sp, sm = copy(cs), copy(cs)
                _, okp = update_state(mat, sp, cs, strain+step)
                _, okm = update_state(mat, sm, cs, strain-step)
                @test okp.successful && okm.successful
                Dfd[:,j] = (sp.σ-sm.σ)/(2dh)
            end
            D = Serendip.calcD(mat,state)
            @test norm(D-Dfd)/max(1,norm(Dfd)) < 2e-5
            @test norm(D-D')/max(1,norm(D)) < 1e-10
        end
        for ex in (0.001, 0.03)
            strain = [ex,ex,0.]+cs.εp
            state = copy(cs)
            _, status = update_state(mat,state,cs,strain)
            @test status.successful
            dh = 1e-7
            Dfd = zeros(3,3)
            for j in 1:3
                step = zeros(3);step[j]=dh
                sp,sm=copy(cs),copy(cs)
                update_state(mat,sp,cs,strain+step)
                update_state(mat,sm,cs,strain-step)
                Dfd[:,j]=(sp.σ-sm.σ)/(2dh)
            end
            @test norm(Serendip.calcD(mat,state)-Dfd)/norm(Dfd) < 2e-5
        end
    end
    @test all(get(counts,k,0)>0 for k in (:slack,:wrinkled_elastic,:wrinkled_plastic,:taut_elastic,:taut_plastic))
end

@testset "TensionFieldSoftening loading cycles, energy and rupture" begin
    for beta in (1.,2.), h in (0.5,1.,2.)
        GF = 0.2
        mat = TensionFieldSoftening(E=200,nu=0.3,fy=2,GF=GF,beta=beta)
        cs = Serendip.TensionFieldSofteningState(Context(ndim=3));cs.h=h
        crit = (beta+1)*GF/(h*mat.fy)
        # Independent integration of the prescribed plastic-work envelope.
        n=1000
        energy = h*sum((Serendip.tension_field_strength(mat,crit*i/n,crit)[1]+
            Serendip.tension_field_strength(mat,crit*(i+1)/n,crit)[1])*crit/(2n) for i in 0:n-1)
        @test energy ≈ GF rtol=1e-6
        # Both a uniaxial wrinkled path and a biaxial taut path reach rupture.
        for biaxial in (false,true)
            committed=copy(cs)
            for ex in range(0.001,(biaxial ? 0.4 : 0.9)*crit,length=30)
                ey = biaxial ? ex : -2ex
                state=copy(committed)
                _,status=update_state(mat,state,committed,[ex,ey,0.]-committed.ε)
                @test status.successful
                status.successful || break
                @test !state.ruptured
                Serendip.commit_state(committed,state)
            end
            # Discard a failed trial to verify that rollback cannot heal or
            # accidentally commit a material point that has not yet failed.
            state=copy(committed)
            _,status=update_state(mat,state,committed,[2crit,biaxial ? 2crit : -4crit,0.]-committed.ε)
            @test status.successful && state.ruptured
            @test !committed.ruptured
            @test iszero(state.σ) && iszero(Serendip.calcD(mat,state))
            @test state.εpa ≈ crit
            Serendip.reset_state(state,committed)
            @test !state.ruptured
            update_state(mat,state,committed,[2crit,biaxial ? 2crit : -4crit,0.]-committed.ε)
            Serendip.commit_state(committed,state)
            terminal=committed.εp
            for ex in (0.,-crit,3crit)
                _,status=update_state(mat,state,committed,[ex,ex,0.]-committed.ε)
                @test status.successful && state.ruptured
                @test iszero(state.σ) && iszero(Serendip.calcD(mat,state))
                @test state.εp == terminal
            end
        end
        # Elastic unloading, slackening and reloading before rupture.
        committed=copy(cs)
        for (k,ex) in enumerate((0.04,0.035,-0.02,0.,0.04,0.06))
            state=copy(committed)
            oldp=committed.εpa
            _,status=update_state(mat,state,committed,[ex,-0.1,0.]-committed.ε)
            @test status.successful
            k > 1 && ex in (0.035,-0.02,0.,0.04) && @test state.εpa ≈ oldp
            Serendip.commit_state(committed,state)
        end
    end
    # This softening law snaps back before a continuous strain-controlled
    # rupture can be followed. Reject the update without marking rupture.
    mat=TensionFieldSoftening(E=200,nu=0.3,fy=2,GF=0.001,beta=1)
    cs=Serendip.TensionFieldSofteningState(Context(ndim=3));cs.h=1
    state=copy(cs)
    _,status=update_state(mat,state,cs,[0.02,-0.04,0.])
    @test !status.successful
    @test !state.ruptured && !cs.ruptured
end

@testset "TensionFieldSoftening one-sided transitions" begin
    for beta in (0.,1.,2.)
        mat = beta==0 ? TensionFieldSoftening(E=200,nu=0.3,fy=2) :
            TensionFieldSoftening(E=200,nu=0.3,fy=2,GF=0.4/(beta+1),beta=beta)
        cs=Serendip.TensionFieldSofteningState(Context(ndim=3));cs.h=1.
        dh=1e-8
        for (axial,left,right) in ((0.,0.,200.),(0.01,200.,beta==0 ? 0. : 200*(-10beta)/(200-10beta)))
            center=copy(cs);minus=copy(cs);plus=copy(cs)
            _,ok=update_state(mat,center,cs,[axial,-0.1,0.])
            _,okm=update_state(mat,minus,cs,[axial-dh,-0.1,0.])
            _,okp=update_state(mat,plus,cs,[axial+dh,-0.1,0.])
            @test ok.successful && okm.successful && okp.successful
            @test (center.σ[1]-minus.σ[1])/dh ≈ left atol=1e-5
            @test (plus.σ[1]-center.σ[1])/dh ≈ right atol=1e-5
        end
        center=copy(cs)
        update_state(mat,center,cs,[0.002,-0.0006,0.])
        for sign in (-1,1)
            side=copy(cs)
            _,status=update_state(mat,side,cs,[0.002,-0.0006+sign*dh,0.])
            @test status.successful
            @test (side.σ-center.σ)/(sign*dh) ≈ Serendip.calcD(mat,side)[:,2] atol=1e-6
        end
    end
end
