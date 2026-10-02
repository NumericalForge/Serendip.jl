using Serendip
using Test
using LinearAlgebra

# ❱❱❱ Geometry and mesh

geo = GeoModel()
bl1  = add_block(geo, [0, 0], 0.1, 0.1, 0.0; nx=1, ny=1, shape=:quad4, tag="solid")
bl2  = add_block(geo, [0.1, 0], 0.1, 0.1, 0.0; nx=1, ny=1, shape=:quad4, tag="solid")
mesh = Mesh(geo)

add_cohesive_elements(mesh, tag="interface")

left_elem = select(mesh, :element, :bulk, x<=0.1)
select(get_nodes(left_elem), tag="left")

right_elem = select(mesh, :element, :bulk, x>=0.1)
select(get_nodes(right_elem), tag="right")

E    = 27.e6
nu   = 0.2
fc   = -24e3
ft   = 2.4e3
wc   = 1.7e-4
zeta = 5.0

# ❱❱❱ Finite element analyses

trajectories = [ 
    "pure extension",
    "extension with shear",
    "pure shear",
    "compression with shear" 
]

models = (
    (cmodel=MohrCoulombCohesive, props=(E=E, nu=nu, ft=ft, mu=1.4, zeta=zeta, wc=wc)),
    # (cmodel=PowerYieldCohesive, props=(E=E, nu=nu, fc=fc, ft=ft, zeta=zeta, wc=wc, alpha=1.5, gamma=0.05, theta=1.5)),
    # (cmodel=AsinhPowerCohesive, props=(E=E, nu=nu, fc=fc, ft=ft, zeta=zeta, wc=wc, alpha=0.5, beta=1.0, theta=1.0, psi=1.4, B=ft)),
    (cmodel=AsinhYieldCohesive, props=(E=E, nu=nu, fc=fc, ft=ft, zeta=zeta, wc=wc, alpha=0.33, beta=0.2, theta=1.0, psi=1.4)),
)

for trajectory in trajectories
    for model in models
        @announced_testset "$(trajectory) / $(model.cmodel)" begin
            cmodel = model.cmodel
            props = model.props

            mapper = RegionMapper()
            add_mapping(mapper, "solid", MechSolid, LinearElastic, E=E, nu=nu)
            add_mapping(mapper, "interface", MechCohesive, cmodel; props...)

            fe_model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=1.0)
            ana = MechAnalysis(fe_model)

            select(fe_model, :element, "interface", :ip, tag="jips")
            log1 = add_logger(ana, :ip, "jips", string(cmodel) * ".dat")
            add_monitor(ana, :ip, "jips", (:σn, :τ))

            stage = add_stage(ana, nincs=80, nouts=20)
            add_bc(stage, :node, "left", ux=0, uy=0)

            if trajectory == "pure extension"
                add_bc(stage, :node, "right", ux=0.0002)
            elseif trajectory == "extension with shear"
                add_bc(stage, :node, "right", ux=0.00001, uy=0.0001)
            elseif trajectory == "pure shear"
                add_bc(stage, :node, "right", ux=0.0, uy=0.001)
            elseif trajectory == "compression with shear"
                add_bc(stage, :node, "right", ux=-0.000001, uy=0.0002)
            end

            status = run(ana, autoinc=true, tol=0.1, rspan=0.03, dTmax=0.1, quiet=true)
            @test status.successful
            @test size(log1.table, 1) > 0
        end
    end
end


@announced_testset "AsinhYieldCohesive regularized algorithmic tangent" begin
    material_args = (
        E=31.0e6,
        nu=0.20,
        fc=-35.0e3,
        ft=3.0e3,
        GF=0.050,
        ft_law=:hordijk,
        alpha=0.50,
        beta=0.20,
        theta=1.00,
        psi=1.40,
        zeta=5.00,
    )
    continuum_mat = AsinhYieldCohesive(; material_args...)
    mat = AsinhYieldCohesive(; material_args..., tangent=:consistent)

    @test continuum_mat.tangent == :continuum
    @test mat.tangent == :consistent
    @test_throws ArgumentError AsinhYieldCohesive(; material_args..., tangent=:invalid)

    function initial_asinh_state()
        state = Serendip.AsinhYieldCohesiveState(Context(ndim=3))
        state.h = 0.006
        return state
    end

    function integrate_asinh(cstate, Δw)
        state = copy(cstate)
        _, status = update_state(mat, state, cstate, Δw)
        @test status.successful
        return state
    end

    function preload_asinh(target_ratio, Δw)
        state = initial_asinh_state()
        for _ in 1:400
            state = integrate_asinh(state, Δw)
            state.up >= target_ratio*mat.wc && return state
        end
        error("Could not reach up/wc=$target_ratio during preloading.")
    end

    function tangent_error(cstate, Δw, h)
        state = integrate_asinh(cstate, Δw)
        @test state.Δλ > 0.0
        Dalg = Matrix(Serendip.calcD(mat, state))
        Dfd = zeros(3, 3)

        for j in 1:3
            Δwp = copy(Δw)
            Δwm = copy(Δw)
            Δwp[j] += h
            Δwm[j] -= h

            statep = integrate_asinh(cstate, Δwp)
            statem = integrate_asinh(cstate, Δwm)
            @test statep.Δλ > 0.0
            @test statem.Δλ > 0.0
            Dfd[:, j] = (statep.σ-statem.σ)/(2h)
        end

        return (
            error=norm(Dalg-Dfd)/norm(Dfd),
            state=state,
        )
    end

    @testset "continuum tangent selection" begin
        state = integrate_asinh(initial_asinh_state(), [-2.0e-7, 2.0e-6, -1.0e-6])
        @test state.Δλ > 0.0

        kn, ks = state.kn, state.ks
        De = diagm([kn, ks, ks])
        σmax = Serendip.calc_σmax(continuum_mat, state.up)
        n, ∂f∂σmax = Serendip.yield_derivs(continuum_mat, state.σ, σmax)
        m = Serendip.potential_derivs(continuum_mat, state.σ)
        H = Serendip.calc_tangent_softening_modulus(continuum_mat, state)
        expected = De - (De*m)*(n'*De)/(dot(n, De*m) - ∂f∂σmax*H*norm(m))

        @test Matrix(Serendip.calcD(continuum_mat, state)) ≈ expected
        @test Matrix(Serendip.calcD(mat, state)) != expected
    end

    cases = (
        (
            name="fresh tension with oblique shear",
            cstate=initial_asinh_state(),
            Δw=[2.0e-7, 8.0e-7, 3.0e-7],
            h=5.0e-11,
            capped=true,
            tolerance=0.10,
        ),
        (
            name="fresh compression with sliding",
            cstate=initial_asinh_state(),
            Δw=[-2.0e-7, 2.0e-6, -1.0e-6],
            h=5.0e-11,
            capped=true,
            tolerance=0.10,
        ),
        (
            name="moderate softening with rotated shear",
            cstate=preload_asinh(0.25, [0.0, 8.0e-7, 4.0e-7]),
            Δw=[1.5e-7, -3.0e-7, 7.0e-7],
            h=2.0e-9,
            capped=false,
            tolerance=1.0e-5,
        ),
        (
            name="advanced softening in tension",
            cstate=preload_asinh(0.60, [0.0, -4.0e-7, 8.0e-7]),
            Δw=[1.5e-7, 6.0e-7, 5.0e-7],
            h=3.0e-9,
            capped=false,
            tolerance=1.0e-4,
        ),
    )

    for case in cases
        @testset "$(case.name)" begin
            result = tangent_error(case.cstate, case.Δw, case.h)
            H = Serendip.deriv_σmax_up(mat, result.state.up)
            Hcap = -mat.ft/(0.5*mat.wc)
            Hregularized = Serendip.calc_tangent_softening_modulus(mat, result.state)

            if case.capped
                @test H < Hcap
                @test Hregularized == Hcap
                @test 1.0e-2 < result.error < case.tolerance
            else
                @test H >= Hcap
                @test Hregularized == H
                @test result.error < case.tolerance
            end
        end
    end

    @testset "terminal residual stiffness" begin
        state = initial_asinh_state()
        state.w = typeof(state.w)(mat.wc, 0.0, 0.0)
        state.up = mat.wc
        state.Δλ = 1.0

        kn, ks = Serendip.calc_kn_ks(mat, state)
        De = diagm([kn, ks, ks])
        Dterminal = Matrix(Serendip.calcD(mat, state))

        @test kn ≈ mat.E*(0.1*mat.ζ)/state.h
        @test ks ≈ mat.E/(2*(1 + mat.ν))*(0.1*mat.ζ)/state.h
        @test Dterminal ≈ De*1.0e-6
    end
end
