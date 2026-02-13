# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export LinearElastic

"""
    LinearElastic(; E=1.0, nu=0.0)

Linear elastic isotropic material model.
This model can be used with bulk elements (e.g. solids), beams, bars, and shells.
It assumes small deformations and linear stress-strain response.

# Parameters
- `E`: Young's modulus. Must be positive.
- `nu`: Poisson's ratio. Must satisfy `0 ≤ ν < 0.5`.

The shear correction factor `alpha_s` (defaulting to 5/6 for some formulations) is not a direct parameter of this constitutive model but is handled by the element formulations (e.g., `MechBeam`, `MechShell`) that use it.
"""
mutable struct LinearElastic<:Constitutive
    E ::Float64
    ν::Float64

    function LinearElastic(;
        E::Float64=NaN,
        nu::Float64=0.0,
        )
        @check E > 0.0
        @check nu >= 0.0 && nu < 0.5
        return new(E, nu)
    end
end


mutable struct LinearElasticState<:ConstState
    ctx::Context
    σ::Vec6
    ε::Vec6
    αs::Float64

    function LinearElasticState(ctx::Context)
        this = new(ctx)
        this.σ = zeros(Vec6)
        this.ε = zeros(Vec6)
        this.αs = 1.0 # Shear correction factor: to be set from shell element
        return this
    end
end


mutable struct ElasticBeamState<:ConstState
    ctx::Context
    σ::Vec3
    ε::Vec3
    αs::Float64

    function ElasticBeamState(ctx::Context)
        this = new(ctx)
        this.σ = zeros(Vec3)
        this.ε = zeros(Vec3)
        this.αs = 1.0 # Shear correction factor: to be set from beam element
        return this
    end
end


mutable struct ElasticBarState<:ConstState
    ctx::Context
    σ::Float64
    ε::Float64
    function ElasticBarState(ctx::Context; σ::Float64=0.0)
        this = new(ctx)
        this.σ = σ
        this.ε = 0.0
        return this
    end
end


mutable struct ElasticFrameState<:ConstState
    ctx::Context
    σ::Float64
    τ::Float64
    ε::Float64
    function ElasticFrameState(ctx::Context)
        this = new(ctx)
        this.σ = 0.0
        this.τ = 0.0
        this.ε = 0.0
        return this
    end
end


compat_state_type(::Type{LinearElastic}, ::Type{MechBulk})   = LinearElasticState
compat_state_type(::Type{LinearElastic}, ::Type{MechShell})  = LinearElasticState
compat_state_type(::Type{LinearElastic}, ::Type{MechBeam})   = ElasticBeamState
compat_state_type(::Type{LinearElastic}, ::Type{MechBar})    = ElasticBarState
compat_state_type(::Type{LinearElastic}, ::Type{MechEmbBar}) = ElasticBarState
compat_state_type(::Type{LinearElastic}, ::Type{MechFrame})  = ElasticFrameState


function calcDe(E::Real, ν::Real, stress_state::Symbol=:auto, αs::Float64=1.0)
    if stress_state==:plane_stress
        c = E/(1-ν^2)
        return @SArray [
            c     c*ν   0.0   0.0           0.0           0.0
            c*ν   c     0.0   0.0           0.0           0.0
            0.0   0.0   0.0   0.0           0.0           0.0
            0.0   0.0   0.0   αs*c*(1.0-ν)  0.0           0.0
            0.0   0.0   0.0   0.0           αs*c*(1.0-ν)  0.0
            0.0   0.0   0.0   0.0           0.0           c*(1.0-ν) ]
        ezz = -ν/E*(sxx+syy)
    else
        c = E/((1+ν)*(1-2*ν))
        return @SArray [
            c*(1-ν) c*ν     c*ν     0.0         0.0         0.0
            c*ν     c*(1-ν) c*ν     0.0         0.0         0.0
            c*ν     c*ν     c*(1-ν) 0.0         0.0         0.0
            0.0     0.0     0.0     c*(1-2*ν)   0.0         0.0
            0.0     0.0     0.0     0.0         c*(1-2*ν)   0.0
            0.0     0.0     0.0     0.0         0.0         c*(1-2*ν) ]
    end
end


# ❱❱❱ LinearElastic model for 3D and 2D bulk elements under plain-strain state


function calcD(mat::LinearElastic, state::LinearElasticState)
    stress_state = state.αs==1.0 ? state.ctx.stress_state : :plane_stress
    return calcDe(mat.E, mat.ν, stress_state, state.αs)
end


function update_state(mat::LinearElastic, state::LinearElasticState, cstate::LinearElasticState, Δε::AbstractArray, αs::Float64=1.0)
    De = calcD(mat, state)
    Δσ = De*Δε
    state.ε = cstate.ε + Δε
    state.σ = cstate.σ + Δσ
    return Δσ, success()
end


function state_values(mat::LinearElastic, state::LinearElasticState)
    stress_state = state.αs==1.0 ? state.ctx.stress_state : :plane_stress
    return stress_strain_dict(state.σ, state.ε, stress_state)
end


# ❱❱❱ LinearElastic for beam elements


function calcD(mat::LinearElastic, state::ElasticBeamState)
    E, ν = mat.E, mat.ν
    G    = state.αs*E/2/(1+ν)
    De = @SMatrix [ E    0.0  0.0
                    0.0  2*G  0.0
                    0.0  0.0  2*G ]
    return De

end


function update_state(mat::LinearElastic, state::ElasticBeamState, cstate::ElasticBeamState, Δε::Vector{Float64})
    D = calcD(mat, state)
    Δσ = D*Δε
    state.ε = cstate.ε + Δε
    state.σ = cstate.σ + Δσ
    return Δσ, success()
end


function state_values(mat::LinearElastic, state::ElasticBeamState)
    # vals = OrderedDict{Symbol,Float64}(
    #   :σx´   => state.σ[1],
    #   :εx´   => state.ε[1],
    #   :σx´y´ => state.σ[2]/SR2
    # )
    # σ´ = Vec6( state.σ[1], 0.0, 0.0, 0.0, 0.0, state.σ[2] )

    if state.ctx.ndim==2
        vals = OrderedDict{Symbol,Float64}(
            :σx´   => state.σ[1],
            :εx´   => state.ε[1],
            :σx´y´ => state.σ[2]/SR2
        )
        σ´ = Vec6( state.σ[1], 0.0, 0.0, 0.0, 0.0, state.σ[2] )
    else
        vals = OrderedDict{Symbol,Float64}(
            :σx´   => state.σ[1],
            :εx´   => state.ε[1],
            :σx´z´ => state.σ[2]/SR2,
            :σx´y´ => state.σ[3]/SR2
        )
        σ´ = Vec6( state.σ[1], 0.0, 0.0, 0.0, state.σ[2], state.σ[3] )
        # vals[:σx´z´] = state.σ[2]/SR2
        # vals[:σx´y´] = state.σ[3]/SR2
        # σ´[5] = state.σ[2]
        # σ´[6] = state.σ[3]
    end

    vals[:σvm] = √(3*J2(σ´))
    return vals
end


# ❱❱❱ LinearElastic model for bar elements


function calcD(mat::LinearElastic, ips::ElasticBarState)
    return mat.E
end


function update_state(mat::LinearElastic, state::ElasticBarState, cstate::ElasticBarState, Δε::Float64)
    Δσ = mat.E*Δε
    state.ε = cstate.ε + Δε
    state.σ = cstate.σ + Δσ
    return Δσ, success()
end


function state_values(mat::LinearElastic, state::ElasticBarState)
    return OrderedDict(
      :σx´ => state.σ,
      :εx´ => state.ε,
      )
end

# ❱❱❱ LinearElastic model for frame elements

function calcD(mat::LinearElastic, ips::ElasticFrameState)
    return mat.E
end


function update_state(mat::LinearElastic, state::ElasticFrameState, cstate::ElasticFrameState)
    return 0.0, success()
end


function state_values(mat::LinearElastic, state::ElasticFrameState)
    return OrderedDict(
        :σx´   => state.σ,
        :σx´y´ => state.τ,
        :εx´   => state.ε,
    )
end