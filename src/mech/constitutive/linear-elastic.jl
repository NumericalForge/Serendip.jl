# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export LinearElastic

"""
    LinearElastic(; E=1.0, nu=0.0, alpha_s=5/6)

Linear elastic isotropic material model.
This model can be used with bulk elements (e.g. solids), beams, bars, and shells.
It assumes small deformations and linear stress-strain response.

# Parameters
- `E`: *Young's modulus*. Must be positive.
- `nu = 0.0`: *Poisson's ratio*. Must satisfy `0 ≤ ν < 0.5`.
- `alpha_s = 5/6`: *Shear correction factor*, used in beam/shell formulations. Must be positive.
"""
mutable struct LinearElastic<:Constitutive
    E ::Float64
    ν::Float64
    αs::Float64

    function LinearElastic(;
        E::Float64=NaN,
        nu::Float64=0.0,
        alpha_s::Float64=5/6
        )
        @check E > 0.0
        @check nu >= 0.0 && nu < 0.5
        @check alpha_s > 0.0
        return new(E, nu, alpha_s)
    end
end


mutable struct ElasticSolidState<:IpState
    ctx::Context
    σ::Vec6
    ε::Vec6

    function ElasticSolidState(ctx::Context)
        this = new(ctx)
        this.σ = zeros(Vec6)
        this.ε = zeros(Vec6)
        return this
    end
end


mutable struct ElasticPlaneStressState<:IpState
    ctx::Context
    σ::Vec6
    ε::Vec6

    function ElasticPlaneStressState(ctx::Context)
        this = new(ctx)
        this.σ = zeros(Vec6)
        this.ε = zeros(Vec6)
        return this
    end
end


mutable struct ElasticBeamState<:IpState
    ctx::Context
    σ::Vec3
    ε::Vec3

    function ElasticBeamState(ctx::Context)
        this = new(ctx)
        this.σ = zeros(Vec3)
        this.ε = zeros(Vec3)
        return this
    end
end

mutable struct ElasticBarState<:IpState
    ctx::Context
    σ::Float64
    ε::Float64
    function ElasticBarState(ctx::Context)
        this = new(ctx)
        this.σ = 0.0
        this.ε = 0.0
        return this
    end
end

mutable struct ElasticFrameState<:IpState
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


compat_state_type(::Type{LinearElastic}, ::Type{MechBulk}, ctx::Context)  = ctx.stress_state==:plane_stress ? ElasticPlaneStressState : ElasticSolidState
compat_state_type(::Type{LinearElastic}, ::Type{MechShell}, ctx::Context)  = ElasticPlaneStressState
compat_state_type(::Type{LinearElastic}, ::Type{MechBeam}, ctx::Context)   = ElasticBeamState
compat_state_type(::Type{LinearElastic}, ::Type{MechBar}, ctx::Context)    = ElasticBarState
compat_state_type(::Type{LinearElastic}, ::Type{MechEmbBar}, ctx::Context) = ElasticBarState
compat_state_type(::Type{LinearElastic}, ::Type{MechFrame}, ctx::Context) = ElasticFrameState


function calcDe(E::Real, ν::Real, stress_state::Symbol=:d3, αs::Float64=1.0)
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


# LinearElastic model for 3D and 2D bulk elements under plain-strain state


function calcD(mat::LinearElastic, state::ElasticSolidState)
    return calcDe(mat.E, mat.ν)
end


function update_state(mat::LinearElastic, state::ElasticSolidState, dε::AbstractArray)
    De = calcDe(mat.E, mat.ν)

    dσ = De*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end


function state_values(mat::LinearElastic, state::ElasticSolidState)
    return stress_strain_dict(state.σ, state.ε, state.ctx.stress_state)
end


# LinearElastic model for 2D bulk elements under plane-stress state and shell elements


function calcD(mat::LinearElastic, state::ElasticPlaneStressState; is_shell::Bool=false)
    αs = is_shell ? mat.αs : 1.0
    return calcDe(mat.E, mat.ν, :plane_stress, αs)
end


function update_state(mat::LinearElastic, state::ElasticPlaneStressState, dε::AbstractArray; is_shell::Bool=false)
    αs = is_shell ? mat.αs : 1.0
    De = calcDe(mat.E, mat.ν, :plane_stress, αs)
    dσ = De*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end


function state_values(mat::LinearElastic, state::ElasticPlaneStressState)
    return stress_strain_dict(state.σ, state.ε, :plane_stress)
end


# LinearElastic for beam elements



function calcD(mat::LinearElastic, state::ElasticBeamState)
    E, ν = mat.E, mat.ν

    αs   = mat.αs
    De = @SMatrix [ E    0.0         0.0
                    0.0  αs*E/(1+ν)  0.0
                    0.0  0.0         αs*E/(1+ν) ]
    return De

    # return @SMatrix [
        # E  0.0  0.0
        # 0.0  E/(1+ν)  0.0
        # 0.0  0.0  E/(1+ν)
    # ]
end


function update_state(mat::LinearElastic, state::ElasticBeamState, dε::Vector{Float64})
    D = calcD(mat, state)
    dσ = D*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
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


# LinearElastic model for bar elements

function calcD(mat::LinearElastic, ips::ElasticBarState)
    return mat.E
end


function update_state(mat::LinearElastic, state::ElasticBarState, Δε::Float64)
    Δσ = mat.E*Δε
    state.ε += Δε
    state.σ += Δσ
    return Δσ, success()
end


function state_values(mat::LinearElastic, state::ElasticBarState)
    return OrderedDict(
      :σx´ => state.σ,
      :εx´ => state.ε,
      )
end

# LinearElastic model for frame elements

function calcD(mat::LinearElastic, ips::ElasticFrameState)
    return mat.E
end


function update_state(mat::LinearElastic, state::ElasticFrameState)
    return 0.0, success()
end


function state_values(mat::LinearElastic, state::ElasticFrameState)
    return OrderedDict(
        :σx´   => state.σ,
        :σx´y´ => state.τ,
        :εx´   => state.ε,
    )
end