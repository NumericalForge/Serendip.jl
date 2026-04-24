# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export LinearCohesive

mutable struct LinearCohesiveState<:ConstState
    ctx::Context
    σ  ::Vec3
    w  ::Vec3
    h  ::Float64
    function LinearCohesiveState(ctx::Context)
        this = new(ctx)
        this.σ = zeros(Vec3)
        this.w = zeros(Vec3)
        return this
    end
end


mutable struct LinearCohesive<:Constitutive
    E::Float64
    ν::Float64
    ζ::Float64

    function LinearCohesive(; E::Real=NaN, nu::Real=0.0, zeta::Real=5.0)
        @check E > 0.0 "LinearCohesive: Young's modulus E must be positive"
        @check nu >= 0.0 "LinearCohesive: Poisson's ratio nu must be non-negative"
        @check zeta > 0.0 "LinearCohesive: Thickness ratio zeta must be positive"
        return new(E, nu, zeta)
    end
end


# Type of corresponding state structure
compat_state_type(::Type{<:LinearCohesive}, ::Type{MechCohesive}) = LinearCohesiveState

# LinearCohesive
function calcD(mat::LinearCohesive, state::LinearCohesiveState)
    kn = mat.E*mat.ζ/state.h
    G  = mat.E/(2*(1 + mat.ν))
    ks = G*mat.ζ/state.h

    return @SMatrix [  kn  0.0  0.0
                      0.0   ks  0.0
                      0.0  0.0   ks ]
end


function update_state(mat::LinearCohesive, state::LinearCohesiveState, cstate::LinearCohesiveState, Δu::AbstractArray)
    Δw = length(Δu) == 2 ? Vec3(Δu[1], Δu[2], 0.0) : Vec3(Δu[1], Δu[2], Δu[3])
    D  = calcD(mat, cstate)
    Δσ = D*Δw

    state.w = cstate.w + Δw
    state.σ = cstate.σ + Δσ
    return Δσ, success()
end


function state_values(::LinearCohesive, state::LinearCohesiveState)
    ndim = state.ctx.ndim
    σn, τ1, τ2 = state.σ
    τ = √(τ1^2 + τ2^2)
    if ndim == 3
        return Dict(
            :w => state.w[1],
            :σn => σn,
            :τ  => τ,
            :σ2 => τ1,
            :σ3 => τ2,
          )
    else
        return Dict(
            :w => state.w[1],
            :σn => σn,
            :τ  => τ,
            :σ2 => τ1,
        )
    end
end


function output_keys(::LinearCohesive)
    return Symbol[:w, :σn, :τ]
end
