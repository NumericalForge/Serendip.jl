# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export LinearCohesive

mutable struct LinearCohesiveState<:IpState
    ctx::Context
    σ  ::Array{Float64,1}
    w  ::Array{Float64,1}
    h  ::Float64
    function LinearCohesiveState(ctx::Context)
        this = new(ctx)
        this.σ = zeros(ctx.ndim)
        this.w = zeros(ctx.ndim)
        return this
    end
end


mutable struct LinearCohesive<:Constitutive
    kn::Float64 # Normal stiffness
    ks::Float64 # Shear stiffness
    zeta::Float64 # Thickness ratio

    function LinearCohesive(; E::Real=NaN, nu::Real=0.0, zeta::Real=5.0)
        @check E > 0.0 "LinearCohesive: Young's modulus E must be positive"
        @check nu >= 0.0 "LinearCohesive: Poisson's ratio nu must be non-negative"
        @check zeta > 0.0 "LinearCohesive: Thickness ratio zeta must be positive"
        return new(E, nu, zeta)
    end
end


# Type of corresponding state structure
compat_state_type(::Type{<:LinearCohesive}, ::Type{MechInterface}, ::Context) = LinearCohesiveState

# LinearCohesive
function calcD(mat::LinearCohesive, state::LinearCohesiveState)
    ndim = state.ctx.ndim
    kn = mat.E*mat.ζ/state.h
    G  = mat.E/(2*(1 + mat.ν))
    ks = G*mat.ζ/state.h

    if ndim==2
        return [  kn  0.0
                 0.0   ks ]
    else
        return  [  kn  0.0  0.0
                  0.0   ks  0.0
                  0.0  0.0   ks ]
    end
end


function update_state(mat::LinearCohesive, state::LinearCohesiveState, Δu)
    D  = calcD(mat, state)
    Δσ = D*Δu

    state.w += Δu
    state.σ += Δσ
    return Δσ, success()
end


function state_values(::LinearCohesive, state::LinearCohesiveState)
    ndim = state.ctx.ndim
    if ndim == 3
       return Dict(
          :jw  => state.w[1],
          :jw2  => state.w[2],
          :jw3  => state.w[3],
          :jσn  => state.σ[1],
          :js2  => state.σ[2],
          :js3  => state.σ[3],
          )
    else
        return Dict(
          :jw  => state.w[1],
          :jw2  => state.w[2],
          :jσn  => state.σ[1],
          :js2  => state.σ[2],
          )
    end
end


function output_keys(::LinearCohesive)
    return Symbol[:jw, :jw]
end