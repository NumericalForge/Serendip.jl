# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export LinearInterface

mutable struct LinearInterfaceState<:IpState
    ctx::Context
    σ   ::Array{Float64,1}
    w   ::Array{Float64,1}
    h   ::Float64
    function LinearInterfaceState(ctx::Context)
        this = new(ctx)
        this.σ = zeros(ctx.ndim)
        this.w = zeros(ctx.ndim)
        this.h = 0.0
        return this
    end
end


mutable struct LinearInterface<:Constitutive
    kn::Float64 # Normal stiffness
    ks::Float64 # Shear stiffness

    function LinearInterface(; ks::Float64=NaN, kn::Float64=NaN)
        @check kn > 0.0
        @check ks >= 0.0
        return new(ks, kn)
    end
end


# Type of corresponding state structure
compat_state_type(::Type{LinearInterface}, ::Type{MechInterface}, ::Context) = LinearInterfaceState


function calcD(mat::LinearInterface, state::LinearInterfaceState)
    ndim = state.ctx.ndim
    kn   = mat.kn
    ks   = mat.ks

    if ndim==2
        return [  kn  0.0
                 0.0   ks ]
    else
        return  [  kn  0.0  0.0
                  0.0   ks  0.0
                  0.0  0.0   ks ]
    end
end


function update_state(mat::LinearInterface, state::LinearInterfaceState, Δu)
    ndim = state.ctx.ndim
    D  = calcD(mat, state)
    Δσ = D*Δu

    state.w += Δu
    state.σ += Δσ
    return Δσ, success()
end


function state_values(::LinearInterface, state::LinearInterfaceState)
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


function output_keys(::LinearInterface)
    return Symbol[:jw, :jw]
end