# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export LinearContact

mutable struct LinearContactState<:IpState
    ctx::Context
    σ   ::Array{Float64,1}
    w   ::Array{Float64,1}
    h   ::Float64
    function LinearContactState(ctx::Context)
        this   = new(ctx)
        this.σ = zeros(3)
        this.w = zeros(3)
        this.h = 0.0
        return this
    end
end


mutable struct LinearContact<:Material
    kn::Float64 # Normal stiffness
    ks::Float64 # Shear stiffness

    function LinearContact(; ks=NaN, kn=NaN)
        @check kn > 0.0
        @check ks >= 0.0
        return new(ks, kn)
    end
end

# Type of corresponding state structure
compat_state_type(::Type{LinearContact}, ::Type{MechInterface}, ::Context) = LinearContactState


function calcD(mat::LinearContact, state::LinearContactState)
    ndim = state.ctx.ndim
    state.w[1] > 0.0 && return zeros(ndim, ndim)

    kn = mat.kn
    ks = mat.ks

    if ndim==2
        return [  kn  0.0
                 0.0   ks ]
    else
        return  [  kn  0.0  0.0
                  0.0   ks  0.0
                  0.0  0.0   ks ]
    end
end


function update_state(mat::LinearContact, state::LinearContactState, Δu)
    ndim = state.ctx.ndim
    state.w[1:ndim] += Δu

    if state.w[1] > 0.0
        Δσ = -state.σ
    else
        D  = calcD(mat, state)
        Δσ = D*Δu
    end

    state.σ[1:ndim] += Δσ
    return Δσ, success()
end


function state_values(::LinearContact, state::LinearContactState)
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


function output_keys(::LinearContact)
    return Symbol[:jw, :jσn]
end