# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export ExpDecayBondSlip

mutable struct ExpoDecayBondSlipState<:IpState
    ctx::Context
    σ  ::Array{Float64,1}
    u  ::Array{Float64,1}
    τy ::Float64      # max stress
    s ::Float64      # accumulated relative displacement
    elastic::Bool
    function ExpoDecayBondSlipState(ctx::Context)
        this = new(ctx)
        ndim = ctx.ndim
        this.σ = zeros(ndim)
        this.u = zeros(ndim)
        this.τy = 0.0
        this.s = 0.0
        this.elastic = false
        return this
    end
end


# CebBondSlip_params = [
#     FunInfo(:ExpDecayBondSlip, "Consitutive model for a rod-solid interface according to CEB."),
#     KwArgInfo(:taumax, "Shear strength", cond=:(taumax>0)),
#     KwArgInfo(:taures, "Residual shear stress", cond=:(taures>=0)),
#     KwArgInfo(:speak, "Slip for the peak value", cond=:(speak>0)),
#     KwArgInfo(:sc, "Characteristic slip where τ remain constant", cond=:(sc>0)),
#     KwArgInfo(:alpha, "Ascending curvature parameter", 0.4, cond=:(0.0<=alpha<=1.0)),
#     # KwArgInfo(:beta, "Descending curvature parameter", 1.0, cond=:(0.0<=beta<=1.0)),
#     KwArgInfo(:kn, "Normal stiffness", cond=:(kn>0)),
#     KwArgInfo(:ks, "Shear stiffness", nothing),
# ]
# @doc docstring(CebBondSlip_params) ExpDecayBondSlip(; kwargs...)


mutable struct ExpDecayBondSlip<:Material
    τmax:: Float64
    τres:: Float64
    speak:: Float64
    sc  :: Float64
    α   :: Float64
    # β   :: Float64
    ks  :: Float64
    kn  :: Float64

    function ExpDecayBondSlip(; kwargs...)
        args = checkargs(kwargs, CebBondSlip_params)
        ks = args.ks === nothing ? args.taumax/args.speak : args.ks
        @check args.taumax > args.taures

        this = new(args.taumax, args.taures, args.speak, args.sc, args.alpha, ks, args.kn)
        return this
    end
end


compat_state_type(::Type{ExpDecayBondSlip}, ::Type{MechBondSlip}, ctx::Context) = ExpoDecayBondSlipState

# Type of corresponding state structure
compat_state_type(::Type{ExpDecayBondSlip}) = ExpoDecayBondSlipState

# Element types that work with this material
compat_elem_types(::Type{ExpDecayBondSlip}) = (MechBondSlip,)


function Tau(mat::ExpDecayBondSlip, s::Float64)
    if s<mat.speak
        return mat.τmax*(s/mat.speak)^mat.α
    else
        w = (s - mat.speak)/(mat.sc - mat.speak)
        z = (1 + 27*w^3)*exp(-6.93*w) - 28*w*exp(-6.93)
        return mat.τres + (mat.τmax - mat.τres)*z
    end
end


function deriv(mat::ExpDecayBondSlip, state::ExpoDecayBondSlipState, s::Float64)
    s_factor = 0.01
    if s < s_factor*mat.speak
        s = s_factor*mat.speak   # to avoid undefined derivative
    end

    if s<=mat.speak
        return mat.α*mat.τmax/mat.speak*(s/mat.speak)^(mat.α-1)
    elseif s<mat.sc
        w = (s - mat.speak)/(mat.sc - mat.speak)
        dwds = 1/(mat.sc - mat.speak)
        dzdw = ( 81*w^2 - 6.93*(1 + 27*w^3))*exp(-6.93*w) - 28*exp(-6.93)
        return (mat.τmax - mat.τres)*dzdw*dwds
    else
        return mat.ks*1e-3
    end
end


function calcD(mat::ExpDecayBondSlip, state::ExpoDecayBondSlipState)
    ks = mat.ks

    if !state.elastic
        dτydsy = deriv(mat, state, state.s)
        ks = dτydsy
    end

    kn = mat.kn
    if state.ctx.ndim==2
        return [  ks  0.0
                 0.0   kn ]
    else
        return  [  ks  0.0  0.0
                  0.0   kn  0.0
                  0.0  0.0   kn ]
    end
end


function yield_func(mat::ExpDecayBondSlip, state::ExpoDecayBondSlipState, τ::Float64)
    return abs(τ) - state.τy
end

function update_state(mat::ExpDecayBondSlip, state::ExpoDecayBondSlipState, Δu::Vect)
    ks = mat.ks
    kn = mat.kn
    Δs = Δu[1]      # relative displacement
    τini = state.σ[1] # initial shear stress
    τtr  = τini + ks*Δs # elastic trial

    ftr  = yield_func(mat, state, τtr)

    if ftr<0.0
        τ = τtr
        state.elastic = true
    else
        Δs = (abs(τtr)-state.τy)/ks # only plastic part

        # if state.elastic
        #     @show Δs
        # end

        if state.s<mat.speak && Δs>0.2*mat.speak
            # @show Δs
            return state.σ, failure("ExpDecayBondSlip: Plastic slip is too large")
        end

        state.s  += Δs
        state.τy  = Tau(mat, state.s)
        τ         = state.τy*sign(τtr)

        state.elastic  = false
    end

    # calculate Δσ
    Δτ = τ - τini
    Δσ = kn*Δu
    Δσ[1] = Δτ

    # update u and σ
    state.u .+= Δu
    state.σ .+= Δσ

    return Δσ, success()
end


function state_values(mat::ExpDecayBondSlip, state::ExpoDecayBondSlipState)
    return OrderedDict(
      :s   => state.u[1] ,
      :tau  => state.σ[1] ,
      )
end

