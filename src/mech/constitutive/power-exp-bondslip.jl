# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export PowerExpBondSlip

"""
    PowerExpBondSlip(; taumax, taures=0.0, speak, sc=20*speak, alpha=0.5, ks=taumax/speak, kn=100*ks)

Construct an exponential decay bond–slip material model.

# Arguments
- `taumax`: Peak bond stress τmax (must be positive).
- `taures=0.0`: Residual bond stress τres (non-negative, ≤ τmax).
- `speak`: Slip at which τmax is reached.
- `sc=20*speak`: Critical slip where softening branch tend to stabilize.
- `alpha=0.5`: Shape parameter controlling decay (must lie in (0,1)).
- `ks=taumax/speak`: Tangential stiffness in slip direction (defaults to τmax/speak).
- `kn=100*ks`: Normal stiffness (defaults to 100·ks).

# Returns
- `PowerExpBondSlip`: A material instance representing the exponential decay bond–slip law.

# Notes
- Checks enforce positivity of parameters and ordering `τmax ≥ τres`.
- Suitable for bond–slip interface modeling with exponential decay to residual strength.
"""
mutable struct PowerExpBondSlip<:Constitutive
    τmax:: Float64
    τres:: Float64
    speak:: Float64
    sc  :: Float64
    α   :: Float64
    ks  :: Float64
    kn  :: Float64

    function PowerExpBondSlip(; 
        taumax::Real=NaN,
        taures::Real=0.0,
        speak::Real,
        sc::Real=20*speak, # critical slip (bend point)
        alpha::Real=0.5,
        ks::Real=NaN,
        kn::Real=NaN,
    )
        @check taumax > 0 "PowerExpBondSlip: τmax must be positive"
        @check taures >= 0 "PowerExpBondSlip: τres must be positive"
        ks = isnan(ks) ? taumax/speak : ks
        kn = isnan(kn) ? 100*ks : kn

        @check 0<alpha<=1 "PowerExpBondSlip: α must be in (0,1]"
        @check ks > 0 "PowerExpBondSlip: ks must be positive"
        @check kn > 0 "PowerExpBondSlip: kn must be positive"
        @check taumax >= taures "PowerExpBondSlip: τmax must be greater than τres"


        this = new(taumax, taures, speak, sc, alpha, ks, kn)
        return this
    end
end

mutable struct PowerExpBondSlipState<:IpState
    ctx::Context
    σ  ::Array{Float64,1}
    u  ::Array{Float64,1}
    τy ::Float64      # max stress
    s ::Float64      # accumulated relative displacement
    elastic::Bool
    function PowerExpBondSlipState(ctx::Context)
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


compat_state_type(::Type{PowerExpBondSlip}, ::Type{MechBondSlip}, ctx::Context) = PowerExpBondSlipState

# Type of corresponding state structure
compat_state_type(::Type{PowerExpBondSlip}) = PowerExpBondSlipState

# Element types that work with this material
compat_elem_types(::Type{PowerExpBondSlip}) = (MechBondSlip,)


function Tau(mat::PowerExpBondSlip, s::Float64)
    if s<mat.speak
        # return mat.τmax*(s/mat.speak)^mat.α
        return mat.τmax*(1 - (1 - (s/mat.speak))^(1/mat.α))
    else
        w = (s - mat.speak)/(mat.sc - mat.speak)
        z = (1 + 27*w^3)*exp(-6.93*w) - 28*w*exp(-6.93)
        return mat.τres + (mat.τmax - mat.τres)*z
    end
end


function deriv(mat::PowerExpBondSlip, state::PowerExpBondSlipState, s::Float64)
    s_factor = 0.01
    if s < s_factor*mat.speak
        s = s_factor*mat.speak   # to avoid undefined derivative
    end

    if s<=mat.speak
        # return mat.α*mat.τmax/mat.speak*(s/mat.speak)^(mat.α-1)
        return mat.τmax/(mat.speak*mat.α)*(1 - s/mat.speak)^(1/mat.α - 1)
    elseif s<mat.sc
        w = (s - mat.speak)/(mat.sc - mat.speak)
        dwds = 1/(mat.sc - mat.speak)
        dzdw = ( 81*w^2 - 6.93*(1 + 27*w^3))*exp(-6.93*w) - 28*exp(-6.93)
        return (mat.τmax - mat.τres)*dzdw*dwds
    else
        return mat.ks*1e-3
    end
end


function calcD(mat::PowerExpBondSlip, state::PowerExpBondSlipState)
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


function yield_func(mat::PowerExpBondSlip, state::PowerExpBondSlipState, τ::Float64)
    return abs(τ) - state.τy
end

function update_state(mat::PowerExpBondSlip, state::PowerExpBondSlipState, Δu::Vect)
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

        if state.s<mat.speak && Δs>0.2*mat.speak
            return state.σ, failure("PowerExpBondSlip: Plastic slip is too large")
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


function state_values(mat::PowerExpBondSlip, state::PowerExpBondSlipState)
    return OrderedDict(
      :s => state.u[1] ,
      :τ => state.σ[1] ,
      )
end

