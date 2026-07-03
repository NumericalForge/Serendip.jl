# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export PowerExpBondSlip

"""
    PowerExpBondSlip(; taumax, taures=0.0, speak, sc=20*speak, alpha=0.5, beta=20.0,
                     eta=0.0, mu=0.0, eta_s=0.0, ks=10*taumax/speak, kn=100*ks)

Construct an exponential decay bond–slip material model.

# Arguments
- `taumax`: Base peak bond stress at near-zero confinement (must be positive).
- `taures=0.0`: Base residual bond stress at near-zero confinement
  (non-negative, ≤ `taumax`).
- `speak`: Slip at which τmax is reached.
- `sc=20*speak`: Base critical slip at near-zero confinement where the softening
  branch tends to stabilize.
- `alpha=0.5`: Shape parameter controlling decay (must lie in `(0,1]`).
- `beta=20.0`: Shape parameter for the post-peak transition (`beta > 0`).
- `eta=0.0`: Confinement sensitivity of the peak bond stress,
  `τmax_eff = taumax - η p`.
- `mu=0.0`: Confinement sensitivity of the residual bond stress,
  `τres_eff = taures - μ p`.
- `eta_s=0.0`: Compliance-like confinement sensitivity of the residual slip,
  `sc_eff = sc - ηs p`.
- `ks=10*taumax/speak`: Tangential stiffness in slip direction.
- `kn=100*ks`: Normal stiffness (defaults to 100·ks).

# Returns
- `PowerExpBondSlip`: A material instance representing the exponential decay bond–slip law.

# Notes
- Checks enforce positivity of parameters and ordering `τmax ≥ τres`.
- Suitable for bond–slip interface modeling with exponential decay to residual strength.
- A practical calibration estimate is `speak ≈ 0.06ϕ`, where `ϕ` is the bar diameter.
- Compression is assumed negative, so `p < 0` increases `τmax`, `τres`, and `sc`
  when `eta`, `mu`, and `eta_s` are positive.
- For confinement-aware analyses, the constructor inputs `taumax`, `taures`, and `sc`
  are interpreted as base values at near-zero confinement pressure.
- In that confinement-aware setting, a practical starting point is
  `taumax ≈ 7.5 MPa` for rebars up to about `16 mm` diameter,
  `taures ≈ 0.07*taumax`, and `sc ≈ 1.6*lr` or `sc ≈ 10 mm` when
  rib spacing data are not available.
- The effective confined values used internally are `τmax_eff = taumax - η p`,
  `τres_eff = taures - μ p`, and `sc_eff = sc - ηs p`.
"""
mutable struct PowerExpBondSlip<:Constitutive
    τmax:: Float64
    τres:: Float64
    speak:: Float64
    sc  :: Float64
    α   :: Float64
    β   :: Float64
    η   :: Float64
    μ   :: Float64
    ηs  :: Float64
    ks  :: Float64
    kn  :: Float64

    function PowerExpBondSlip(; 
        taumax::Real=NaN,
        taures::Real=0.0,
        speak::Real,
        sc::Real=20*speak, # critical slip (bend point)
        alpha::Real=0.5,
        beta::Real=20.0,
        eta::Real=0.0,
        mu::Real=0.0,
        eta_s::Real=0.0,
        ks::Real=NaN,
        kn::Real=NaN,
    )
        @check taumax > 0 "PowerExpBondSlip: τmax must be positive"
        @check taures >= 0 "PowerExpBondSlip: τres must be non-negative"
        @check speak > 0 "PowerExpBondSlip: speak must be positive"
        @check sc > speak "PowerExpBondSlip: sc must be greater than speak"
        ks = isnan(ks) ? 10*taumax/speak : ks
        kn = isnan(kn) ? 100*ks : kn

        @check 0<alpha<=1 "PowerExpBondSlip: α must be in (0,1]"
        @check 0<beta "PowerExpBondSlip: β must be positive"
        @check eta >= 0 "PowerExpBondSlip: η must be non-negative"
        @check mu >= 0 "PowerExpBondSlip: μ must be non-negative"
        @check eta_s >= 0 "PowerExpBondSlip: ηs must be non-negative"
        @check ks > 0 "PowerExpBondSlip: ks must be positive"
        @check kn > 0 "PowerExpBondSlip: kn must be positive"
        @check taumax >= taures "PowerExpBondSlip: τmax must be greater than τres"


        this = new(taumax, taures, speak, sc, alpha, beta, eta, mu, eta_s, ks, kn)
        return this
    end
end


mutable struct PowerExpBondSlipState<:ConstState
    ctx::Context
    σ  ::Vector{Float64}
    u  ::Vector{Float64}
    τy ::Float64      # max stress
    s ::Float64      # accumulated relative displacement
    p ::Float64      # confinement pressure
    elastic::Bool
    function PowerExpBondSlipState(ctx::Context; p::Real=NaN)
        this = new(ctx)
        ndim = ctx.ndim
        this.σ = zeros(ndim)
        this.u = zeros(ndim)
        this.τy = 0.0
        this.s = 0.0
        this.p = p
        this.elastic = false
        return this
    end
end


compat_state_type(::Type{PowerExpBondSlip}, ::Type{MechBondSlip}) = PowerExpBondSlipState

# Type of corresponding state structure
compat_state_type(::Type{PowerExpBondSlip}) = PowerExpBondSlipState

# Element types that work with this material
compat_elem_types(::Type{PowerExpBondSlip}) = (MechBondSlip,)


function effective_params(mat::PowerExpBondSlip, p::Float64)
    τmax = mat.τmax - mat.η*p
    τres = mat.τres - mat.μ*p
    sc   = mat.sc   - mat.ηs*p

    τmax = max(τmax, eps(Float64))
    τres = clamp(τres, 0.0, τmax)
    sc   = max(sc, mat.speak*(1.0 + 1.0e-6))

    return τmax, τres, sc
end


@inline function confinement_pressure(state::PowerExpBondSlipState)
    return isnan(state.p) ? 0.0 : state.p
end


function Tau(mat::PowerExpBondSlip, s::Float64, p::Float64=0.0)
    τmax, τres, sc = effective_params(mat, p)

    if s<mat.speak
        return τmax*(1 - ((s-mat.speak)/mat.speak)^2)^mat.α
    elseif s<=sc
        β=mat.β

        w = (s - mat.speak)/(sc - mat.speak)
        z = (2*(1-w)*(1 + β*w)) / (1 + (1 + β*w)^2)
        return τres + (τmax - τres)*z
    else
        return τres
    end
end


function deriv(mat::PowerExpBondSlip, state::PowerExpBondSlipState, s::Float64)
    p = confinement_pressure(state)
    τmax, τres, sc = effective_params(mat, p)
    s_factor = 0.01
    if s < s_factor*mat.speak
        s = s_factor*mat.speak   # to avoid undefined derivative
    end

    if s<=mat.speak
        s̄ = s/mat.speak
        g = max(2*s̄ - s̄^2, eps(Float64))
        return 2*mat.α*τmax/mat.speak*(1 - s̄)*g^(mat.α - 1)
    elseif s<=sc
        β=mat.β
        w = (s - mat.speak)/(sc - mat.speak)
        num = -2*( β^3*w^2 + β^2*w^2 + 2*β^2*w + 4*β*w + 2 )
        den = (1 + (1+β*w)^2)^2
        dzdw = num/den
        return (τmax - τres)/(sc - mat.speak)*dzdw
    else
        return mat.ks*1e-4
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


function yield_func(mat::PowerExpBondSlip, state::PowerExpBondSlipState, τ::Float64, τy::Float64)
    return abs(τ) - τy
end


function plastic_update(mat::PowerExpBondSlip, state::PowerExpBondSlipState, cstate::PowerExpBondSlipState, τtr::Float64, τyini::Float64)
    ks = mat.ks
    p  = confinement_pressure(state)
    _, _, sc = effective_params(mat, p)

    Δs = (abs(τtr) - τyini)/ks

    if cstate.s < mat.speak && Δs > 0.2*mat.speak
        return state.σ, failure("PowerExpBondSlip: Plastic slip is too large")
    end

    state.s  = cstate.s + Δs
    if state.s < mat.speak
        state.τy = Tau(mat, state.s, p)
    else
        # Clamp the history variable to the active softening range only for the envelope evaluation.
        state.τy = Tau(mat, min(state.s, sc), p)
    end
    state.elastic = false

    return state.τy*sign(τtr), success()
end


function update_state(mat::PowerExpBondSlip, state::PowerExpBondSlipState, cstate::PowerExpBondSlipState, Δu::Vect)
    ks = mat.ks
    kn = mat.kn
    p = confinement_pressure(state)
    Δs = Δu[1]         # relative displacement
    τini = cstate.σ[1] # initial bond stress
    τtr  = τini + ks*Δs # elastic trial
    state.p = p
    τyini = Tau(mat, cstate.s, p)

    ftr  = yield_func(mat, state, τtr, τyini)

    if ftr<0.0
        τ = τtr
        state.s = cstate.s
        state.τy = τyini
        state.elastic = true
    else
        τ, status = plastic_update(mat, state, cstate, τtr, τyini)
        failed(status) && return state.σ, status
    end

    # calculate Δσ
    Δτ = τ - τini
    Δσ = kn*Δu
    Δσ[1] = Δτ

    # update u and σ
    state.u = cstate.u + Δu
    state.σ = cstate.σ + Δσ

    return Δσ, success()
end


function state_values(mat::PowerExpBondSlip, state::PowerExpBondSlipState)
    return OrderedDict(
        :sl => state.u[1] ,
        :τl => state.σ[1] ,
        :p  => state.p ,
    )
end
