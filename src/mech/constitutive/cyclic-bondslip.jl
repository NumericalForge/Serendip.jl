# This file ips part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export CyclicRSJoint, CyclicBondSlip

mutable struct CyclicBondSlipState<:IpState
    ctx::Context
    σ      ::Vector{Float64}
    u      ::Vector{Float64}
    τnl    ::Float64    # max stress
    τmax   ::Float64
    τres   ::Float64
    speak  ::Float64
    srev   ::Float64    # accumulated relative displacement
    sacum  ::Float64    # accumulated relative displacement
    sneg   ::Float64
    spos   ::Float64
    elastic::Bool
    function CyclicBondSlipState(ctx::Context)
        this         = new(ctx)
        ndim         = ctx.ndim
        this.σ       = zeros(ndim)
        this.u       = zeros(ndim)
        this.τnl     = 0.0
        this.τmax    = 0.0
        this.τres    = 0.0
        this.speak   = 0.0
        this.srev    = 0.0
        this.sacum   = 0.0
        this.sneg    = 0.0
        this.spos    = 0.0
        this.elastic = false
        return this
    end
end


"""
    CyclicBondSlip(; taumax, taures, speak, sres, alpha=0.4, beta=1.0, ks, kn, p)

Constitutive model for a cyclic rod–solid bond interface with degradation.
The shear bond law uses a piecewise envelope: power-law hardening up to `speak`,
a short plateau, then softening toward `taures`. Under cyclic loading, the
envelope evolves with the reversal amplitude.

# Arguments
- `taumax::Real`
  Initial peak bond stress τ_max (> 0).
- `taures::Real`
  Residual bond stress τ_res (≥ 0, < τ_max).
- `speak::Real`
  Slip at initial peak stress (> 0).
- `sres::Real`
  Slip where the envelope reaches τ_res (> 0).
- `alpha::Real=0.4`
  Exponent for the ascending branch (0 ≤ α ≤ 1).
- `beta::Real=1.0`
  Exponent for the descending branch (0 ≤ β ≤ 1).
- `ks::Real`
  Shear stiffness (> 0) with constraint `ks ≥ taumax/speak`.
- `kn::Real`
  Normal stiffness (> 0).
- `p::Real`
  Interface perimeter (> 0).

# Behavior
- Cyclic effects: peak τ_max decays, `speak` shifts, and `taures` evolves with a
  measure of reversal amplitude and accumulated slip (handled via internal state).
- Tangent shear stiffness equals the envelope derivative in inelastic updates,
  or `ks` when elastic.

# Returns
- A `CyclicBondSlip` object.

"""
mutable struct CyclicBondSlip<:Constitutive
    τmax:: Float64
    τres:: Float64
    speak:: Float64
    s2  :: Float64
    sres:: Float64
    α   :: Float64
    β   :: Float64
    ks  :: Float64
    kn  :: Float64
    p   :: Float64

    function CyclicBondSlip(; 
            taumax::Real = NaN,
            taures::Real = NaN,
            speak::Real  = NaN,
            sres::Real   = NaN,
            alpha::Real  = 0.4,
            beta::Real   = 1.0,
            ks::Real     = NaN,
            kn::Real     = NaN,
            p::Real      = NaN,
    )
        @check taumax > 0.0 "CyclicBondSlip: Shear strength taumax must be positive."
        @check taures >= 0.0 "CyclicBondSlip: Residual shear stress taures must be non-negative."
        @check speak > 0.0 "CyclicBondSlip: Peak slip speak must be positive."
        @check sres > 0.0 "CyclicBondSlip: Residual slip sres must be non-negative."
        @check 0.0 <= alpha <= 1.0 "CyclicBondSlip: Ascending curvature parameter alpha must be in [0, 1]."
        @check 0.0 <= beta <= 1.0 "CyclicBondSlip: Descending curvature parameter beta must be in [0, 1]."
        @check ks > 0.0 "CyclicBondSlip: Shear stiffness ks must be non-negative."
        @check kn > 0.0 "CyclicBondSlip: Normal stiffness kn must be positive."
        @check p > 0.0 "CyclicBondSlip: Perimeter p must be positive."
        @check taumax > taures "CyclicBondSlip: Shear strength taumax must be greater than residual shear stress taures."
        @check ks >= taumax/speak "CyclicBondSlip: Shear stiffness ks must be greater than or equal to taumax/speak."
        
        return new(taumax, taures, speak, 1.1*speak, sres, alpha, beta, ks, kn, p)
    end
end


# Type of corresponding state structure
compat_state_type(::Type{CyclicBondSlip}, ::Type{MechBondSlip}, ctx::Context) = CyclicBondSlipState


function tau(mat::CyclicBondSlip, ips::CyclicBondSlipState, s::Float64)
    s = abs(s)
    s2 = ips.speak*1.1
    if s<ips.speak
        return ips.τmax*(s/ips.speak)^mat.α
    elseif s<s2
        return ips.τmax
    elseif s<mat.sres
        return ips.τmax - (ips.τmax-ips.τres)*((s-s2)/(mat.sres-s2))^mat.β
    else
        return ips.τres
    end
end


function tau_deriv(mat::CyclicBondSlip, ips::CyclicBondSlipState, s::Float64)
    s = abs(s)
    s2 = ips.speak*1.1

    if s==0.0
        s1_factor = 0.01
        s = s1_factor*ips.speak   # to avoid undefined derivative
    end

    if s<=ips.speak
        return mat.α*ips.τmax/ips.speak*(s/ips.speak)^(mat.α-1)
    elseif s<s2
        return mat.ks*1e-3
    elseif s<mat.sres
        return -mat.β*(ips.τmax-ips.τres)/(mat.sres-s2)*((s-s2)/(mat.sres-s2))^(mat.β-1)
    else
        return mat.ks*1e-3
    end
end


function calcD(mat::CyclicBondSlip, ips::CyclicBondSlipState)
    ks = mat.ks
    kn = mat.kn

    if !ips.elastic
        s = ips.u[1]
        τ = ips.σ[1]

        if ips.τmax==0.0 && ips.τres==0.0 && ips.speak==0.0
            ips.τmax  = mat.τmax
            ips.τres  = mat.τres
            ips.speak = mat.speak
        end

        if s*τ<0.0 || abs(τ)>ips.τmax
            dτydsy = 1.0
        else
            dτydsy = tau_deriv(mat, ips, s)
        end
        ks = dτydsy
    end

    if ips.ctx.ndim==2
        return [  ks  0.0
                 0.0   kn ]
    else
        return  [  ks  0.0  0.0
                  0.0   kn  0.0
                  0.0  0.0   kn ]
    end
end


function update_state(mat::CyclicBondSlip, ips::CyclicBondSlipState, Δu::Vect)
    ks = mat.ks
    kn = mat.kn
    s  = ips.u[1]   # relative displacement
    Δs = Δu[1]      # relative displacement increment
    τini = ips.σ[1] # initial shear stress
    τtr  = τini + ks*Δs # elastic trial
    str  = s + Δs
    ips.sacum += abs(Δs)
    s>ips.spos && (ips.spos=s)
    s<ips.sneg && (ips.sneg=s)

    # amplitude
    sa = ips.spos-ips.sneg
    @assert sa>=0.0
    sh = sa*ips.srev/mat.speak^2

    ips.τmax = mat.τmax*exp(-0.0152*sh)

    ips.speak = mat.speak*(1 + 0.33*sh^0.5)

    if ips.sacum<mat.speak
        ips.τres = mat.τres*(ips.sacum/mat.speak)^mat.α
    else
        ips.τres = mat.τres*exp(-0.0027*sh)
    end


    if str*τtr<0
        τnl = ips.τres
    else
        τnl = max(ips.τres, tau(mat, ips, str))
    end

    ftr = abs(τtr) - τnl

    if ftr<0.0
        τ = τtr
        ips.elastic = true
    else
        if str*τtr<0
            Δsrev = (abs(τtr)-τnl)/ks
            @assert Δsrev>0.0
            ips.srev += Δsrev
        end
        τ = sign(τtr)*τnl
        Δτ = τ - τini
        ips.elastic = false
    end

    # calculate Δσ
    Δτ = τ - τini
    Δσ = kn*Δu
    Δσ[1] = Δτ

    # update u and σ
    ips.u .+= Δu
    ips.σ .+= Δσ

    return Δσ, success()
end


function stress_update2(mat::CyclicBondSlip, ips::CyclicBondSlipState, Δu::Vect)
    ks = mat.ks
    kn = mat.kn
    s  = ips.u[1]   # relative displacement
    Δs = Δu[1]      # relative displacement increment
    τini = ips.σ[1] # initial shear stress
    τtr  = τini + ks*Δs # elastic trial
    str  = s + Δs
    ips.sacum += abs(Δs)
    s>ips.spos && (ips.spos=s)
    s<ips.sneg && (ips.sneg=s)

    # amplitude
    sa = ips.spos-ips.sneg
    @assert sa>=0.0
    sh = sa*ips.srev/mat.sres^2

    # τmax = mat.τmax*(1 - (min(ips.srev, mat.sres)/mat.sres)^0.8)
    # τmax = mat.τmax*exp(-1.2*(ips.srev/mat.sres))
    # τmax = mat.τmax*min(1, 1.2*exp(-1.8*(ips.spos-ips.sneg)/mat.sres))
    # ips.τmax = mat.τmax*exp(-1.25*sh)
    ips.τmax = mat.τmax*exp(-1.02*sh)

    # ips.speak = mat.speak*(1 + log(1 + 5*ips.srev/mat.sres))
    ips.speak = mat.speak*(1 + 2.8*sh^0.5)

    if ips.sacum<mat.speak
        ips.τres = mat.τres*(ips.sacum/mat.speak)^mat.α
    else
        ips.τres = mat.τres*exp(-0.17*sh)
    end


    if str*τtr<0
        τnl = ips.τres
    else
        τnl = max(ips.τres, tau(mat, ips, str))
    end

    ftr = abs(τtr) - τnl

    if ftr<0.0
        τ = τtr
        ips.elastic = true
    else
        if str*τtr<0
            Δsrev = (abs(τtr)-τnl)/ks
            @assert Δsrev>0.0
            ips.srev += Δsrev
        end
        τ = sign(τtr)*τnl
        Δτ = τ - τini
        ips.elastic = false
    end

    # calculate Δσ
    Δτ = τ - τini
    Δσ = kn*Δu
    Δσ[1] = Δτ

    # update u and σ
    ips.u .+= Δu
    ips.σ .+= Δσ

    return Δσ, success()
end


function state_values(mat::CyclicBondSlip, ips::CyclicBondSlipState)
    return OrderedDict(
      :s => ips.u[1] ,
      :τ => ips.σ[1] ,
      )
end

