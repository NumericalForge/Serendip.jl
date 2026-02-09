# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export CebBondSlip, CebBondSlip

mutable struct CebBondSlipState<:ConstState
    ctx::Context
    σ  ::Vector{Float64}
    u  ::Vector{Float64}
    τy ::Float64      # max stress
    sy ::Float64      # accumulated relative displacement
    elastic::Bool
    function CebBondSlipState(ctx::Context)
        this = new(ctx)
        ndim = ctx.ndim
        this.σ = zeros(ndim)
        this.u = zeros(ndim)
        this.τy = 0.0
        this.sy = 0.0
        this.elastic = false
        return this
    end
end


"""
    CebBondSlip(; taumax, taures=0.0, s1, s2, s3, alpha=0.4, ks=taumax/s1, kn)

Constitutive model for bond–slip behavior of a reinforcing bar embedded in a solid,
according to the CEB (Comité Euro-International du Béton) formulation.

# Parameters
- `taumax::Float64` : Peak shear strength (> 0).
- `taures::Float64` : Residual shear stress (≥ 0, must be < τmax).
- `s1::Float64` : Characteristic slip defining the end of the ascending branch (> 0).
- `s2::Float64` : Characteristic slip at the start of the softening branch (> 0).
- `s3::Float64` : Characteristic slip where residual shear stress is reached (> 0).
- `alpha::Float64` : Curvature parameter for the ascending branch (0 ≤ α ≤ 1, default 0.4).
- `ks::Float64` : Initial shear stiffness (default = τmax/s1, must satisfy ks ≥ τmax/s₁).
- `kn::Float64` : Normal stiffness of the interface (> 0).

# Notes
- The model defines the shear stress–slip relation in three branches:
  ascending (0–s1), softening (s2–s3), and residual plateau (≥ s3).
- Default `ks` ensures consistency with the initial slope of the bond law.
- The normal stiffness `kn` penalizes opening displacement at the interface.

# References
CEB-FIP Model Code recommendations for bond–slip laws in reinforced concrete.
"""
mutable struct CebBondSlip<:Constitutive
    τmax:: Float64
    τres:: Float64
    s1  :: Float64
    s2  :: Float64
    s3  :: Float64
    α   :: Float64
    ks  :: Float64
    kn  :: Float64

    function CebBondSlip(;
        taumax::Real=NaN,
        taures::Real=0.0,
        s1::Real=NaN,
        s2::Real=NaN,
        s3::Real=NaN,
        alpha::Real=0.4,
        ks::Real=NaN,
        kn::Real=NaN
    )
        @check taumax > 0.0 "CebBondSlip: taumax must be > 0.0. Got $(taumax)"
        @check taures >= 0.0 "CebBondSlip: taures must be >= 0.0. Got $(taures)"
        @check taumax > taures "CebBondSlip: taures must be < taumax. Got taumax=$(taumax), taures=$(taures)"
        @check s1 > 0.0 "CebBondSlip: s1 must be > 0.0. Got $(s1)"
        @check s2 > 0.0 "s2 must be > 0.0. Got $(s2)"
        @check s3 > 0.0 "s3 must be > 0.0. Got $(s3)"
        @check 0.0 <= alpha <= 1.0 "CebBondSlip: alpha must be in [0.0, 1.0]. Got $(alpha)"
        
        if isnan(ks)
            ks = taumax/s1
        end
        if isnan(kn)
            kn = min(1e3*ks, 1e9)
        end
        @check kn > 0.0  "CebBondSlip: kn must be > 0.0. Current value $(kn)"
        @check ks >= taumax/s1 "ks must satisfy ks >= taumax/s1. Current value $(ks), taumax/s1=$(taumax/s1)"

        return new(taumax, taures, s1, s2, s3, alpha, ks, kn)
    end
end


compat_state_type(::Type{CebBondSlip}, ::Type{MechBondSlip}) = CebBondSlipState

# Type of corresponding state structure
compat_state_type(::Type{CebBondSlip}) = CebBondSlipState

# Element types that work with this material
compat_elem_types(::Type{CebBondSlip}) = (MechBondSlip,)


function Tau(mat::CebBondSlip, sy::Float64)
    if sy<mat.s1
        return mat.τmax*(sy/mat.s1)^mat.α
    elseif sy<mat.s2
        return mat.τmax
    elseif sy<mat.s3
        return mat.τmax - (mat.τmax-mat.τres)*(sy-mat.s2)/(mat.s3-mat.s2)
    else
        return mat.τres
    end
end


function deriv(mat::CebBondSlip, state::CebBondSlipState, sy::Float64)
    if sy==0.0
        s1_factor = 0.01
        sy = s1_factor*mat.s1   # to avoid undefined derivative
    end

    if sy<=mat.s1
        return mat.α*mat.τmax/mat.s1*(sy/mat.s1)^(mat.α-1)
    elseif sy<mat.s2
        return mat.ks*1e-3
    elseif sy<mat.s3
        return -(mat.τmax-mat.τres)/(mat.s3-mat.s2)
    else
        return mat.ks*1e-3
    end
end


function calcD(mat::CebBondSlip, state::CebBondSlipState)
    ks = mat.ks

    if !state.elastic
        dτydsy = deriv(mat, state, state.sy)
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


function yield_func(::CebBondSlip, state::CebBondSlipState, τ::Float64, τy::Float64)
    return abs(τ) - τy
end

function update_state(mat::CebBondSlip, state::CebBondSlipState, cstate::CebBondSlipState, Δu::Vect)
    ks = mat.ks
    kn = mat.kn
    Δs = Δu[1]      # relative displacement
    τini = cstate.σ[1] # initial shear stress
    τtr  = τini + ks*Δs # elastic trial

    ftr  = yield_func(mat, state, τtr, cstate.τy)

    if ftr<0.0
        τ = τtr
        state.elastic = true
    else
        # dτydsy = deriv(mat, state, state.sy)
        # @s ks
        # @s dτydsy
        # Δsy     = (abs(τtr)-state.τy)/(ks+dτydsy)
        # Δsy     = (abs(τtr)-state.τy)/abs(ks+dτydsy)
        Δsy       = (abs(τtr) - cstate.τy)/ks
        state.sy  = cstate.sy + Δsy
        state.τy  = Tau(mat, state.sy)
        τ  = state.τy*sign(τtr)
        Δτ = τ - τini
        state.elastic = false
    end

    # calculate Δσ
    Δτ = τ - τini
    Δσ = kn*Δu
    Δσ[1] = Δτ

    # update u and σ
    state.u .= cstate.u .+ Δu
    state.σ .= cstate.σ .+ Δσ

    return Δσ, success()
end


function state_values(mat::CebBondSlip, state::CebBondSlipState)
    return OrderedDict(
        :s => state.u[1] ,
        :τ => state.σ[1] ,
    )
end

