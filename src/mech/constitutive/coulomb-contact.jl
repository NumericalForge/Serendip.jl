# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export CoulombContact

"""
    CoulombContact(; kn, ks, mu)

Pure frictional constitutive model for interface/contact elements.

This model uses a Coulomb type friction criterion with no cohesion.

Yield function:
`f(σ) = tau + mu*σn`, where `tau = sqrt(tau1^2 + tau2^2)`.

# Keyword arguments
- `kn::Real`
  Normal stiffness per unit area (> 0).
- `ks::Real`
  Shear stiffness per unit area (> 0).
- `mu::Real`
  Friction coefficient (> 0).

# Returns
A `CoulombContact` object.
"""
mutable struct CoulombContact<:Constitutive
    kn ::Float64
    ks ::Float64
    μ  ::Float64

    function CoulombContact(; 
        kn::Real = NaN,
        ks::Real = NaN,
        mu::Real = NaN,
        )

        @check kn > 0 "CoulombContact: Normal stiffness per area kn must be > 0. Got $(repr(kn))."
        @check ks > 0 "CoulombContact: Shear stiffness per area ks must be > 0. Got $(repr(ks))."
        @check mu > 0 "CoulombContact: Friction coefficient mu must be > 0. Got $(repr(mu))."
        return new(kn, ks, mu)
    end
end


mutable struct CoulombContactState<:ConstState
    ctx::Context
    σ  ::Vec3        # stress
    w  ::Vec3        # relative displacements
    up ::Float64     # effective plastic relative displacement
    Δλ ::Float64     # plastic multiplier
    function CoulombContactState(ctx::Context)
        this    = new(ctx)
        this.σ  = zeros(Vec3)
        this.w  = zeros(Vec3)
        this.up = 0.0
        this.Δλ = 0.0
        return this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{CoulombContact}, ::Type{MechContact}) = CoulombContactState


function yield_func(mat::CoulombContact, σ::Vec3)
    σn, τ1, τ2 = σ
    τ = sqrt(τ1^2 + τ2^2)
    return τ + σn*mat.μ
end


function yield_derivs(mat::CoulombContact, σ::Vec3)
    _, τ1, τ2 = σ
    τ = sqrt(τ1^2 + τ2^2 + eps())
    return Vec3(mat.μ, τ1/τ, τ2/τ)
end


# Plastic potential derivatives for pure frictional flow.
function potential_derivs(mat::CoulombContact, σ::Vec3)
    _, τ1, τ2 = σ
    return Vec3(0.0, τ1, τ2)
end


function calc_kn_ks(mat::CoulombContact, state::CoulombContactState)
    return mat.kn, mat.ks
end


function calcD(mat::CoulombContact, state::CoulombContactState)
    kn, ks = calc_kn_ks(mat, state)

    De = @SMatrix [ kn   0.0  0.0
                    0.0  ks   0.0
                    0.0  0.0  ks ]

    if state.Δλ == 0.0
        return De
    elseif state.w[1] >= 0.0
        return @SMatrix zeros(3, 3)
    else
        _, τ1, τ2 = state.σ
        τ = sqrt(τ1^2 + τ2^2 + eps())
        μ = mat.μ

        Dep = @SMatrix [ kn           0.0                    0.0
                        -μ*kn*τ1/τ    ks*(1 - τ1^2/τ^2)     -ks*τ1*τ2/τ^2
                        -μ*kn*τ2/τ   -ks*τ1*τ2/τ^2           ks*(1 - τ2^2/τ^2) ]

        return Dep
    end
end


function update_state(mat::CoulombContact, state::CoulombContactState, cstate::CoulombContactState, Δw::Vector{Float64})
    kn, ks = calc_kn_ks(mat, state)
    De = @SMatrix [ kn   0.0  0.0
                    0.0  ks   0.0
                    0.0  0.0  ks ]


    # Trial state
    wtr = cstate.w + Δw

    σtr = cstate.σ + De*Δw

    # No traction in opening/non-compressive regime.
    # if wtr[1] >= 0.0
    #     state.Δλ = 0.0
    #     state.up = cstate.up + norm(Δw)
    #     state.w  = wtr
    #     state.σ  = Vec3(0.0, 0.0, 0.0)
    #     Δσ       = state.σ - cstate.σ

    #     return Δσ, success()
    # end

    # Trial yield function (compression branch)
    ftr = yield_func(mat, σtr)

    if ftr <= 0.0
        # Pure elastic increment
        state.Δλ = 0.0
        state.σ  = σtr
        state.up = cstate.up
    elseif σtr[1] >= 0.0
        state.Δλ = 1.0
        state.up = cstate.up + norm(Δw)
        state.w  = wtr
        state.σ  = Vec3(0.0, 0.0, 0.0)
        Δσ       = state.σ - cstate.σ

        return Δσ, success()
    else
        # Plastic increment
        σntr, τ1tr, τ2tr = σtr
        σntr += eps()
        τtr = sqrt(τ1tr^2 + τ2tr^2)

        state.Δλ = -1/ks * (τtr/(σntr*mat.μ) + 1)
        @assert state.Δλ >= 0.0
        state.σ  = Vec3(σntr, τ1tr/(1.0 + ks*state.Δλ), τ2tr/(1.0 + ks*state.Δλ))
        τ        = τtr/(1.0 + ks*state.Δλ)
        state.up = cstate.up + state.Δλ*τ
    end

    state.w = wtr
    Δσ      = state.σ - cstate.σ

    return Δσ, success()
end


function state_values(mat::CoulombContact, state::CoulombContactState)
    σn, τ1, τ2 = state.σ
    wn, s1, s2 = state.w
    τ = sqrt(τ1^2 + τ2^2)
    s = sqrt(s1^2 + s2^2)

    return Dict(
        :w  => wn,
        :s  => s,
        :σn => σn,
        :τ  => τ,
        :up => state.up,
    )
end


function output_keys(::CoulombContact)
    return Symbol[:w, :s, :σn, :τ, :up]
end
