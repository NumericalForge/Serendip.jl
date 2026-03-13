# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export LinearInterface, LinearContact



# Super type for LinearInterface and LinearContact
abstract type ElasticInterface <: Constitutive end


"""
    LinearInterface(; kn, ks)

Bilateral linear elastic interface model.

Applies a penalty-type relation between relative displacement and traction
in both normal and shear directions. The normal response is linear in both
tension and compression; the shear response is linear elastic according to `ks`.

# Keyword arguments
- `kn::Float64` 
  Normal stiffness (> 0).
- `ks::Float64` 
  Shear stiffness (≥ 0).

# Returns
A `LinearInterface` constitutive object.
"""
mutable struct LinearInterface<:ElasticInterface
    kn::Float64 # Normal stiffness
    ks::Float64 # Shear stiffness

    function LinearInterface(; ks::Float64=NaN, kn::Float64=NaN)
        @check kn > 0.0
        @check ks >= 0.0
        return new(ks, kn)
    end
end


mutable struct ElasticInterfaceState<:ConstState
    ctx::Context
    σ  ::Vec3
    w  ::Vec3
    function ElasticInterfaceState(ctx::Context)
        this = new(ctx)
        this.σ = zeros(Vec3)
        this.w = zeros(Vec3)
        return this
    end
end


"""
    LinearContact(; kn, ks)

Unilateral linear elastic interface (contact) model.

Applies a penalty-type relation between relative displacement and traction
in compression and shear. In the normal direction, the response is linear
elastic in compression but traction is zero in tension. The shear response
is linear elastic according to `ks`.

# Keyword arguments
- `kn::Float64`
  Normal stiffness (> 0).
- `ks::Float64`
  Shear stiffness (≥ 0).

# Returns
A `LinearContact` constitutive object.
"""
mutable struct LinearContact<:ElasticInterface
    kn::Float64 # Normal stiffness
    ks::Float64 # Shear stiffness

    function LinearContact(; ks::Float64=NaN, kn::Float64=NaN)
        @check kn > 0.0
        @check ks >= 0.0
        return new(ks, kn)
    end
end


# Type of corresponding state structure
compat_state_type(::Type{<:ElasticInterface}, ::Type{MechContact}) = ElasticInterfaceState


# Functions for both LinearInterface and LinearContact

# LinearInterface
function calcD(mat::LinearInterface, state::ElasticInterfaceState)
    kn, ks = mat.kn, mat.ks
    return @SMatrix [  kn  0.0  0.0
                       0.0   ks  0.0
                       0.0  0.0   ks ]
end

# LinearContact
function calcD(mat::LinearContact, state::ElasticInterfaceState)
    ndim = state.ctx.ndim
    state.w[1] > 0.0 && return @SMatrix zeros(ndim, ndim)
    kn, ks = mat.kn, mat.ks
    return @SMatrix [  kn  0.0  0.0
                       0.0   ks  0.0
                       0.0  0.0   ks ]
end


function update_state(mat::ElasticInterface, state::ElasticInterfaceState, cstate::ElasticInterfaceState, Δu::AbstractArray)
    D  = calcD(mat, cstate)
    Δσ = D*Δu

    state.w = cstate.w + Δu
    state.σ = cstate.σ + Δσ
    return Δσ, success()
end


function state_values(mat::ElasticInterface, state::ElasticInterfaceState)
    σn, τ1, τ2 = state.σ
    τ = √(τ1^2 + τ2^2)
    wn, s1, s2 = state.w
    s = √(s1^2 + s2^2)

    return Dict(
        :w => state.w[1],
        :s  => s,
        :σn => state.σ[1],
        :τ  => τ,
    )
end


function output_keys(::ElasticInterface)
    return Symbol[:w, :σn, :s, :τ]
end