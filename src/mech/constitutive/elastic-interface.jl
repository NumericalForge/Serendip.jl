# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export LinearInterface, LinearContact

mutable struct ElasticInterfaceState<:ConstState
    ctx::Context
    σ   ::Vector{Float64}
    w   ::Vector{Float64}
    function ElasticInterfaceState(ctx::Context)
        this = new(ctx)
        this.σ = zeros(ctx.ndim)
        this.w = zeros(ctx.ndim)
        return this
    end
end

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


function elastic_interface_D(ndim::Int, kn::Float64, ks::Float64)
    if ndim==2
        return [  kn  0.0
                 0.0   ks ]
    else
        return  [  kn  0.0  0.0
                  0.0   ks  0.0
                  0.0  0.0   ks ]
    end
end


# Functions for both LinearInterface and LinearContact

# LinearInterface
function calcD(mat::LinearInterface, state::ElasticInterfaceState)
    return elastic_interface_D(state.ctx.ndim, mat.kn, mat.ks)
end

# LinearContact
function calcD(mat::LinearContact, state::ElasticInterfaceState)
    ndim = state.ctx.ndim
    state.w[1] > 0.0 && return zeros(ndim, ndim)
    return elastic_interface_D(ndim, mat.kn, mat.ks)
end


function update_state(mat::ElasticInterface, state::ElasticInterfaceState, Δu)
    ndim = state.ctx.ndim
    D  = calcD(mat, state)
    Δσ = D*Δu

    state.w += Δu
    state.σ += Δσ
    return Δσ, success()
end


function state_values(mat::ElasticInterface, state::ElasticInterfaceState)
    ndim = state.ctx.ndim
    τ = norm(state.σ[2:ndim])
    s = norm(state.w[2:ndim])

    return Dict(
        :w => state.w[1],
        :s  => s,
        :σn => state.σ[1],
        :τ  => τ,
    )
end


function output_keys(::ElasticInterface)
    return Symbol[:w, :s, :σn, :τ]
end