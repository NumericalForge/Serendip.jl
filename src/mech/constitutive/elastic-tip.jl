# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export LinearTip, LinearTipContact

abstract type ElasticTip <: Constitutive end

mutable struct ElasticTipState<:IpState
    ctx::Context
    f ::Float64
    w ::Float64
    function ElasticTipState(ctx::Context)
        this = new(ctx)
        this.f = 0.0
        this.w = 0.0
        return this
    end
end


"""
    LinearTip(; k)

Constitutive model for a linear tip spring. Used to model the
elastic response at the tips of bar or beam elements inside bulk material.
Transmits tension and compression.

# Arguments
- `k::Real`
  Stiffness (≥ 0). Reaction update: `Δf = k·Δw`.

# Notes
- `w` is the relative tip displacement along the element axis.

# Returns
- A `LinearTip` object.
"""
mutable struct LinearTip<:ElasticTip
    k::Float64

    function LinearTip(;
            k::Real=NaN,
        )
        @check !isnan(k) "LinearTip: Stiffness k must be provided."
        @check k>=0.0 "LinearTip: Stiffness k must be non-negative."
        
        return new(k)
    end
end


"""
    LinearTipContact(; k)

Constitutive model for a tip contact. Used to model the
penetration response at the tips of bar or beam elements inside bulk material.
Transmits only compression using a linear law.

# Arguments
- `k::Real`
  Contact stiffness (≥ 0). Controls the reaction force (in the direction of the bar/beam element) per unit relative displacement.

# Returns
- A `LinearTipContact` object.
"""
mutable struct LinearTipContact<:ElasticTip
    k::Float64

    function LinearTipContact(;
            k::Real=NaN,
        )
        @check !isnan(k) "LinearTipContact: Stiffness k must be provided."
        @check k>=0.0 "LinearTipContact: Stiffness k must be non-negative."
        
        return new(k)
    end
end


# Type of corresponding state structure
compat_state_type(::Type{<:ElasticTip}, ::Type{MechBondTip}, ctx::Context) = ElasticTipState


# LinearTip
function calcD(mat::LinearTip, ::ElasticTipState)
    return mat.k
end


function update_state(mat::LinearTip, state::ElasticTipState, Δw::Float64)
    Δf = mat.k*Δw
    state.f += Δf
    state.w += Δw
    return Δf, success()
end


# LinearTipContact
function calcD(mat::LinearTipContact, state::ElasticTipState)
    if state.w>0.0 # penetration
        return mat.k
    else
        return 0.0
    end
end


function update_state(mat::LinearTipContact, state::ElasticTipState, Δw::Float64)
    fini = state.f
    ftr  = fini + mat.k*Δw

    if ftr>0.0
        f = ftr
    else
        f = 0.0
    end
    
    Δf = f - fini
    state.f += Δf
    state.w += Δw
    return Δf, success()
end


function state_values(::ElasticTip, state::ElasticTipState)
    return OrderedDict(
        :s => state.w,
        :τ => state.f
    )
end
