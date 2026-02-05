# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export LinearTip

mutable struct LinearTipState<:ConstState
    ctx::Context
    f ::Float64
    w ::Float64
    function LinearTipState(ctx::Context)
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
mutable struct LinearTip<:Constitutive
    k::Float64
    fixed::Bool

    function LinearTip(;
            k::Real=NaN,
            fixed::Bool=false,
        )
        isnan(k) && throw(SerendipException("LinearTip: Stiffness k must be provided."))
        k<0.0 && throw(SerendipException("LinearTip: Stiffness k must be non-negative."))
        return new(k, fixed)
    end
end


# # Type of corresponding state structure
# compat_state_type(::Type{<:LinearElasticTip}, ::Type{MechBondTip}) = LinearElasticTipState
compat_state_type(::Type{LinearTip}, ::Type{MechBondTip}) = LinearTipState


# LinearTip
function calcD(mat::LinearTip, state::LinearTipState)
    if mat.fixed || state.w>0.0 # fixed or penetration
        return mat.k
    else
        return 0.0
    end
end


function update_state(mat::LinearTip, state::LinearTipState, Δw::Float64)
    fini = state.f

    if mat.fixed || state.w + Δw > 0.0
        state.f = mat.k*Δw
    else
        state.f = 0.0
    end

    Δf = state.f - fini
    state.w += Δw
    return Δf, success()
end


function state_values(::LinearTip, state::LinearTipState)
    return OrderedDict(
        :s => state.w,
        :τ => state.f
    )
end
