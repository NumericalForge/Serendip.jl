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
    LinearTip(; k, fixed=false)

Unilateral linear elastic tip spring.

Used with `MechBondTip` to model the elastic response at the tips of bar or
beam elements embedded in bulk material. By default, the spring is active only
for positive relative tip displacement (`w > 0`), which represents penetration.
Set `fixed=true` for a bilateral spring that remains active for both signs of
`w`.

# Keyword arguments
- `k::Real`
  Tip stiffness (≥ 0).
- `fixed::Bool=false`
  Keep the spring active for both positive and negative relative tip
  displacements.

# Notes
- `w` is the relative tip displacement along the element axis.

# Returns
A `LinearTip` constitutive object.
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


function update_state(mat::LinearTip, state::LinearTipState, cstate::LinearTipState, Δw::Float64)
    if mat.fixed || cstate.w + Δw > 0.0
        state.f = mat.k*Δw
    else
        state.f = 0.0
    end

    Δf = state.f - cstate.f
    state.w = cstate.w + Δw
    return Δf, success()
end


function state_values(::LinearTip, state::LinearTipState)
    return OrderedDict(
        :stip => state.w,
        :τtip => state.f
    )
end
