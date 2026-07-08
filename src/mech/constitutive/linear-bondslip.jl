# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export LinearBondSlip

mutable struct LinearBondSlipState<:ConstState
    ctx::Context
    σ ::Vector{Float64}
    u ::Vector{Float64}
    function LinearBondSlipState(ctx::Context)
        this = new(ctx)
        this.σ = zeros(ctx.ndim)
        this.u = zeros(ctx.ndim)
        return this
    end
end

"""
    LinearBondSlip(; ks, kn)

Linear elastic bond-slip interface model.

Applies a linear relation between relative displacement and traction for
rod-solid bond-slip interfaces. The first local direction uses the tangential
stiffness `ks`, while the remaining directions use the normal stiffness `kn`.

# Keyword arguments
- `ks::Float64`
  Tangential bond-slip stiffness (≥ 0).
- `kn::Float64`
  Normal stiffness (> 0).

# Returns
A `LinearBondSlip` constitutive object.

# Notes
- Compatible with `MechBondSlip` elements.
"""
mutable struct LinearBondSlip<:Constitutive
    ks::Float64
    kn::Float64

    function LinearBondSlip(; ks=NaN, kn=NaN)
        @check kn > 0.0
        @check ks >= 0.0
        return new(ks, kn)
    end
end


# Type of corresponding state structure
compat_state_type(::Type{LinearBondSlip}, ::Type{MechBondSlip}) = LinearBondSlipState


function calcD(mat::LinearBondSlip, state::LinearBondSlipState)
    ks = mat.ks
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


function update_state(mat::LinearBondSlip, state::LinearBondSlipState, cstate::LinearBondSlipState, Δu::AbstractVector{Float64})
    D  = calcD(mat, state)
    Δσ = D*Δu

    state.u .= cstate.u + Δu
    state.σ .= cstate.σ + Δσ
    return Δσ, success()
end


function state_values(mat::LinearBondSlip, state::LinearBondSlipState)
    return OrderedDict(
        :sl => state.u[1],
        :τl => state.σ[1]
    )
end
