# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export LinearBondTip

mutable struct LinearBondTipState<:IpState
    ctx::Context
    f ::Float64
    w ::Float64
    function LinearBondTipState(ctx::Context)
        this = new(ctx)
        this.f = 0.0
        this.w = 0.0
        return this
    end
end


mutable struct LinearBondTip<:Constitutive
    k::Float64

    function LinearBondTip(;k=NaN)
        @check k>=0
        this = new(k)
        return this
    end
end


compat_state_type(::Type{LinearBondTip}, ::Type{MechBondTip}, ctx::Context) = LinearBondTipState


# # Type of corresponding state structure
# compat_state_type(::Type{LinearBondTip}) = LinearBondTipState

# # Element types that work with this material
# compat_elem_types(::Type{LinearBondTip}) = (MechBondTip,)


function calcD(mat::LinearBondTip, state::LinearBondTipState)
    return mat.k
end


function update_state(mat::LinearBondTip, state::LinearBondTipState, Δw)
    Δf = mat.k*Δw
    state.f += Δf
    state.w += Δw
    return Δf, success()
end


function state_values(mat::LinearBondTip, state::LinearBondTipState)
    return OrderedDict(
      :ur   => state.w ,
      :tau  => state.f )
end
