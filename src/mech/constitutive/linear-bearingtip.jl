# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export LinearBearingTip

mutable struct LinearBearingTipState<:IpState
    ctx::Context
    f ::Float64
    w ::Float64
    function LinearBearingTipState(ctx::Context)
        this = new(ctx)
        this.f = 0.0
        this.w = 0.0
        return this
    end
end

# LinearBearingTip_params = [
#     FunInfo( :LinearBearingTip, "A model for a bar tip contact."),
#     KwArgInfo( :k, "Elastic stiffness", 1.0, cond=:(k>=0) ),
#     KwArgInfo( :fixed, "Flag to control if the tip is fixed", false, type=Bool),
# ]
# @doc docstring(LinearBearingTip_params) LinearBearingTip

mutable struct LinearBearingTip<:Constitutive
    k::Float64
    fixed::Bool

    function LinearBearingTip(;
        k::Real=1.0,
        fixed::Bool=false
        )
        @check k>=0.0 "Stiffness k must be non-negative."
        this = new(k, fixed)
        return this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{LinearBearingTip}, ::Type{MechBondTip}, evn::Context) = LinearBearingTipState


function calcD(mat::LinearBearingTip, state::LinearBearingTipState)
    if state.w>0.0 || mat.fixed
        return mat.k
    else
        return 0.0
    end
end


function update_state(mat::LinearBearingTip, state::LinearBearingTipState, Δw)
    fini = state.f
    ftr  = fini + mat.k*Δw

    if ftr>0.0 || mat.fixed
        f = ftr
    else
        f = 0.0
    end
    Δf = f - fini
    state.f += Δf
    state.w += Δw
    return Δf, success()
end


function state_values(mat::LinearBearingTip, state::LinearBearingTipState)
    return OrderedDict(
      :ur   => state.w ,
      :tau  => state.f )
end
