# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export LinearBondSlip, ElasticRSJoint, ElasticBondSlip

mutable struct LinearBondSlipState<:IpState
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

# LinearBondSlip_params = [
#     FunInfo(:LinearBondSlip, "Elastic material for a rod-solid interface."),
#     KwArgInfo(:ks, "Shear stiffness", cond=:(ks>=0)),
#     KwArgInfo(:kn, "Normal stiffness", cond=:(kn>0)),
# ]
# @doc docstring(LinearBondSlip_params) LinearBondSlip(; kwargs...)

mutable struct LinearBondSlip<:Constitutive
    ks::Float64
    kn::Float64

    function LinearBondSlip(; ks=NaN, kn=NaN)
        @check kn > 0.0
        @check ks >= 0.0
        return new(ks, kn)
    end
end

const ElasticBondSlip = LinearBondSlip
const ElasticRSJoint = LinearBondSlip


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


function update_state(mat::LinearBondSlip, state::LinearBondSlipState, Δu::AbstractVector{Float64})
    D = calcD(mat, state)
    Δσ = D*Δu

    state.u .+= Δu
    state.σ .+= Δσ
    return Δσ, success()
end


function state_values(mat::LinearBondSlip, state::LinearBondSlipState)
    return OrderedDict(
        :s => state.u[1],
        :τ => state.σ[1]
    )
end
