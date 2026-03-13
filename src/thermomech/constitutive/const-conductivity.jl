# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export ConstConductivity


mutable struct ConstConductivityState<:ConstState
    ctx::Context
    ut::Float64
    Q::Vector{Float64}
    function ConstConductivityState(ctx::Context)
        this = new(ctx)
        this.ut = 0.0
        this.Q  = zeros(ctx.ndim)
        return this
    end
end


mutable struct ConstConductivity<:Constitutive
    k ::Float64 # thermal conductivity with/m/K
    cv::Float64

    function ConstConductivity(; k, cv=0.0)
        k >= 0 || error("ConstConductivity: `k` must be nonnegative.")
        cv >= 0 || error("ConstConductivity: `cv` must be nonnegative.")
        return new(k, cv)
    end
end


# Type of corresponding state structure
compat_state_type(::Type{ConstConductivity}, ::Type{ThermoSolid}) = ConstConductivityState
if isdefined(@__MODULE__, :ThermoShell)
    @eval compat_state_type(::Type{ConstConductivity}, ::Type{ThermoShell}) = ConstConductivityState
end


function calc_cv(mat::ConstConductivity, ut::Float64) # Specific heat
    return mat.cv
end

function calcK(mat::ConstConductivity, state::ConstConductivityState) # Thermal conductivity matrix
    if state.ctx.ndim==2
        return mat.k*eye(2)
    else
        return mat.k*eye(3)
    end
end


function update_state(mat::ConstConductivity, state::ConstConductivityState, Δut::Float64, G::Vector{Float64}, Δt::Float64)
    K = calcK(mat, state)
    q = -K*G
    state.Q  += q*Δt
    state.ut += Δut
    return q
end


function state_values(mat::ConstConductivity, state::ConstConductivityState)
    D = OrderedDict{Symbol, Float64}()
    #D[:qx] = state.q[1]
    #D[:qy] = state.q[2]
    #if state.ctx.ndim==3
        #D[:qz] = state.q[3]
    #end
    return D
end
