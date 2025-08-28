# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

# Rectangular Reissner Mindlin Plate FEM

export ElasticPlateRM

mutable struct ElasticPlateRMState<:IpState
    ctx::Context
    function ElasticPlateRMState(ctx::Context)
        return new(ctx)
    end
end

mutable struct ElasticPlateRM<:Constitutive
    E::Float64
    ν::Float64
    ρ::Float64

    function ElasticPlateRM(prms::Dict{Symbol,Float64})
        return  ElasticPlateRM(;prms...)
    end

    function ElasticPlateRM(;E=NaN, nu=NaN, ρ=0.0)
        E>0.0 || error("Invalid value for E: $E")
        (0<=nu<0.5) || error("Invalid value for nu: $nu")

        this = new(E, nu, ρ)
        return this
    end
end



# Type of corresponding state structure
compat_state_type(::Type{ElasticPlateRM}) = ElasticPlateRMState


function state_values(mat::ElasticPlateRM, state::ElasticPlateRMState)
    return OrderedDict{Symbol, Float64}()
end
