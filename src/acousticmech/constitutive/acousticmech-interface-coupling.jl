# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export AcousticMechInterfaceCoupling


mutable struct AcousticMechInterfaceState<:ConstState
    ctx::Context

    function AcousticMechInterfaceState(ctx::Context)
        return new(ctx)
    end
end


struct AcousticMechInterfaceCoupling<:Constitutive
    function AcousticMechInterfaceCoupling(; kwargs...)
        return new()
    end
end


compat_state_type(::Type{AcousticMechInterfaceCoupling}, ::Type{<:AcousticMechInterface}) = AcousticMechInterfaceState


function state_values(mat::AcousticMechInterfaceCoupling, state::AcousticMechInterfaceState)
    return OrderedDict{Symbol, Float64}()
end
