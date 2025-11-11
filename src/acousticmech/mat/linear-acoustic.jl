# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export LinearAcousticFluid

mutable struct LinearAcousticFluidState<:IpState
    ctx::Context
    up ::Float64          # pressure
    # V  ::Vector{Float64} # fluid velocity?

    function LinearAcousticFluidState(ctx::Context)
        this = new(ctx)
        # this.V  = zeros(ctx.ndim)
        this.up = 0.0
        return this
    end
end


mutable struct LinearAcousticFluid<:Constitutive
    # μ::Float64 # viscocity

    function LinearAcousticFluid(; kwargs...)
        return new()
    end
end


# Type of corresponding state structure
compat_state_type(::Type{LinearAcousticFluid}, ::Type{AcousticFluid}) = LinearAcousticFluidState


function update_state(mat::LinearAcousticFluid, state::LinearAcousticFluidState, Δup::Float64, G::Vector{Float64}, Δt::Float64)
    state.up += Δup
    return nothing
end


function state_values(mat::LinearAcousticFluid, state::LinearAcousticFluidState)
    D = OrderedDict{Symbol, Float64}()
    # D[:up] = state.up

    return D
end
