# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

abstract type Constitutive end


"""
`ConstState`

Abstract type for objects to store the state at integration points.
"""
abstract type ConstState
    #ctx::Context
    #other data
end


function init_state(::Constitutive, ::ConstState; args...)
end


# Default implementation for state commit
function commit_state(cstate::ConstState, state::ConstState)
    names = fieldnames(typeof(state))
    for name in names
        val = getfield(state, name)
        if isbits(val)
            setfield!(cstate, name, val)
        elseif typeof(val)<:AbstractArray
            copyto!(getfield(cstate,name), val)
        else
            setfield!(cstate, name, val)
        end
    end 
end


# Default implementation for state reset
function reset_state(state::ConstState, cstate::ConstState)
    commit_state(state, cstate)
end
   

function Base.copy(src::ConstState)
    T = typeof(src)
    dst = ccall(:jl_new_struct_uninit, Any, (Any,), T)
    names = fieldnames(T)
    for name in names
        val = getfield(src, name)
        if hasmethod(copy, (typeof(val),))
            setfield!(dst, name, copy(val))
        else
            setfield!(dst, name, val)
        end
    end
    return dst
end