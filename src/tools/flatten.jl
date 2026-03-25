# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

function _flatten(x, out=Any[])
    if x isa Tuple || x isa AbstractArray
        for item in x
            _flatten(item, out)
        end
    else
        push!(out, x)
    end
    return out
end


function _flatten_selectors(x, out=Any[])
    if x isa Tuple
        for item in x
            _flatten_selectors(item, out)
        end
    else
        push!(out, x)
    end
    return out
end


function _flatten(x, ::Type{T}, out=T[]) where {T}
    if x isa Tuple || x isa AbstractArray
        for item in x
            _flatten(item, T, out)
        end
    else
        x isa T || error("_flatten: expected components of type $(T), got $(typeof(x))")
        push!(out, x)
    end
    return out
end
