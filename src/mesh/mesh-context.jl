# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

mutable struct MeshContext
    ndim::Int                  # Mesh dimension

    function MeshContext(ndim=3)
        return new(ndim)
    end
end