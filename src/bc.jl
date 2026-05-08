# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

mutable struct BoundaryCondition{T}
    kind  ::Symbol  # :node, :face, :curve, (:body for body loads)
    selector::Any
    conds ::AbstractDict
    target::Vector{T}
end


function configure_bc_dofs(bc::BoundaryCondition)
    if bc.kind == :node
        nodes = bc.target
    else
        nodes = [ node for facet in bc.target for node in facet.nodes ]
    end

    # Set prescribed essential bcs
    for node in nodes
        for (key, cond) in bc.conds
            dof = get_dof(node, key)
            dof===nothing && continue
            dof.name==key && (dof.prescribed=true)
        end
    end
end
