# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

abstract type MechFormulation<:ElementFormulation end
# abstract type Element{<:MechFormulation})<:Element{MechBulk} end


"""
`elem_config_dofs(elem)`

Sets up the dofs for all nodes in `elem` according with its type.
This function can be specialized by concrete types.
"""
function elem_config_dofs(elem::Element{<:MechFormulation})
    for node in elem.nodes
        add_dof(node, :ux, :fx)
        elem.ctx.ndim>=2 && add_dof(node, :uy, :fy)
        elem.ctx.ndim==3 && add_dof(node, :uz, :fz)
    end
end


function reset_displacements(model::AbstractDomain)
    for n in model.nodes
        for dof in n.dofs
            for key in (:ux, :uy, :uz, :rx, :ry, :rz)
                haskey(dof.vals, key) && (dof.vals[key]=0.0)
            end
        end
    end
end