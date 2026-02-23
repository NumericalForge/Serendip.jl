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


"""
    dof_map(elem)

Returns the displacement DOF map for an element in local ordering expected by
its element routines.
"""
function dof_map(elem::Element{<:MechFormulation})
    ndim = elem.ctx.ndim
    keys = (:ux, :uy, :uz)[1:ndim]
    return Int[get_dof(node, key).eq_id for node in elem.nodes for key in keys]
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
