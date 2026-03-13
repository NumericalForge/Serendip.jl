# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

abstract type ThermoMech <: ElementFormulation end


"""
`elem_config_dofs(elem)`

Sets up the dofs for all nodes in `elem` according to its thermo-mechanical type.
This function can be overloaded by concrete types.
"""
function elem_config_dofs(elem::Element{<:ThermoMech})
    for node in elem.nodes
        add_dof(node, :ut, :ft)
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        elem.ctx.ndim == 3 && add_dof(node, :uz, :fz)
    end
end


"""
`elem_init(elem)`

Sets up `elem` according to its thermo-mechanical type.
"""
function elem_init(elem::Element{<:ThermoMech})
    return nothing
end


"""
`update_elem!(elem, dU)`

Returns the force increment vector for `elem` due to the global dof increment `dU`.
This function must be redefined by concrete thermo-mechanical formulations.
"""
function update_elem!(elem::Element{<:ThermoMech}, dU::Vector{Float64}, Δt::Float64)
    error("update_elem! not defined for formulation $(typeof(elem.etype))")
end


"""
`elem_vals(elem)`

Returns element-level values intended to be constant along the element.
"""
function elem_vals(elem::Element{<:ThermoMech})
    return Dict{Symbol, Float64}()
end


"""
`elem_internal_forces(elem, F)`

Gets internal nodal forces from the current thermo-mechanical state.
This function must be defined by each concrete type.
"""
function elem_internal_forces(elem::Element{<:ThermoMech}, F::Vector{Float64})
end
