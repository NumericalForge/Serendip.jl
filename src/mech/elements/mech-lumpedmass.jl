# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

mutable struct MechLumpedMass<:MechFormulation
    id    ::Int
    shape ::CellShape

    nodes ::Vector{Node}
    ips   ::Vector{Ip}
    tag   ::String
    mat::Constitutive
    active::Bool
    couplings::Vector{Element}
    ctx::Context

    function MechLumpedMass()
        return new()
    end
end

compat_role(::Type{MechLumpedMass}) = VERTEX_CELL


function dof_map(elem::MechLumpedMass)
    ndim = elem.ctx.ndim
    keys = (:ux, :uy, :uz)[1:ndim]
    return Int[node.dofdict[key].eq_id for node in elem.nodes for key in keys]
end


function elem_stiffness(elem::MechLumpedMass)
    ndim = elem.ctx.ndim
    mat  = elem.cmodel
    K = zeros(ndim, ndim)

    map = dof_map(elem)
    return K, map, map
end


function elem_mass(elem::MechLumpedMass)
    ndim = elem.ctx.ndim
    mat  = elem.cmodel

    M = mat.m*Matrix{Float64}(I, ndim, ndim)

    map = dof_map(elem)

    return M, map, map
end


function update_elem!(elem::MechLumpedMass, U::Vector{Float64}, Δt::Float64)
    return Float64[], Int[], success()
end


function elem_vals(elem::MechLumpedMass)
    vals = OrderedDict()
    return vals
end
