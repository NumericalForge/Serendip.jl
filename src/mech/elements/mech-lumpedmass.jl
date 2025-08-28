# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

mutable struct MechLumpedMass<:MechFormulation
    id    ::Int
    shape ::CellShape

    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat::Constitutive
    active::Bool
    couplings::Array{Element,1}
    ctx::Context

    function MechLumpedMass()
        return new()
    end
end

compat_role(::Type{MechLumpedMass}) = VERTEX_CELL


function elem_stiffness(elem::MechLumpedMass)
    ndim = elem.ctx.ndim
    mat  = elem.cmodel
    K = zeros(ndim, ndim)

    keys = [:ux, :uy, :uz][1:ndim]
    map  = Int[ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    return K, map, map
end


function elem_mass(elem::MechLumpedMass)
    ndim = elem.ctx.ndim
    mat  = elem.cmodel

    M = mat.m*Matrix{Float64}(I, ndim, ndim)

    keys = [:ux, :uy, :uz][1:ndim]
    map  = Int[ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    return M, map, map
end


function update_elem!(elem::MechLumpedMass, U::Array{Float64,1}, Δt::Float64)
    return Float64[], Int[], success()
end


function elem_vals(elem::MechLumpedMass)
    vals = OrderedDict()
    return vals
end
