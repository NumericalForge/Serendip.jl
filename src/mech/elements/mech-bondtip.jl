# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export MechBondTip

struct MechBondTip<:MechFormulation
    MechBondTip() = new()
end

compat_role(::Type{MechBondTip}) = :tip


function set_quadrature(elem::Element{MechBondTip}, n::Int=1; state::NamedTuple=NamedTuple())
    ip = Ip([0.0, 0.0, 0.0], 0.0)
    ip.id = 1
    ip.state = compat_state_type(typeof(elem.cmodel), MechBondTip)(elem.ctx; state...)
    ip.owner = elem
    ip.coord = elem.nodes[end].coord

    elem.ips = [ ip ]
end


function mountB(elem::Element{MechBondTip}, Ch, Cm)
    # Calculates the matrix that relates nodal displacements with relative displacements


    # B = T* [-MM'  I]      ndim x ndim*(m+n)

    # where
    # T is a unit vector pointing outwards the member
    # MM is a matrix containing tresspased element shape functions
    # evaluated at the tip node coords

    # MM' = [ M_1*I M_2*I ... M_m*I]
    # I is a ndim x ndim identity matrix


    ndim    = elem.ctx.ndim
    tip     = elem.nodes[end]
    host    = elem.couplings[1]
    member  = elem.couplings[2]
    nsnodes = length(host.nodes)

    if hash(tip)==hash(member.nodes[1])
        R = [-1.0, 0, 0]
        D = member.shape.deriv(R)
        J = Cm'*D
        T = -normalize(J)
    else
        R = [+1.0, 0, 0]
        D = member.shape.deriv(R)
        J = Cm'*D
        T = normalize(J)
    end

    # Mount MM matrix
    stack = Array{Float64,2}[]
    X = tip.coord
    R = inverse_map(host.shape, Ch, X)
    M = host.shape.func(R)
    for Mi in M
        push!(stack, Mi*eye(ndim))
    end

    MM = hvcat(nsnodes, stack...)
    B = T'*[ -MM  eye(ndim) ]
    return B
end


function elem_stiffness(elem::Element{MechBondTip})
    ndim  = elem.ctx.ndim
    host  = elem.couplings[1]
    member = elem.couplings[2]
    Ch    = get_coords(host)
    Cm    = get_coords(member)

    B = mountB(elem, Ch, Cm)
    k = calcD(elem.cmodel, elem.ips[1].state)
    coef = k
    K = coef*B'*B

    keys = (:ux, :uy, :uz)[1:ndim]
    map = [ get_dof(node,key).eq_id for node in elem.nodes for key in keys ]
    return K, map, map
end

function elem_internal_forces(elem::Element{MechBondTip}, ΔUg::Vector{Float64}=Float64[], Δt::Float64=0.0)
    ndim   = elem.ctx.ndim
    host   = elem.couplings[1]
    member = elem.couplings[2]
    Ch     = get_coords(host)
    Cm     = get_coords(member)

    B    = mountB(elem, Ch, Cm)
    keys = (:ux, :uy, :uz)[1:ndim]
    map = [ get_dof(node,key).eq_id for node in elem.nodes for key in keys ]

    update = !isempty(ΔUg)

    if update
        ΔU = ΔUg[map]
        Δw = dot(B, ΔU)
        Δf, _ = update_state(elem.cmodel, elem.ips[1].state, Δw)
    else
        Δf = elem.ips[1].state.f
    end

    ΔF = Δf*B'

    return ΔF, map, success()
end
