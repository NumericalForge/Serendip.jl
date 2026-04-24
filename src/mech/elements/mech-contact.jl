# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export MechContact

"""
    MechContact()

Zero-thickness contact interface formulation for mechanical analyses.
"""
struct MechContact<:MechFormulation
end


# Return the shape family that works with this element
compat_role(::Type{MechContact}) = :contact

mutable struct MechContactCache <: ElementCache
    depth_list::FixedSizeVector{Float64}
end


function elem_init(elem::Element{MechContact})
    nips = length(elem.ips)
    depth_list = FixedSizeVector{Float64}(undef, nips)

    if elem.ctx.ndim == 2 && elem.ctx.stress_state != :axisymmetric
        owners = [owner for owner in elem.couplings if owner isa Element{MechSolid}]

        for (i, ip) in enumerate(elem.ips)
            if isempty(owners)
                depth_list[i] = elem.ctx.thickness
            else
                depth_list[i] = minimum(
                    evaluate(owner.etype.thickness, x=ip.coord.x, y=ip.coord.y, z=ip.coord.z)
                    for owner in owners
                )
            end
        end
    else
        depth_list .= 1.0
    end

    elem.cache = MechContactCache(depth_list)
    return nothing
end


function elem_stiffness(elem::Element{MechContact})
    ndim   = elem.ctx.ndim
    nnodes = length(elem.nodes)
    hnodes = div(nnodes, 2) # half the number of total nodes
    nstr   = 3

    C = get_coords(elem)[1:hnodes,:]
    B = fzeros(nstr, nnodes*ndim)
    K = fzeros(nnodes*ndim, nnodes*ndim)

    DB = fzeros(nstr, nnodes*ndim)
    J  = fzeros(ndim, ndim-1)
    NN = fzeros(nstr, nnodes*ndim)

    for (i, ip) in enumerate(elem.ips)
        depth = elem.ctx.stress_state == :axisymmetric ? 2*pi*ip.coord.x : ndim == 2 ? elem.cache.depth_list[i] : 1.0

        # compute shape Jacobian
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)

        @mul J = C'*dNdR
        detJ = norm2(J)

        # compute B matrix
        for i in 1:hnodes
            for dof in 1:ndim
                NN[dof, (i-1)*ndim + dof              ] = -N[i]
                NN[dof, hnodes*ndim + (i-1)*ndim + dof] =  N[i]
            end
        end
        
        T = calc_interface_rotation(J)
        @mul B = T'*NN

        # compute K
        coef = detJ*ip.w*depth
        D    = calcD(elem.cmodel, ip.state)

        @mul DB = D*B
        @mul K += coef*B'*DB
    end

    keys = (:ux, :uy, :uz)[1:ndim]
    map  = Int[ get_dof(node, key).eq_id for node in elem.nodes for key in keys ]
    return K, map, map
end


function elem_internal_forces(elem::Element{MechContact}, ΔU::Vector{Float64}=Float64[], dt::Float64=0.0)
    ndim   = elem.ctx.ndim
    nnodes = length(elem.nodes)
    hnodes = div(nnodes, 2) # half the number of total nodes
    map    = dof_map(elem)
    nstr   = 3

    update = !isempty(ΔU)
    if update
        Δω = zeros(nstr)
    end

    ΔF = fzeros(nnodes*ndim)
    C  = get_coords(elem)[1:hnodes,:]
    B  = fzeros(nstr, nnodes*ndim)

    J  = fzeros(ndim, ndim-1)
    NN = fzeros(nstr, nnodes*ndim)

    for (i, ip) in enumerate(elem.ips)
        depth = elem.ctx.stress_state == :axisymmetric ? 2*pi*ip.coord.x : ndim == 2 ? elem.cache.depth_list[i] : 1.0

        # compute shape Jacobian
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = norm2(J)

        # compute B matrix
        for i in 1:hnodes
            for dof in 1:ndim
                NN[dof, (i-1)*ndim + dof              ] = -N[i]
                NN[dof, hnodes*ndim + (i-1)*ndim + dof] =  N[i]
            end
        end

        T = calc_interface_rotation(J)
        @mul B = T'*NN

        if update
            @mul Δω = B*ΔU
            
            Δσ, status = update_state(elem.cmodel, ip.state, ip.cstate, Δω)
            failed(status) && return ΔF, map, status
        else
            Δσ = ip.state.σ
        end

        # internal force
        coef = detJ*ip.w*depth
        @mul ΔF += coef*B'*Δσ
    end

    return ΔF, map, success()
end


function elem_recover_nodal_values(elem::Element{MechContact})
    nips = length(elem.ips)

    keys = output_keys(elem.cmodel)
    vals = zeros(nips, length(keys))
    for (i,ip) in enumerate(elem.ips)
        dict = state_values(elem.cmodel, ip.state)
        vals[i,:] = [ dict[key] for key in keys ]
    end

    node_vals = OrderedDict{Symbol, Vector{Float64}}()
    E = extrapolator(elem.shape, nips)
    for (i,key) in enumerate(keys)
        V = E*vals[:,i]
        node_vals[key] = [ V; V ]
    end

    return node_vals
end
