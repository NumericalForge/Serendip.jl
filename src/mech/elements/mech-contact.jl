# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export MechContact

struct MechContact<:MechFormulation
end

# Return the shape family that works with this element
compat_role(::Type{MechContact}) = :interface


function elem_init(elem::Element{MechContact})
    # Computation of characteristic length 'h' for cohesive elements
    # and set it in the integration point state

    hasfield(typeof(elem.ips[1].state), :h) || return

    ndim = elem.ctx.ndim

    # Avg volume of linked elements
    V = 0.0
    for elem in elem.couplings
        V += cell_extent(elem)
    end
    V /= length(elem.couplings)

    # Area of cohesive element
    A = 0.0
    C = get_coords(elem)
    n = div(length(elem.nodes), 2)
    C = C[1:n, :]
    J = fzeros(ndim, ndim-1)

    for ip in elem.ips
        # compute shape Jacobian
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = norm2(J)
        detJ <= 0 && error("Invalid Jacobian norm for cohesive element")
        A += detJ*ip.w
    end

    # Calculate and save h at cohesive element's integration points
    h = V/A
    for ip in elem.ips
        ip.state.h = h
    end

end


function elem_stiffness(elem::Element{MechContact})
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    hnodes = div(nnodes, 2) # half the number of total nodes

    C = get_coords(elem)[1:hnodes,:]
    B = fzeros(ndim, nnodes*ndim)
    K = fzeros(nnodes*ndim, nnodes*ndim)

    DB = fzeros(ndim, nnodes*ndim)
    J  = fzeros(ndim, ndim-1)
    T  = fzeros(ndim, ndim)
    NN = fzeros(ndim, nnodes*ndim)

    for ip in elem.ips
        if elem.ctx.stress_state==:axisymmetric
            th = 2*pi*ip.coord.x
        end

        # compute shape Jacobian
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)

        @mul J = C'*dNdR
        detJ = norm2(J)

        # compute B matrix
        for i in 1:hnodes
            for dof=1:ndim
                NN[dof, (i-1)*ndim + dof              ] = -N[i]
                NN[dof, hnodes*ndim + (i-1)*ndim + dof] =  N[i]
            end
        end
        
        set_interface_rotation(J, T)
        @mul B = T'*NN

        # compute K
        coef = detJ*ip.w*th
        D    = calcD(elem.cmodel, ip.state)
        @mul DB = D*B
        @mul K += coef*B'*DB
    end

    keys = (:ux, :uy, :uz)[1:ndim]
    map  = Int[ get_dof(node, key).eq_id for node in elem.nodes for key in keys ]
    return K, map, map
end


function elem_internal_forces(elem::Element{MechContact}, ΔUg::Vector{Float64}=Float64[], dt::Float64=0.0)
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    hnodes = div(nnodes, 2) # half the number of total nodes
    keys   = (:ux, :uy, :uz)[1:ndim]
    map    = Int[ get_dof(node, key).eq_id for node in elem.nodes for key in keys ]

    update = !isempty(ΔUg)
    if update
        ΔU = ΔUg[map]
        Δω = zeros(ndim)
    end

    ΔF = fzeros(nnodes*ndim)
    C  = get_coords(elem)[1:hnodes,:]
    B  = fzeros(ndim, nnodes*ndim)

    J  = fzeros(ndim, ndim-1)
    T  = fzeros(ndim, ndim)
    NN = fzeros(ndim, nnodes*ndim)


    for ip in elem.ips
        if elem.ctx.stress_state==:axisymmetric
            th = 2*pi*ip.coord.x
        end

        # compute shape Jacobian
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = norm2(J)

        # compute B matrix
        for i in 1:hnodes
            for dof=1:ndim
                NN[dof, (i-1)*ndim + dof              ] = -N[i]
                NN[dof, hnodes*ndim + (i-1)*ndim + dof] =  N[i]
            end
        end

        set_interface_rotation(J, T)
        @mul B = T'*NN

        if update
            @mul Δω = B*ΔU
            Δσ, status = update_state(elem.cmodel, ip.state, Δω)
            failed(status) && return ΔF, map, status
        else
            Δσ = ip.state.σ
        end

        # internal force
        coef = detJ*ip.w*th
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
