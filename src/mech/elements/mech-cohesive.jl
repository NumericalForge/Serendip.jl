# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export MechCohesive

mutable struct MechCohesive<:MechFormulation
    _state::Symbol # :idle, :leading, :split 
    MechCohesive() = new(:idle)
end

# Return the shape family that works with this element
compat_role(::Type{MechCohesive}) = :interface


function elem_init(elem::Element{MechCohesive})
    # Computation of characteristic length 'h' for cohesive elements
    # and set it in the integration point state

    hasfield(typeof(elem.ips[1].state), :h) || return

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

    for ip in elem.ips
        # compute shape Jacobian
        dNdR = elem.shape.deriv(ip.R)
        J    = C'*dNdR
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


function check_open_condition(elem::Element{MechCohesive})
    # Check if the element should be opened
    # according to normal and shear stresses at integration points of bulk elements
    # coupled to the cohesive element
    @assert elem.etype._state != :split
    ndim = elem.ctx.ndim

    ips = [ ip for elem in elem.couplings for ip in elem.ips ]
    Xs  = [ ip.coord for ip in ips ]
    nips = length(ips)

    # compute element center
    m  = div(length(elem.nodes), 2)
    C  = get_coords(elem.nodes[1:m])
    Xc = sum(C, dims=1)/m

    # number of selected ips
    n = min(2^ndim, nips)
    
    # find the n closest ips to the element center
    dists = [ norm(Xc - Xs[i]) for i in 1:nips ]
    perm  = partialsortperm(dists, 1:n)
    ips   = ips[perm]

    # average stress at selected ips
    σ = sum( ip.state.σ for ip in ips )/n

    # compute Jacobian and rotation matrix at element center
    R    = elem.shape.base_shape==TRI3 ? [1/3, 1/3] : [0.0, 0.0]
    dNdR = elem.shape.deriv(R)
    J = C'*dNdR
    T = fzeros(ndim, ndim)
    set_joint_rotation(J, T)

    # compute normal and shear stresses
    n1 = T[:,1]
    n2 = T[:,2]
    if ndim==3
        n3  = T[:,3]
        σn1 = σ*n1
        σn  = dot(σn1, n1)
        τ1  = dot(σn1, n2)
        τ2  = dot(σn1, n3)
        τ   = √(τ1*τ1 + τ2*τ2)
    else
        σn = dot(σ, n1)
        τ  = dot(σ, n2)
    end

    # check against maximum stresses if the element is leading
    if elem.etype._state == :leading
        for ip in elem.ips
            σn = max(σn, ip.state.σ[1])
            τ  = max(τ,  norm(ip.state.σ[2:end]))
        end
    end

    # check if the element should be opened
    if σn>0.8*elem.cmodel.ft || τ>0.5*elem.cmodel.ft
        elem.etype._state = :split
    end
end


function elem_stiffness(elem::Element{MechCohesive})
    elem.etype._state == :idle || return zeros(0,0), Int[], Int[] # return empty arrays if the element is idle

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
        
        set_joint_rotation(J, T)
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


function elem_internal_forces(elem::Element{MechCohesive}, ΔUg::Vector{Float64}=Float64[], Δt::Float64=0.0)
    elem.etype._state == :idle || return zeros(0), Int[], success() # return empty arrays if the element is idle

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

        set_joint_rotation(J, T)
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


function elem_recover_nodal_values(elem::Element{MechCohesive})
    nips = length(elem.ips)

    keys = output_keys(elem.cmodel)
    vals = zeros(nips, length(keys))
    for (i,ip) in enumerate(elem.ips)
        dict = state_values(elem.cmodel, ip.state)
        vals[i,:] = [ dict[key] for key in keys ]
    end

    node_vals = OrderedDict{Symbol, Array{Float64,1}}()
    E = extrapolator(elem.shape, nips)
    for (i,key) in enumerate(keys)
        V = E*vals[:,i]
        node_vals[key] = [ V; V ]
    end

    return node_vals
end
