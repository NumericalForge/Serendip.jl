# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export MechCohesive

struct MechCohesive<:MechFormulation
    MechCohesive() = new()
end


# Return the shape family that works with this element
compat_role(::Type{MechCohesive}) = :cohesive


mutable struct MechCohesiveCache <: ElementCache
    open_state::Symbol
end


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

    elem.cache = MechCohesiveCache(:intact)
end



function elem_stiffness(elem::Element{MechCohesive})
    elem.cache.open_state == :intact && return zeros(0,0), Int[], Int[] # return empty arrays if the element is intact

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


function elem_internal_forces(elem::Element{MechCohesive}, ΔUg::Vector{Float64}=Float64[], Δt::Float64=0.0)
    elem.cache.open_state == :intact && return zeros(0), Int[], success() # return empty arrays if the element is intact

    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    hnodes = div(nnodes, 2) # half the number of total nodes
    keys   = (:ux, :uy, :uz)[1:ndim]
    map    = Int[ get_dof(node, key).eq_id for node in elem.nodes for key in keys ]

    update = !isempty(ΔUg)
    if update
        ΔU = ΔUg[map]
        # @show ΔU
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

    # @show ΔF

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

    node_vals = OrderedDict{Symbol, Vector{Float64}}()
    E = extrapolator(elem.shape, nips)
    for (i,key) in enumerate(keys)
        V = E*vals[:,i]
        node_vals[key] = [ V; V ]
    end

    return node_vals
end


# function calc_cohesive_stress(elem::Element{MechCohesive})
#     # Check if the element should be opened
#     # according to normal and shear stresses at integration points of bulk elements
#     # coupled to the cohesive element

#     ndim = elem.ctx.ndim

#     # compute element center
#     m  = div(length(elem.nodes), 2)
#     C  = get_coords(elem.nodes[1:m])
#     Xc = vec(sum(C, dims=1)/m)

#     # gather all integration points from coupled bulk elements
#     ips = [ ip for elem in elem.couplings for ip in elem.ips ]
#     Xip  = [ ip.coord for ip in ips ]
#     nips = length(ips)

#     # number of selected ips
#     n = min(2^ndim, nips)
    
#     # find the n closest ips to the element center
#     dists = [ norm(Xc - Xip[i]) for i in 1:nips ]
#     perm  = partialsortperm(dists, 1:n)
#     ips   = ips[perm]

#     # average stress at selected ips
#     σ = sum( ip.state.σ for ip in ips )/n

#     # compute Jacobian and rotation matrix at element center
#     R    = elem.shape.base_shape==TRI3 ? [1/3, 1/3] : [0.0, 0.0]
#     dNdR = elem.shape.deriv(R)
#     J = C'*dNdR
#     T = fzeros(ndim, ndim)
#     set_interface_rotation(J, T)

#     # compute normal and shear stresses
#     n1 = T[:,1]
#     n2 = T[:,2]
#     if ndim==3
#         n3 = T[:,3]
#         # @show typeof(σ)
#         # @show σ
#         t1 = dott(σ, n1)
#         σn = dot(t1, n1)
#         τ1 = dot(t1, n2)
#         τ2 = dot(t1, n3)
#         return [ σn, τ1, τ2 ]
#     else
#         n1 = Vec3(n1[1], n1[2], 0.0)
#         n2 = Vec3(n2[1], n2[2], 0.0)
#         t1 = dott(σ, n1)
#         σn = dot(t1, n1)
#         τ  = dot(t1, n2)
#         return [ σn, τ ]
#     end

# end

"""
    calc_cohesive_stress2(elem::Element{MechCohesive})

Calculates the stress tensor at each integration point of a cohesive element using quadratic regression.

The stress is estimated by performing a quadratic least-squares regression on the stress values
from all integration points of the neighboring bulk elements. For each integration point of the
cohesive element, the resulting polynomial is evaluated at its coordinate. Finally, the stress
tensor is projected onto the element's local coordinate system to obtain normal and shear components.

# Arguments
- `elem::Element{MechCohesive}`: The cohesive element for which to calculate the stress.

# Returns
- A vector of vectors, where each inner vector contains the normal and shear stress components
  for one integration point in the element's local coordinate system.

"""
function calc_cohesive_stress(elem::Element{MechCohesive})
    ndim = elem.ctx.ndim
    nnodes = length(elem.nodes)
    hnodes = div(nnodes, 2) # half the number of total nodes
    C = get_coords(elem.nodes[1:hnodes], ndim)
    # 1. Get all integration points from neighboring bulk elements
    ips = [ip for bulk_elem in elem.couplings for ip in bulk_elem.ips]
    nips = length(ips)

    # Helper for polynomial terms
    function reg_terms(x::Float64, y::Float64, nterms::Int64)
        nterms==4 && return ( 1.0, x, y, x*y )
        nterms==3 && return ( 1.0, x, y )
        return (1.0,)
    end

    function reg_terms(x::Float64, y::Float64, z::Float64, nterms::Int64)
        nterms==4  && return ( 1.0, x, y, z )
        return (1.0,)
    end

    if ndim==3
        nterms = nips>=4 ? 4 : 1
    else
        nterms = nips>=4 ? 4 : nips>=3 ? 3 : 1
    end

    # 2. Build the regression matrix M
    M = FMatrix{Float64}(undef, nips, nterms)
    for (i, ip) in enumerate(ips)
        x, y, z = ip.coord
        M[i, :] .= (ndim == 3) ? reg_terms(x, y, z, nterms) : reg_terms(x, y, nterms)
    end

    # 2. Build the coefficients matrix N
    nstr  = length(ips[1].state.σ)
    ncips = length(elem.ips)
    N     = FMatrix{Float64}(undef, ncips, nterms)
    for (i, ip) in enumerate(elem.ips)
        x, y, z = ip.coord
        N[i, :] .= (ndim == 3) ? reg_terms(x, y, z, nterms) : reg_terms(x, y, nterms)
    end

    N_inv_M = N*pinv(M)

    V = FixedSizeMatrix{Float64}(undef, nstr, ncips) # to store stress components at selected ips

    # 3. Perform regression for each stress component
    for i in 1:nstr
        W = [ ip.state.σ[i] for ip in ips ]
        # @show W
        V[i, :] = N_inv_M * W
    end

    projected_stress = Vector{Float64}[]
    T = fzeros(ndim, ndim)
    for (i,ip) in enumerate(elem.ips)
        # compute Jacobian and rotation matrix at element center
        dNdR = elem.shape.deriv(ip.R)
        J = C'*dNdR
        set_interface_rotation(J, T)
        
        # compute normal and shear stresses
        n1 = T[:,1]
        n2 = T[:,2]
        σ  = V[:, i]
        if ndim==3
            n3 = T[:,3]
            t1 = dott(σ, n1)
            σn = dot(t1, n1)
            τ1 = dot(t1, n2)
            τ2 = dot(t1, n3)
            push!(projected_stress, [ σn, τ1, τ2 ])
        else
            n1 = Vec3(n1[1], n1[2], 0.0)
            n2 = Vec3(n2[1], n2[2], 0.0)
            t1 = dott(σ, n1)
            σn = dot(t1, n1)
            τ  = dot(t1, n2)
            push!(projected_stress, [ σn, τ ])
        end
    end

    return projected_stress
end


function update_model_cohesive_elems(model::FEModel, dofs::Vector{Dof})
    cohesive_elems = [ elem for elem in model.elems if elem isa Element{MechCohesive} && elem.active && elem.cache.open_state != :split ]

    # dictionary to be used to find the reference node for each cohesive element for later node reordering
    cohe_node_d = Dict{Int, Int}()
    for cohe in cohesive_elems
        ref_node_id = maximum( n.id for n in cohe.nodes )
        for node in cohe.nodes
            cohe_node_d[cohe.id] = ref_node_id
        end
    end

    # list of new nodes ids for cohesive elements
    new_node_ids = [ Int[] for _ in 1:length(model.nodes) ]
    
    ndim = model.ctx.ndim
    n_id = length(model.nodes)
    nnodes_ini = length(model.nodes)

    for elem in cohesive_elems
        # check open condition
        # σ = calc_cohesive_stress(elem)
        # η = strength_utilization(elem.cmodel, σ)
        # η > 0.95 || continue

        σs = calc_cohesive_stress(elem)
        η = maximum( strength_utilization(elem.cmodel, σ) for σ in σs )
        # @show [ ip.state.σ for e in elem.couplings for ip in e.ips ]
        # @show [ strength_utilization(elem.cmodel, σ) for σ in σs ]
        η > 0.95 || continue

        # Sync stress
        for (i,ip) in enumerate(elem.ips)
            ip.state.σ = σs[i]
            # ip.state.σ = σ
        end

        # @show calc_cohesive_stress(elem)
        # @show σs
        
        # error()

        # function from mesh/gen-interface-elems.jl
        new_nodes_elem = split_cohesive_element(model, elem)
        # set the current cohesive element as split
        
        length(new_nodes_elem) == 0 && continue

        # select as reference node
        if length(new_nodes_elem)>0
            ref_node = cohe_node_d[elem.id]
            append!(new_node_ids[ref_node], [ n.id for n in new_nodes_elem ])
        end
        
    end # elem

    length(model.nodes) == nnodes_ini  && return false, Int[], Int[], 0
    new_nodes = model.nodes[nnodes_ini+1:end]
    
    new_dofs = [ dof for node in new_nodes for dof in node.dofs ]
    append!(dofs, new_dofs)    
    
    # ❱❱❱ Rebuild node list

    # @show new_node_ids
    
    _nodes = Node[]
    for (i, node) in enumerate(model.nodes)
        i > nnodes_ini && break
        push!(_nodes, node)
        i_nodes = [ model.nodes[id] for id in new_node_ids[i] ]
        append!(_nodes, i_nodes )
    end

    if length(_nodes) != length(model.nodes)
        # @show length(_nodes)
        # @show length(model.nodes)
        error("Error updating model nodes after splitting cohesive elements")
    end

    model.nodes = _nodes
    for (i, node) in enumerate(_nodes)
        node.id = i
    end

    # ❱❱❱ Split dofs
    presc = [ dof.prescribed for dof in dofs ]
    pdofs = dofs[presc]
    udofs = dofs[.!presc]

    resize!(dofs, 0) # used to keep the same reference instead of doing dofs = [ udofs; pdofs ]
    append!(dofs, udofs)
    append!(dofs, pdofs)
    nu = length(udofs)

    # ❱❱❱ maps from old dofs to new dofs
    map1 = zeros(Int, length(dofs)) # map for displacements
    map2 = zeros(Int, length(dofs)) # map for forces
    
    for (i,dof) in enumerate(dofs)
        map1[i] = dof.eq_id
        dof.eq_id = i # set new eq_id
    end
    
    map2 = copy(map1)
    for dof in new_dofs
        dof.prescribed && continue
        map2[dof.eq_id] = 0 # zero to avoid adding extra forces for new dofs
    end

    return true, map1, map2, nu
end