# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export MechCohesive

struct MechCohesive<:MechFormulation
    MechCohesive() = new()
end


# Return the shape family that works with this element
compat_role(::Type{MechCohesive}) = :interface


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
    # elem.cache.open_state == :intact && return zeros(0,0), Int[], Int[] # return empty arrays if the element is intact

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
    # elem.cache.open_state == :intact && return zeros(0), Int[], success() # return empty arrays if the element is intact

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


function calc_cohesive_σ(elem::Element{MechCohesive})
    # Check if the element should be opened
    # according to normal and shear stresses at integration points of bulk elements
    # coupled to the cohesive element

    ndim = elem.ctx.ndim

    # compute element center
    m  = div(length(elem.nodes), 2)
    C  = get_coords(elem.nodes[1:m])
    Xc = vec(sum(C, dims=1)/m)

    # gather all integration points from coupled bulk elements
    ips = [ ip for elem in elem.couplings for ip in elem.ips ]
    Xip  = [ ip.coord for ip in ips ]
    nips = length(ips)

    # number of selected ips
    n = min(2^ndim, nips)
    
    # find the n closest ips to the element center
    dists = [ norm(Xc - Xip[i]) for i in 1:nips ]
    perm  = partialsortperm(dists, 1:n)
    ips   = ips[perm]

    # average stress at selected ips
    σ = sum( ip.state.σ for ip in ips )/n

    # compute Jacobian and rotation matrix at element center
    R    = elem.shape.base_shape==TRI3 ? [1/3, 1/3] : [0.0, 0.0]
    dNdR = elem.shape.deriv(R)
    J = C'*dNdR
    T = fzeros(ndim, ndim)
    set_interface_rotation(J, T)

    # compute normal and shear stresses
    n1 = T[:,1]
    n2 = T[:,2]
    if ndim==3
        n3 = T[:,3]
        t1 = σ*n1
        σn = dot(t1, n1)
        τ1 = dot(t1, n2)
        τ2 = dot(t1, n3)
        return [ σn, τ1, τ2 ]
    else
        n1 = Vec3(n1[1], n1[2], 0.0)
        n2 = Vec3(n2[1], n2[2], 0.0)
        t1 = dot(σ, n1)
        σn = dot(t1, n1)
        τ  = dot(t1, n2)
        return [ σn, τ ]
    end

end


function update_cohesive_elems(model::FEModel, dofs::Vector{Dof})
    cohesive_elems = [ elem for elem in model.elems if elem isa Element{MechCohesive} && elem.active && elem.cache.open_state != :split ]
    new_nodes = Node[]
    
    nnodes = length(model.nodes)
    n_id   = nnodes

    for elem in cohesive_elems
        # check open condition
        σ = calc_cohesive_σ(elem)
        η = strength_utilization(elem.cmodel, σ)
        η > 0.98 || continue

        # check_open_condition(elem) == :split || continue
        m = div(length(elem.nodes),2) # half number of nodes
        new_nodes_subset = Node[]

        # iterate along half the nodes (since the other half are duplicates)
        for (k, node) in enumerate(elem.nodes[1:m])
            # only duplicate if the "pair" still references the same node
            node.id != elem.nodes[m+k].id && continue

            # bulk neighbors touching this node
            bulks = Element{MechBulk}[ e for e in node.elems if e.role == :bulk ]

            # get cohesive elements that share the node and are not split
            cohes = Element{MechCohesive}[ e for e in node.elems if e.role == :interface && e.cache.open_state != :split ]

            # get bulk elements linked to the cohesive elements
            buffer = Element{MechBulk}[]
            for co in cohes
                append!(buffer, co.couplings)
            end
            # keep only elements with two or more neighbors
            buffer = [ b for (b,cnt) in countmap(buffer) if cnt>=2 ]

            locked_bulks = Element{MechBulk}[]
            if length(buffer) != length(bulks)
                n = length(buffer)
                for i in 1:n
                    for j in i+1:n
                        # add to locked bulks if they share a facet
                        length(intersect(buffer[i].nodes, buffer[j].nodes)) < 2 && continue  # skip if they share less than two nodes
                        length(intersect(get_facets(buffer[i]), get_facets(buffer[j]))) > 0 && push!(locked_bulks, buffer[i])
                    end
                end

                # update the bulk list removing locked bulks
                bulks = setdiff(bulks, locked_bulks)
            end

            # duplicate node for each selected bulk
            skip_first = length(locked_bulks) == 0 # skip first bulk if no locked bulks
            for (i, bulk) in enumerate(bulks)
                
                i==1 && skip_first && continue # skip first bulk if no locked bulks
                
                # locate position of this node inside the bulk
                pos = findfirst( n->n.id==node.id, bulk.nodes )
                pos === nothing && error("Node $(node.id) not found in bulk element $(bulk.id)")
                
                new_node = copy(node)
                n_id += 1
                new_node.id = n_id
                # new_node_old_node_d[new_node.id] = node.id
                new_node.elems = [bulk] # start membership with this bulk only
                bulk.nodes[pos] = new_node
                push!(new_nodes_subset, new_node)

                # keep incidence consistent: remove bulk from old node
                filter!(e -> e.id !== bulk.id, node.elems)
            end

            # rewire intact cohesive neighbors that touch the duplicated bulk side
            for cohe in cohes
                cohe.cache.open_state = :leading # assume all as leading
                face1_nodes = cohe.nodes[1:m]
                face2_nodes = cohe.nodes[m+1:end]

                for (i, face_nodes) in enumerate((face1_nodes, face2_nodes))
                    bulk = cohe.couplings[i]  # bulk coupled to each face
                    pos_b = findfirst( n->hash(n)==hash(node), bulk.nodes )
                    pos_b === nothing && error("Node $(node.id) not found in bulk element $(bulk.id)")
                    bulk_node = bulk.nodes[pos_b]
                    bulk_node.id > nnodes || continue # only consider new nodes
                    pos_c = findfirst(n->hash(n)==hash(node), face_nodes)
                    push!(bulk_node.elems, cohe) # update nodal reference
                    face_nodes[pos_c] = bulk_node # update face nodes
                end
                cohe.nodes = vcat( face1_nodes, face2_nodes )

                # remove references to cohe from node
                if !(node.id in [ n.id for n in cohe.nodes ])
                    filter!(e -> e.id !== cohe.id, node.elems)
                end
            end

        end

        append!(new_nodes, new_nodes_subset)

        # set the current cohesive element as split
        elem.cache.open_state = :split
        
        # Sync stress
        for ip in elem.ips
            ip.state.σ = σ
        end
    
    end

    length(new_nodes) == 0  && return false, Int[], Int[], 0
    new_dofs = [ dof for node in new_nodes for dof in node.dofs ]

    # update global node list and dofs
    append!(model.nodes, new_nodes)
    append!(dofs, new_dofs)
    
    # Split dofs
    presc = [ dof.prescribed for dof in dofs ]
    pdofs = dofs[presc]
    udofs = dofs[.!presc]

    resize!(dofs, 0) # used to keep the same reference instead of doing dofs = [ udofs; pdofs ]
    append!(dofs, udofs)
    append!(dofs, pdofs)
    nu = length(udofs)

    # maps from old dofs to new dofs
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