
function fast_calc_cohesive_stress(elem::Element{MechCohesive})
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
    calc_interface_rotation(J, T)

    # compute normal and shear stresses
    n1 = T[:,1]
    n2 = T[:,2]
    if ndim==3
        n3 = T[:,3]
        t1 = dott(σ, n1)
        σn = dot(t1, n1)
        τ1 = dot(t1, n2)
        τ2 = dot(t1, n3)
        return [ σn, τ1, τ2 ]
    else
        n1 = Vec3(n1[1], n1[2], 0.0)
        n2 = Vec3(n2[1], n2[2], 0.0)
        t1 = dott(σ, n1)
        σn = dot(t1, n1)
        τ  = dot(t1, n2)
        return [ σn, τ ]
    end

end


"""
    update_model_cohesive_elems(model, dofs)

Main update loop for the Extrinsic Cohesive Zone Model.
1. Checks stress criteria on inactive cohesive elements.
2. Activates elements that exceed the threshold.
3. specific topology check: If an element activates but is blocked by a 
   previously activated element, it triggers a node split ("unzipping").
4. Rebuilds the global node list to preserve bandwidth (Cuthill-McKee order).
5. Updates the global DOF vector and mapping.

# Returns
- `flag`: True if mesh was modified.
- `dof_map`: Map from old DOF indices to new DOF indices.
- `new_free_idxs`: List of indices of the newly created DOFs (to initialize them to 0).
- `nu`: Number of unprescribed DOFs.
- `new_nodes_d`: Dict mapping Parent Node -> Vector{Child Nodes}.
"""
function update_model_cohesive_elems(model::FEModel, dofs::Vector{Dof})
    
    # Identify candidates
    candidates = [e for e in model.elems 
                  if e isa Element{MechCohesive} && e.active && !e.cache.mobilized]

    # Dictionary (Node -> Children)
    # This prevents ID collision for new nodes (which all share ID=-1)
    new_nodes = Node[]
    new_nodes_d = Dict{Node, Vector{Node}}()
    
    has_changes = false

    # ❱❱❱ Activate cohesive elements and split nodes as needed
    for elem in candidates
        # σs = calc_cohesive_stress(elem)
        
        # for (i, ip) in enumerate(elem.ips); ip.state.σ = σs[i]; end
        
        # η  = maximum(stress_strength_ratio(elem.cmodel, σ) for σ in σs)
        η = stress_strength_ratio(elem)
        # η >= 0.99 || continue
        η >= 0.99 || continue

        elem.cache.mobilized = true
        
        unique_nodes = unique(elem.nodes)
        bulk1, bulk2 = elem.couplings[1], elem.couplings[2]

        for node in unique_nodes
            if !can_reach(bulk1, bulk2, node)
                # Split
                new_node = split_node_cluster(model, node, bulk1)
                
                # Register using the NODE OBJECT as key
                if !haskey(new_nodes_d, node)
                    new_nodes_d[node] = Node[]
                end
                push!(new_nodes_d[node], new_node)
                push!(new_nodes, new_node)

                has_changes = true
            end
        end
    end

    # Early exit if no changes
    if !has_changes; return Int[], Int[], Int[], 0; end

    # ❱❱❱ Rebuild the node list to preserving bandwidth
    
    nodes_tmp = Vector{Node}()
    
    # Helper to recursively add node and its children keeping matrix bandwidth low.
    # This ensures that if a new node (child) is split again (grandchild), 
    # the grandchild is also added.
    function add_node_and_children(n::Node)
        push!(nodes_tmp, n)
        if haskey(new_nodes_d, n)
            for child in new_nodes_d[n]
                add_node_and_children(child)
            end
        end
    end

    for node in model.nodes
        add_node_and_children(node)
    end

    model.nodes = nodes_tmp

    # ❱❱❱ Renumber nodes
    for (i, node) in enumerate(model.nodes)
        node.id = i 
    end

    new_dofs = Dof[]
    for node in new_nodes
        append!(new_dofs, node.dofs)
    end
    append!(dofs, new_dofs)

    presc = [ dof.prescribed for dof in dofs ]
    pdofs = dofs[presc]
    udofs = dofs[.!presc]
    
    resize!(dofs, 0) # use resize! to keep the same array object
    append!(dofs, udofs)
    append!(dofs, pdofs)
    nu = length(udofs)

    dof_map = zeros(Int, length(dofs))
    for (i, dof) in enumerate(dofs)
        dof_map[i] = dof.eq_id
        dof.eq_id = i
    end
    
    # Collect indices of New (not prescribed) DOFs
    new_free_idxs = Int[]
    
    for dof in new_dofs
        dof.prescribed && continue
        push!(new_free_idxs, dof.eq_id)
    end

    # Compute affected DOF indexes
    affected_nodes = collect(keys(new_nodes_d))  
    append!(affected_nodes, new_nodes)

    affected_idxs = [ dof.eq_id for node in affected_nodes for dof in node.dofs ]

    return dof_map, new_free_idxs, affected_idxs, nu
end


# ❱❱❱ TOPOLOGICAL HELPER FUNCTIONS ------------------------------------------------

"""
    split_node_cluster(model, node, seed_elem)

Splits `node` into two based on connectivity to `seed_elem`.
Updates element connectivity pointers but DOES NOT add the new node to `model.nodes`.
"""
function split_node_cluster(model, node, seed_elem)
    
    new_node = copy(node)
    new_node.id = -1 
    new_node.elems = [] 
    
    # FIX: Deep copy DOFs so the new node has distinct degrees of freedom.
    # We use deepcopy to ensure we replicate the Dof values but create new Dof objects.
    if isdefined(node, :dofs)
        new_node.dofs = [deepcopy(d) for d in node.dofs]
    end
    
    # Identify Bulk Cluster (Group 1)
    group_1 = get_connected_group(seed_elem, node)
    
    # Handle Cohesive Elements (Zipper Logic)
    cohesive_candidates = [e for e in node.elems if e.role == :cohesive] 

    for cohe in cohesive_candidates
        m = div(length(cohe.nodes), 2)
        # Face 1
        if cohe.couplings[1] in group_1
            for i in 1:m
                if cohe.nodes[i] === node; cohe.nodes[i] = new_node; end
            end
        end
        # Face 2
        if cohe.couplings[2] in group_1
            for i in (m+1):length(cohe.nodes)
                if cohe.nodes[i] === node; cohe.nodes[i] = new_node; end
            end
        end
        
        # Fix Adjacency
        is_on_old = any(n -> n === node, cohe.nodes)
        if !is_on_old; filter!(el -> el !== cohe, node.elems); end
        
        is_on_new = any(n -> n === new_node, cohe.nodes)
        if is_on_new && !(cohe in new_node.elems); push!(new_node.elems, cohe); end
    end

    # Handle Special Elements (Rigid Logic)
    special_candidates = [e for e in node.elems if e.role in (:line_interface, :tip)]
    for e in special_candidates
        host = e.couplings[1]
        if host in group_1
            push!(group_1, e)
        end
    end

    # Universal Rewire (Bulk + Special)
    for elem in group_1
        for i in 1:length(elem.nodes)
            if elem.nodes[i] === node; elem.nodes[i] = new_node; end
        end
        push!(new_node.elems, elem)
        filter!(el -> el !== elem, node.elems)
    end

    return new_node
end

"""
    can_reach(start_elem, target_elem, pivot_node)
Performs a BFS to see if `start_elem` and `target_elem` are connected 
by INTACT (inactive) bulk material around `pivot_node`.
"""
function can_reach(start_elem, target_elem, pivot_node)
    if start_elem === target_elem; return true; end
    
    queue = [ start_elem ]
    visited = Set{Int}()
    push!(visited, start_elem.id)
    
    while !isempty(queue)
        current = popfirst!(queue)
        
        for neighbor in each_adjacent(current, pivot_node)
            
            neighbor.id in visited && continue
            
            # Blockage check: Active cohesive elements act as walls
            is_interface_active(current, neighbor, pivot_node) && continue

            # Target Found?
            neighbor.id == target_elem.id && return true
            
            push!(visited, neighbor.id)
            push!(queue, neighbor)
        end
    end
    
    return false
end

"""
    get_connected_group(start_elem, pivot_node)
Returns the Set of all bulk elements reachable from `start_elem` around `pivot_node` 
without crossing active interfaces.
"""
function get_connected_group(start_elem, pivot_node)
    group = Set{AbstractCell}()
    push!(group, start_elem)
    queue = [ start_elem ]
    visited = Set{Int}()
    push!(visited, start_elem.id)
    
    while !isempty(queue)
        current = popfirst!(queue)
        for neighbor in each_adjacent(current, pivot_node)
            if !(neighbor.id in visited) && !is_interface_active(current, neighbor, pivot_node)
                push!(group, neighbor)
                push!(visited, neighbor.id)
                push!(queue, neighbor)
            end
        end
    end
    return group
end


"""
    is_interface_active(e1, e2, pivot_node)
Checks if the interface between `e1` and `e2` is blocked by an ACTIVE cohesive element.
"""
function is_interface_active(e1::AbstractCell, e2::AbstractCell, pivot_node::Node)
    for elem in pivot_node.elems
        elem isa Element{MechCohesive} || continue

        c1 = elem.couplings[1]
        c2 = elem.couplings[2]

        if (c1 === e1 && c2 === e2) || (c1 === e2 && c2 === e1)
            return elem.cache.mobilized 
        end
    end
    return false
end
