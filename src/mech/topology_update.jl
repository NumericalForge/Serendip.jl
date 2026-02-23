
# function fast_calc_cohesive_stress(elem::Element{MechCohesive})
#     # Check if the element should be opened
#     # according to normal and shear stresses at integration points of bulk elements
#     # coupled to the cohesive element

#     ndim = elem.ctx.ndim

#     # compute element center
#     m  = div(length(elem.nodes), 2)
#     C  = get_coords(elem.nodes[1:m])
#     Xc = vec(sum(C, dims=1)/m)

#     # gather all integration points from coupled bulk elements
#     ips = [ ip for coupling in elem.couplings for ip in coupling.ips ]
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
#     calc_interface_rotation(J, T)

#     # compute normal and shear stresses
#     n1 = T[:,1]
#     n2 = T[:,2]
#     if ndim==3
#         n3 = T[:,3]
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
    update_model_cohesive_elems(model, dofs, bcs)

Main update loop for the Extrinsic Cohesive Zone Model.
1. Checks stress criteria on inactive cohesive elements.
2. Activates elements that exceed the threshold.
3. specific topology check: If an element activates but is blocked by a 
   previously activated element, it triggers a node split ("unzipping").
4. Rebuilds the global node list to preserve bandwidth (Cuthill-McKee order).
5. Updates the global DOF vector and mapping.

# Returns
- `dof_map`: Mapping from new index to old index (`new_vec = old_vec[dof_map]`).
- `new_idxs`: Indices of newly created unknown DOFs in the reordered DOF vector.
- `affected_idxs`: DOF indices attached to touched/split nodes.
- `nu`: Number of unknown DOFs after reordering.
- `mobilized_elems`: Cohesive elements mobilized in this call.
"""
function update_model_cohesive_elems(model::FEModel, dofs::Vector{Dof}, bcs::Vector{BoundaryCondition})

    # Identify inactive cohesive candidates that are active in the stage.
    inactive_cohesives = [elem for elem in model.elems
                          if elem isa Element{MechCohesive} && elem.active && !elem.cache.mobilized]

    # Parent-node => children created by splitting.
    split_children_by_parent = Dict{Node, Vector{Node}}()
    split_nodes = Node[]

    topology_changed = false
    mobilized_elems = Element{MechCohesive}[]
    η_threshold = 0.7

    # Activate cohesive elements and split nodes as needed.
    for elem in inactive_cohesives
        η = stress_strength_ratio(elem)
        η >= η_threshold || continue
        # η >= 0.99 || continue

        elem.cache.mobilized = true
        push!(mobilized_elems, elem)
        
        m = div(length(elem.nodes), 2)
        face1_nodes = Set(elem.nodes[1:m])
        face2_nodes = Set(elem.nodes[m+1:end])

        # Explicit cohesive topology (already split) does not require runtime unzipping.
        isempty(intersect(face1_nodes, face2_nodes)) && continue

        unique_nodes = unique(elem.nodes)
        bulk1, bulk2 = elem.couplings[1], elem.couplings[2]

        for node in unique_nodes
            if !can_reach(bulk1, bulk2, node)
                # Split node and keep parent->children relation.
                new_node = split_node_cluster(model, node, bulk1)

                if !haskey(split_children_by_parent, node)
                    split_children_by_parent[node] = Node[]
                end
                push!(split_children_by_parent[node], new_node)
                push!(split_nodes, new_node)

                topology_changed = true
            end
        end
    end

    # Early exit if no changes
    if !topology_changed
        return Int[], Int[], Int[], 0, mobilized_elems
    end

    # Rebuild node list preserving locality (keeps parent/children close in numbering).
    reordered_nodes = Node[]
    
    # Helper to recursively add node and its children keeping matrix bandwidth low.
    # This ensures that if a new node (child) is split again (grandchild), 
    # the grandchild is also added.
    function add_node_and_children(n::Node)
        push!(reordered_nodes, n)
        if haskey(split_children_by_parent, n)
            for child in split_children_by_parent[n]
                add_node_and_children(child)
            end
        end
    end

    for node in model.nodes
        add_node_and_children(node)
    end

    model.nodes = reordered_nodes

    # Renumber nodes.
    for (i, node) in enumerate(model.nodes)
        node.id = i 
    end

    new_dofs = Dof[]
    for node in split_nodes
        append!(new_dofs, node.dofs)
    end
    append!(dofs, new_dofs)

    # Update BC connectivity after node splitting.
    update_node_bcs_for_new_nodes(split_nodes, bcs)
    refresh_boundary_connectivity(model, split_nodes)

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
    
    # Collect indices of new unknown DOFs.
    new_idxs = Int[]
    for dof in new_dofs
        dof.prescribed && continue
        push!(new_idxs, dof.eq_id)
    end

    # Compute affected DOF indexes (parents + children touched by split).
    affected_nodes = collect(keys(split_children_by_parent))
    append!(affected_nodes, split_nodes)

    affected_idxs = [ dof.eq_id for node in affected_nodes for dof in node.dofs ]

    return dof_map, new_idxs, affected_idxs, nu, mobilized_elems
end


# ❱❱❱ TOPOLOGICAL HELPER FUNCTIONS ------------------------------------------------


function update_node_bcs_for_new_nodes(new_nodes, bcs)
    essential = (:ux, :uy, :uz)
    for bc in bcs
        bc.kind == :node || continue
        all(k in essential for (k, _) in bc.conds) || continue
        
        for node in new_nodes
            matched = !isempty(select([node], bc.selector; nearest=false))
            matched || continue
            push!(bc.target, node)
            for (key, _) in bc.conds
                dof = get_dof(node, key)
                dof === nothing && continue
                dof.prescribed = true
            end
        end
    end
end


function refresh_boundary_connectivity(model::FEModel, new_nodes::Vector{Node})
    affected_owners = Set{Element}()
    for node in new_nodes
        for elem in node.elems
            push!(affected_owners, elem)
        end
    end

    for face in model.faces
        owner = face.owner
        owner in affected_owners || continue
        face.nodes = owner.nodes[owner.shape.facet_idxs[face.facet_idx]]
    end

    for edge in model.edges
        owner = edge.owner
        owner in affected_owners || continue
        edge.nodes = owner.nodes[owner.shape.edge_idxs[edge.facet_idx]]
    end
end



"""
    split_node_cluster(model, node, seed_elem)

Splits `node` into two based on connectivity to `seed_elem`.
Updates element connectivity pointers but DOES NOT add the new node to `model.nodes`.
"""
function split_node_cluster(model::FEModel, node::Node, seed_elem::AbstractCell)

    new_node = copy(node)
    new_node.id = -1
    new_node.elems = []

    # Deep copy DOFs so the new node has distinct degrees of freedom.
    if isdefined(node, :dofs)
        new_node.dofs = [deepcopy(d) for d in node.dofs]
    end

    # Identify connected bulk cluster that should move to the new node.
    group1 = get_connected_group(seed_elem, node)

    # Handle cohesive elements ("zipper" logic).
    cohesive_elems = [elem for elem in node.elems if elem.role == :cohesive]
    for cohesive in cohesive_elems
        m = div(length(cohesive.nodes), 2)

        if cohesive.couplings[1] in group1
            for i in 1:m
                cohesive.nodes[i] === node && (cohesive.nodes[i] = new_node)
            end
        end
        if cohesive.couplings[2] in group1
            for i in (m+1):length(cohesive.nodes)
                cohesive.nodes[i] === node && (cohesive.nodes[i] = new_node)
            end
        end

        on_old = any(n -> n === node, cohesive.nodes)
        on_new = any(n -> n === new_node, cohesive.nodes)

        !on_old && filter!(el -> el !== cohesive, node.elems)
        on_new && !(cohesive in new_node.elems) && push!(new_node.elems, cohesive)
    end

    # Handle special elements attached to moved host elements.
    special_elems = [elem for elem in node.elems if elem.role in (:line_interface, :tip)]
    for elem in special_elems
        host = elem.couplings[1]
        host in group1 && push!(group1, elem)
    end

    # Universal rewire (bulk + special).
    for elem in group1
        for i in 1:length(elem.nodes)
            elem.nodes[i] === node && (elem.nodes[i] = new_node)
        end
        elem in new_node.elems || push!(new_node.elems, elem)
        filter!(el -> el !== elem, node.elems)
    end

    return new_node
end

"""
    can_reach(start_elem, target_elem, pivot_node)
Performs a BFS to see if `start_elem` and `target_elem` are connected 
by INTACT (inactive) bulk material around `pivot_node`.
"""
function can_reach(start_elem::AbstractCell, target_elem::AbstractCell, pivot_node::Node)
    start_elem === target_elem && return true

    queue = AbstractCell[start_elem]
    head = 1
    visited = Set{Int}()
    push!(visited, start_elem.id)

    while head <= length(queue)
        current = queue[head]
        head += 1

        for neighbor in each_adjacent(current, pivot_node)
            neighbor.id in visited && continue

            # Active cohesive elements act as walls.
            is_interface_active(current, neighbor, pivot_node) && continue

            neighbor === target_elem && return true

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
function get_connected_group(start_elem::AbstractCell, pivot_node::Node)
    group = Set{AbstractCell}()
    push!(group, start_elem)
    queue = AbstractCell[start_elem]
    head = 1
    visited = Set{Int}()
    push!(visited, start_elem.id)

    while head <= length(queue)
        current = queue[head]
        head += 1

        for neighbor in each_adjacent(current, pivot_node)
            neighbor.id in visited && continue
            is_interface_active(current, neighbor, pivot_node) && continue

            push!(group, neighbor)
            push!(visited, neighbor.id)
            push!(queue, neighbor)
        end
    end

    return group
end


"""
    is_interface_active(e1, e2, pivot_node)
Checks if the interface between `e1` and `e2` is blocked by an ACTIVE cohesive element.
"""
function is_interface_active(e1::AbstractCell, e2::AbstractCell, pivot_node::Node)
    for cohesive in pivot_node.elems
        cohesive isa Element{MechCohesive} || continue

        c1 = cohesive.couplings[1]
        c2 = cohesive.couplings[2]

        if (c1 === e1 && c2 === e2) || (c1 === e2 && c2 === e1)
            return cohesive.cache.mobilized
        end
    end
    return false
end
