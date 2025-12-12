"""
    add_contact_elements(mesh, selectors...; tag="", quiet=false)

Inserts contact elements along interfaces between distinct regions.
Simplifies the removal of conflicting cohesive elements by performing a post-generation check.
"""
function add_contact_elements(
    mesh::Mesh,
    selectors...;
    tag::String="",
    quiet::Bool=false,
)
    quiet || printstyled("Addition of contact elements:\n", bold=true, color=:cyan)

    # --- 1. Identify Target Cells ---
    if isempty(selectors)
        target_cells = select(mesh.elems, :bulk)
    else
        target_cells = Cell[]
        for s in selectors
            append!(target_cells, select(mesh.elems, s, :bulk))
        end
        unique!(target_cells)
    end
    isempty(target_cells) && error("add_contact_elements: No target cells found.")
    
    # Identify "Locked" cells (Special elements, or bulks not involved in contact)
    locked_cells = setdiff(mesh.elems, target_cells)

    # --- 1.5 Clean Up Existing Contacts (Same as before) ---
    active_tags = Set{String}(c.tag for c in target_cells)
    filter!(locked_cells) do cell
        if cell.role == :contact && length(cell.couplings) == 2
            t1 = cell.couplings[1].tag
            t2 = cell.couplings[2].tag
            if t1 in active_tags && t2 in active_tags
                return false 
            end
        end
        return true
    end

    # --- 2. Group by Tag ---
    cells_by_tag = Dict{String, Vector{Cell}}()
    for c in target_cells
        push!(get!(cells_by_tag, c.tag, Cell[]), c)
    end
    
    all_tags = sort(collect(keys(cells_by_tag)))
    if length(all_tags) < 2
        quiet || println("  Warning: Less than 2 tags found. No interfaces to split.")
        mesh.elems = [locked_cells; target_cells]
        return mesh
    end

    # --- 3. Topological Analysis & Geometric Matching ---
    node_to_tags = Dict{Node, Set{String}}()
    
    _pos_key(n) = round.(Tuple(n.coord), digits=8) .+ 0.0
    
    # Map: SortedCoords -> (Tag, FaceCell)
    face_geo_map = Dict{Vector{Tuple{Float64,Float64,Float64}}, Tuple{String, Cell}}()
    interfaces = [] # Stores: (FaceA, FaceB)

    for t in all_tags
        skin_faces = get_outer_facets(cells_by_tag[t])
        
        for face in skin_faces
            for n in face.nodes
                push!(get!(node_to_tags, n, Set{String}()), t)
            end
            
            geo_key = sort([_pos_key(n) for n in face.nodes])
            
            if haskey(face_geo_map, geo_key)
                other_tag, other_face = face_geo_map[geo_key]
                if other_tag != t
                    push!(interfaces, (other_face, face))
                    delete!(face_geo_map, geo_key) 
                end
            else
                face_geo_map[geo_key] = (t, face)
            end
        end
    end

    if isempty(interfaces)
        quiet || println("  No coincident interfaces found.")
        mesh.elems = [locked_cells; target_cells]
        return mesh
    end

    # --- 4. Node Splitting (Topology Update) ---
    replacements = Dict{String, Dict{Node, Node}}()
    generated_nodes = Node[]

    for (node, tags) in node_to_tags
        length(tags) < 2 && continue 
        
        sorted_tags = sort(collect(tags))
        
        for i in 2:length(sorted_tags)
            slave_tag = sorted_tags[i]
            
            new_node = copy(node)
            new_node.id = -1
            push!(generated_nodes, new_node)
            
            if !haskey(replacements, slave_tag)
                replacements[slave_tag] = Dict{Node, Node}()
            end
            replacements[slave_tag][node] = new_node
        end
    end

    # --- 5. Rewire Elements ---
    for cell in target_cells
        if haskey(replacements, cell.tag)
            replace_nodes!(cell, replacements[cell.tag])
        end
    end
    
    for cell in locked_cells
        if cell.role in (:line_interface, :tip)
            host = cell.couplings[1]
            if haskey(replacements, host.tag)
                replace_nodes!(cell, replacements[host.tag])
            end
        end
    end

    # --- 6. Create Contact Elements ---
    contact_cells = Cell[]
    generated_tags = Set{String}()
    
    # Store geometric keys of new contacts for collision detection
    contact_geo_keys = Set{Vector{Tuple{Float64,Float64,Float64}}}()

    for (f1, f2) in interfaces
        t1, t2 = f1.owner.tag, f2.owner.tag
        
        map1 = get(replacements, t1, Dict{Node,Node}())
        map2 = get(replacements, t2, Dict{Node,Node}())
        
        nodes1 = [get(map1, n, n) for n in f1.nodes]
        nodes2 = [get(map2, n, n) for n in f2.nodes]
        
        ctag = isempty(tag) ? "$(t1)-$(t2)" : tag
        push!(generated_tags, ctag)
        
        con = [nodes1; nodes2]
        c_elem = Cell(f1.shape, :contact, con, tag=ctag)
        c_elem.couplings = [f1.owner, f2.owner]
        push!(contact_cells, c_elem)
        
        # Register geometry (Use one face, it's sufficient for matching)
        # Note: A cohesive element at this interface would have nodes matching 
        # the geometric footprint of Face 1 (or Face 2).
        push!(contact_geo_keys, sort([_pos_key(n) for n in f1.nodes]))
    end

    # --- 7. Simplify: Post-Process Removal of Cohesive Elements ---
    # Scan locked_cells for cohesive elements.
    # If a cohesive element's face matches a contact element's face geometrically, delete it.
    
    cohesives_to_remove = Set{Cell}()
    
    # Iterate only over cohesive elements in locked_cells
    for cell in locked_cells
        if cell.role == :cohesive
            # A cohesive element has two faces. We check if its geometry matches
            # any of the contact interfaces we just created.
            # Since cohesive elements are typically 0-thickness or symmetric, 
            # sorting the coordinates of ALL nodes (or just half) gives a unique key.
            #
            # However, our contact_geo_keys stores the footprint of ONE face.
            # A cohesive element bridging that face would have nodes [FaceNodes; FaceNodes].
            # So its footprint is just the FaceNodes footprint (but duplicated).
            # The safest check: Check if the first half of the nodes matches a contact key.
            
            m = div(length(cell.nodes), 2)
            face_coords = sort([_pos_key(n) for n in cell.nodes[1:m]])
            
            if face_coords in contact_geo_keys
                push!(cohesives_to_remove, cell)
            end
        end
    end
    
    if !isempty(cohesives_to_remove)
        filter!(c -> !(c in cohesives_to_remove), locked_cells)
    end

    # --- 8. Finalize ---
    current_max = maximum(n.id for n in mesh.nodes)
    for (i, n) in enumerate(generated_nodes)
        n.id = current_max + i
    end
    
    append!(mesh.nodes, generated_nodes)
    mesh.elems = [locked_cells; target_cells; contact_cells]

    synchronize(mesh, sort=true, cleandata=true)

    if !quiet
        @printf "  %5d new contact cells\n" length(contact_cells)
        if !isempty(cohesives_to_remove)
            @printf "  %5d cohesive cells removed (conflict)\n" length(cohesives_to_remove)
        end
        @printf "  %5d nodes duplicated\n" length(generated_nodes)
        if !isempty(generated_tags)
            println("  Generated contact tags: ", join(sort(collect(generated_tags)), ", "))
        end
    end

    return mesh
end

# --- Helper Functions ---

"""
    replace_nodes!(cell, mapping)

Updates the connectivity of `cell` by replacing nodes found in `mapping`.
"""
function replace_nodes!(cell::AbstractCell, mapping::Dict{Node, Node})
    for (i, node) in enumerate(cell.nodes)
        if haskey(mapping, node)
            cell.nodes[i] = mapping[node]
        end
    end
end



"""
    add_boundary_contact_elements(mesh, selector; tag="", nodes_tag="", quiet=false)

Add contact interface elements to the boundary of `mesh`, coupling the original boundary faces with duplicated support nodes.
Useful for modeling surface contact or elastic supports (e.g., Winkler foundation).

Idempotency:
This function generates the new contact elements first, then checks the existing mesh 
to remove any old contact elements that match the **Owner** (Topology) and 
**Location** (Geometry) of the new elements. This safely handles updates without 
duplicating elements.

# Arguments
- `mesh::Mesh`: The mesh where boundary contact elements will be added.
- `selector`: Region selector (e.g. "bottom_face").
- `tag::String=""`: Tag for generated contact elements.
- `nodes_tag::String=""`: Tag for the duplicated support nodes.
- `quiet::Bool=false`: Suppress output.
"""
function add_boundary_contact_elements(
    mesh        :: Mesh,
    selector    :: Union{Expr,Symbolic,Tuple,String,Nothing}=nothing;
    tag         :: String="",
    nodes_tag   :: String="",
    quiet       :: Bool=false,
)
    quiet || printstyled("Addition of boundary interface elements:\n", bold=true, color=:cyan)

    if selector === nothing 
        throw(SerendipException("add_boundary_contact_elements: selector argument is required."))
    end

    # 1. Get target boundary faces
    faces = select(mesh.faces, selector)
    
    if isempty(faces)
        throw(SerendipException("add_boundary_contact_elements: no target cells found for selector $(repr(selector))."))
    end

    # 2. Duplicate Nodes
    #    If nodes are already split (e.g. Node A and Node B at same location),
    #    both will be in 'nodes_to_dup' because 'faces' contains faces for both sides.
    nodes_to_dup = Set{Node}()
    for face in faces
        for node in face.nodes
            push!(nodes_to_dup, node)
        end
    end
    
    new_nodes_map = Dict{Node, Node}()
    generated_nodes = Node[]
    
    for node in nodes_to_dup
        new_n = copy(node)
        new_n.id = -1 
        
        if !isempty(nodes_tag)
            new_n.tag = nodes_tag
        end
        
        push!(generated_nodes, new_n)
        new_nodes_map[node] = new_n
    end

    # 3. Create Contact Elements
    contact_cells = Cell[]
    
    for face in faces
        mesh_nodes = face.nodes
        support_nodes = [new_nodes_map[n] for n in mesh_nodes]
        con = [mesh_nodes; support_nodes]
        
        # Create Contact Element
        contact_cell = Cell(face.shape, :contact, con, tag=tag)
        
        # Link to Bulk Owner
        if face.owner !== nothing
            contact_cell.couplings = [face.owner]
        end
        
        push!(contact_cells, contact_cell)
    end

    # 4. Clean Up Existing Contacts (Post-Generation Check)
    #    We remove any OLD elements in the mesh that conflict with what we just generated.
    #    Conflict = Same Owner AND Same Geometry (Master side).
    
    _pos_key(n) = round.(Tuple(n.coord), digits=8) .+ 0.0
    
    # Build signatures of the NEW elements
    # Set of (OwnerID, GeometricKey)
    new_signatures = Set{Tuple{Int, Vector{Tuple{Float64,Float64,Float64}}}}()
    
    for c in contact_cells
        if !isempty(c.couplings)
            owner_id = c.couplings[1].id
            
            # Geometry of the "Master" side (first half of nodes)
            m = div(length(c.nodes), 2)
            geo_key = sort([_pos_key(n) for n in c.nodes[1:m]])
            
            push!(new_signatures, (owner_id, geo_key))
        end
    end
    
    # Filter mesh.elems to remove conflicts
    removed_count = 0
    filter!(mesh.elems) do cell
        if cell.role == :contact
            # Check if this old cell matches any new cell's signature
            if !isempty(cell.couplings)
                owner_id = cell.couplings[1].id
                
                m = div(length(cell.nodes), 2)
                geo_key = sort([_pos_key(n) for n in cell.nodes[1:m]])
                
                if (owner_id, geo_key) in new_signatures
                    removed_count += 1
                    return false # Remove
                end
            end
        end
        return true # Keep
    end

    # 5. Finalize IDs and Lists
    current_max = maximum(n.id for n in mesh.nodes)
    for (i, n) in enumerate(generated_nodes)
        n.id = current_max + i
    end

    append!(mesh.nodes, generated_nodes)
    append!(mesh.elems, contact_cells)

    # 6. Cleanup
    synchronize(mesh, sort=true, cleandata=true)

    if !quiet
        @printf "  %5d new contact cells\n" length(contact_cells)
        if removed_count > 0
            @printf "  %5d existing contact cells removed (update)\n" removed_count
        end
        @printf "  %5d support nodes created\n" length(generated_nodes)
        @printf "  %5d total cells\n" length(mesh.elems)
    end

    return mesh
end


"""
    add_boundary_shell_elements(mesh, selector; tag="", contact_tag="", quiet=false)

Add shell (surface) elements to the boundary of `mesh`, optionally coupling them with contact interface elements.

Idempotency:
This function generates the new boundary elements first, then checks the existing mesh 
to remove any old shell/contact elements that match the **Owner** (Topology) and 
**Location** (Geometry) of the new faces. This safely handles updates without 
duplicating elements or accidentally deleting neighbors across a crack.

# Arguments
- `mesh::Mesh`: The mesh where shell elements will be added.
- `selector::Union{Expr,Symbolic,Tuple,String}`: Region selector defining which boundary faces are converted to shells.
- `tag::String=""`: Tag assigned to the created shell elements.
- `contact_tag::String=""`: If provided, also create contact elements linking original faces and new shell nodes, tagged with this value.
- `quiet::Bool=false`: Suppress console output if `true`.

# Returns
- `Mesh`: The updated mesh including the new shell and, if applicable, contact elements.
"""
function add_boundary_shell_elements(
    mesh       :: Mesh,
    selector   :: Union{Expr,Symbolic,Tuple,String};
    tag        :: String="",
    contact_tag:: String="",
    quiet      :: Bool=false,
)
    quiet || printstyled("Addition of boundary shell elements:\n", bold=true, color=:cyan)
    
    faces = select(mesh.faces, selector)
    if isempty(faces)
        quiet || println("  Warning: No faces found for selector $(repr(selector))")
        return mesh
    end

    # --- 1. Generation Logic ---
    gen_contact_elems = !isempty(contact_tag)
    new_cells = Cell[]
    generated_nodes = Node[]

    if gen_contact_elems
        # --- Case A: Shells with Contact Interface ---
        nodes_to_dup = Set{Node}()
        for face in faces
            for node in face.nodes
                push!(nodes_to_dup, node)
            end
        end
        
        # Create Duplicate Nodes
        new_nodes_map = Dict{Node, Node}()
        
        for node in nodes_to_dup
            new_n = copy(node)
            new_n.id = -1 
            if !isempty(tag)
                new_n.tag = tag 
            end
            
            push!(generated_nodes, new_n)
            new_nodes_map[node] = new_n
        end

        # Create Elements
        for face in faces
            # A. Shell Element (Uses NEW nodes)
            shell_nodes = [new_nodes_map[n] for n in face.nodes]
            shell_cell = Cell(face.shape, :surface, shell_nodes, tag=tag)
            
            # Coupling: Shell is coupled to the Bulk Owner (indirectly via contact, but good for tracking)
            # Actually, standard practice: Shell couples to nothing (it's floating), 
            # or couples to the Contact Element? 
            # Let's keep it simple: Shells usually don't need couplings if Contact exists.
            push!(new_cells, shell_cell)

            # B. Contact Element (Links Old -> New)
            contact_con = [face.nodes; shell_nodes]
            contact_cell = Cell(face.shape, :contact, contact_con, tag=contact_tag)
            
            couplings = AbstractCell[]
            if face.owner !== nothing
                push!(couplings, face.owner)
            end
            push!(couplings, shell_cell)
            
            contact_cell.couplings = couplings
            push!(new_cells, contact_cell)
        end
        
    else
        # --- Case B: Attached Shells (Shared Nodes) ---
        for face in faces
            shell_cell = Cell(face.shape, :surface, face.nodes, tag=tag)
            if face.owner !== nothing
                shell_cell.couplings = [face.owner]
            end
            push!(new_cells, shell_cell)
        end
    end

    # --- 2. Clean Up Existing Elements (Post-Generation Check) ---
    # We remove any OLD elements in the mesh that conflict with what we just generated.
    # A conflict is defined as: Same Owner AND Same Geometry.
    
    # A. Collect signatures of the target faces
    _pos_key(n) = (round(n.coord.x, digits=8)+0.0, round(n.coord.y, digits=8)+0.0, round(n.coord.z, digits=8)+0.0)
    
    # Set of (OwnerID, GeometricKey)
    target_signatures = Set{Tuple{Int, Vector{Tuple{Float64,Float64,Float64}}}}()
    
    for f in faces
        if f.owner !== nothing
            geo_key = sort([_pos_key(n) for n in f.nodes])
            push!(target_signatures, (f.owner.id, geo_key))
        end
    end
    
    # B. Filter mesh.elems
    removed_count = 0
    filter!(mesh.elems) do cell
        if cell.role in (:surface, :contact)
            
            # Check Owner (Topology)
            # We assume coupling[1] is the Bulk Owner.
            if !isempty(cell.couplings)
                owner_id = cell.couplings[1].id
                
                # Check Geometry
                # Use the first half of nodes (Master side) for checking against the face
                m = (cell.role == :contact) ? div(length(cell.nodes), 2) : length(cell.nodes)
                cell_geo = sort([_pos_key(n) for n in cell.nodes[1:m]])
                
                # If both match, it's a duplicate/conflict
                if (owner_id, cell_geo) in target_signatures
                    removed_count += 1
                    return false # Remove
                end
            end
        end
        return true # Keep
    end

    # --- 3. Finalize ---
    if !isempty(generated_nodes)
        current_max = maximum(n.id for n in mesh.nodes)
        for (i, n) in enumerate(generated_nodes)
            n.id = current_max + i
        end
        append!(mesh.nodes, generated_nodes)
    end

    append!(mesh.elems, new_cells)
    synchronize(mesh, sort=true, cleandata=true)

    if !quiet
        @printf "  %4d dimensions\n" mesh.ctx.ndim
        @printf "  %5d total nodes\n" length(mesh.nodes)
        @printf "  %5d new cells created\n" length(new_cells)
        if removed_count > 0
            @printf "  %5d old cells removed (replaced)\n" removed_count
        end
        @printf "  %5d total cells\n" length(mesh.elems)
    end

    return mesh
end