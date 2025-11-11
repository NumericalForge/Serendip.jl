# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl


"""
    add_boundary_interface_elements(mesh, selector; tag="", nodes_tag="", quiet=false)

Adds interface elements to the boundary of `mesh`, useful for simulating boundary interactions such as elastic supports in the Winkler foundation model.

This function creates new interface elements at the faces selected by `selector`. Each face gets coupled with duplicated support nodes tagged by `nodes_tag`. The interface elements are assigned the specified `tag`.

# Arguments
- `mesh::Mesh`: The mesh where interface elements will be added.
- `selector::Union{Expr, Symbolic, String}`: Region selector for boundary faces where interface elements will be created.
- `tag::String=""`: Tag to assign to the generated interface elements.
- `nodes_tag::String`: Tag for the duplicated support nodes (required).
- `quiet::Bool=false`: Suppress console output if `true`.

# Returns
- `Mesh`: The updated mesh including the new boundary interface elements.

# Example
Use this to model elastic supports on the boundary:
```julia
add_boundary_interface_elements(mesh, selector="boundary_faces", tag="spring", nodes_tag="foundation")
```
"""
function add_boundary_contact_elements(
    mesh        :: Mesh,
    selector    :: Union{Expr,Symbolic,Tuple,String,Nothing}=nothing;
    tag         :: String="",
    nodes_tag   :: String="",
    quiet       :: Bool=false,
)

    quiet || printstyled("Addition of boundary interface elements:\n", bold=true, color=:cyan)

    selector === nothing && throw(SerendipException("add_boundary_interface_elements: selector argument is required."))

    # Get target cells
    faces = select(mesh.faces, selector)
    isempty(faces) && throw(SerendipException("add_boundary_interface_elements: no target cells found for selector $(repr(selector))."))

    # duplicate nodes
    nodes_to_dup = Set{Node}( node for face in faces for node in face.nodes )
    new_nodes_d  = Dict{UInt64, Node}( hash(p) => Node(p.coord, tag=nodes_tag) for p in nodes_to_dup )

    # Add contact elements
    contact_cells = Cell[]
    for face in faces
        con = copy(face.nodes)
        for (i, node) in enumerate(face.nodes)
            hs = hash(node)
            n  = new_nodes_d[hs]
            push!(con, n)
        end
        contact_cell = Cell(face.shape, :contact, con, tag=tag)
        contact_cell.couplings = [ face.owner ]
        push!(contact_cells, contact_cell)
    end

    # Update
    append!(mesh.elems, contact_cells)
    append!(mesh.nodes, collect(values(new_nodes_d)))

    # Update and reorder mesh
    synchronize(mesh, sort=true, cleandata=true)

    if !quiet
        @printf "  %5d new contact cells\n" length(contact_cells)
        @printf "  %5d total cells\n" length(mesh.elems)
    end

    return mesh
end


function add_boundary_shell_elements(
    mesh       :: Mesh,
    selector   :: Union{Expr,Symbolic,Tuple,String};
    tag        :: String="",
    contact_tag:: String="",
    quiet      :: Bool=false,
)

    quiet || printstyled("Addition of boundary shell elements:\n", bold=true, color=:cyan)
    
    faces = select(mesh.faces, selector)
    gen_contact_elems = contact_tag != ""

    # Add shell elements
    new_cells   = Cell[]
    new_nodes_d = Dict{UInt64, Node}()

    if gen_contact_elems
        nodes_to_dup = Set{Node}( node for face in faces for node in face.nodes )
        new_nodes_d  = Dict{UInt64, Node}( hash(p) => Node(p.coord, tag=p.tag) for p in nodes_to_dup )

        for face in faces
            con_sh = Node[]
            for node in face.nodes
                hs = hash(node)
                n  = new_nodes_d[hs]
                push!(con_sh, n)
            end
            shell_cell = Cell(face.shape, :surface, con_sh, tag=tag)
            push!(new_cells, shell_cell)

            contact_cell = Cell(face.shape, :contact, [face.nodes; con_sh], tag=contact_tag)
            contact_cell.couplings = [ face.owner, shell_cell ]
            push!(new_cells, contact_cell)
        end
        new_nodes = collect(values(new_nodes_d))
        append!(mesh.nodes, new_nodes)
    else
        for face in faces
            shell_cell = Cell(face.shape, :surface, face.nodes, tag=tag)
            push!(new_cells, shell_cell)
        end
    end

    # Update and reorder mesh
    append!(mesh.elems, new_cells)
    synchronize(mesh, sort=true, cleandata=true)

    if !quiet
        @printf "  %4d dimensions                           \n" mesh.ctx.ndim
        @printf "  %5d nodes\n" length(mesh.nodes)
        @printf "  %5d new cells\n" length(new_cells)
        # length(new_contact_cells)>0 && @printf("  %5d new interface cells\n", length(new_contact_cells))
        @printf "  %5d total cells\n" length(mesh.elems)
    end

    return mesh
end



function get_updated_nodes(old_nodes::Vector{Node}, new_nodes::Vector{Node})
    con = similar(old_nodes)
    for (i,p1) in enumerate(old_nodes)
        for p2 in new_nodes
            if hash(p1)==hash(p2)
                con[i] = p2
                break
            end
        end
    end
    return con
end


"""
    add_contact_elements(mesh, selectors...; tag="", quiet=false)

Insert contact elements along coincident interfaces in `mesh`.
Interfaces are identified between bulk regions with distinct tags
within the selected subset of elements.

# Arguments
- `mesh::Mesh`: The input mesh, which will be modified.
- `selectors::Union{Expr,Symbolic,Tuple,String}...`: Optional selectors defining which bulk regions to process. If not provided, all bulk elements are considered.
- `tag::String`: Tag assigned to all generated contact elements. If empty, tags are automatically derived from the connected regions.
- `quiet::Bool`: If `true`, suppresses console output.

Returns
- `Mesh`: The updated mesh including the generated contact elements.
"""
function add_contact_elements(
    mesh         ::Mesh,
    selectors    ::Union{Expr,Symbolic,Tuple,String}...;
    tag          ::String="",
    quiet        ::Bool=false,
)
    quiet || printstyled("Addition of contact elements:\n", bold=true, color=:cyan)

    # Target and locked cells: includes solids, lines, etc.
    if length(selectors)==0
        target_cells = select(mesh.elems, :bulk)
        locked_cells = setdiff(mesh.elems, target_cells)
    else
        target_cells = Cell[]
        for selector in selectors
            tc = select(mesh.elems, selector, :bulk)
            append!(target_cells, tc)
        end

        locked_cells = setdiff(mesh.elems, target_cells)
        # locked_cells = setdiff(locked_cells, select(locked_cells, selector, :cohesive)) # remove previous contact elements in filtered region
    end
    
    length(target_cells)==0 && throw(SerendipException("add_contact_elements: no target_cells found for selector $(repr(selector))."))

    # Get tags
    tag_set = Set{String}(cell.tag for cell in target_cells)
    length(tag_set)==0 && throw(SerendipException("add_contact_elements: no tags found. Tagged regions are required."))

    # Get contact faces
    trial_faces = CellFace[]

    # ❱❱❱ Iterate over tags
    for tag in tag_set
        tag_cells = select(target_cells, tag)
        tag_bulks = select(tag_cells, :bulk)
        tag_faces = get_outer_facets(tag_bulks)
        
        # duplicate nodes
        nodes_to_dup = Set{Node}( node for face in tag_faces for node in face.nodes )
        new_nodes_d  = Dict{UInt64, Node}( hash(p) => Node(p.coord, tag=p.tag) for p in nodes_to_dup )

        
        # update connectivities at outer target cells
        for cell in tag_bulks # do not use tag_outer_cells here
            for (i,p) in enumerate(cell.nodes)
                if p in nodes_to_dup
                    cell.nodes[i] = new_nodes_d[hash(p)]
                end
            end
        end
        
        # get target outer cells
        # tag_outer_cells = Set{Cell}( face.owner for face in tag_faces ) |> collect
        
        # update outer faces
        # tag_faces = get_outer_facets(tag_outer_cells)
        tag_faces = get_outer_facets(tag_bulks)
        append!(trial_faces, tag_faces)
    end

    # ❱❱❱ Get paired faces
    face_pairs  = Tuple{Cell, Cell}[]
    face_d      = Dict{UInt64, Cell}()

    for face in trial_faces
        hs = hash(face)
        f  = get(face_d, hs, nothing)
        if f===nothing
            face_d[hs] = face
        else
            push!(face_pairs, (face, f))
            delete!(face_d, hs)
        end
    end

    # ❱❱❱ Generate contact elements
    contact_cells = Cell[]
    tags = Set{String}()
    for (f1, f2) in face_pairs
        n   = length(f1.nodes)
        con = Array{Node}(undef, 2*n)
        k = 0
        for (i,p1) in enumerate(f1.nodes)
            for p2 in f2.nodes
                if hash(p1)==hash(p2)
                    k += 1
                    con[i]   = p1
                    con[n+i] = p2
                    break
                end
            end
        end
        k==n || error("add_contact_elements: faces f1 and f2 are not coincident.")

        if tag==""
            tagA, tagB = sort([f1.owner.tag, f2.owner.tag])
            tag = tagA*"-"*tagB
            push!(tags, tag)
        end
        cell = Cell(f1.shape, :contact, con, tag=tag)
        cell.couplings = [f1.owner, f2.owner]
        push!(contact_cells, cell)
    end

    # ❱❱❱ Remove overlapping cohesive elements
    contact_cells_d = Dict{UInt64, Cell}( hash(c) => c for c in contact_cells )
    cohesive_cells = select(locked_cells, :cohesive)
    cohesive_cells_to_remove = Cell[]
    for cell in cohesive_cells
        hs = hash(cell)
        if haskey(contact_cells_d, hs)
            push!(cohesive_cells_to_remove, cell)
        end
    end
    setdiff!(locked_cells, cohesive_cells_to_remove)

    # ❱❱❱ Fix cells connectivities for elements with couplings
    for c in locked_cells
        if c.role in (:tip, :line_interface)
            bulk = c.couplings[1]
            nspts = length(bulk.nodes)
            c.nodes[1:nspts] .= bulk.nodes
        elseif c.role == :cohesive
            n = length(c.nodes) ÷ 2
            bulk1, _ = c.couplings
            con = get_updated_nodes(c.nodes[1:n], bulk1.nodes)
            f.nodes = [ con; con ] # face nodes are the same in a mesh cohesive element
        end
    end
    
    # All cells
    mesh.elems = [ locked_cells; target_cells; contact_cells]

    # Include new nodes once
    nodes_d = Dict{Int,Node}()
    idx     = length(mesh.nodes)
    for cell in mesh.elems
        for node in cell.nodes
            if node.id < 0
                idx += 1
                node.id = idx # new id
            end
            nodes_d[node.id] = node
        end
    end

    # All nodes
    mesh.nodes = collect(values(nodes_d))

    # Update and reorder mesh
    synchronize(mesh, sort=true, cleandata=true)

    if !quiet
        @printf "  %5d new contact cells\n" length(contact_cells)
        @printf "  %5d total cells\n" length(mesh.elems)
        @printf "  %5d total nodes\n" length(mesh.nodes)
        if length(tags)>0
            s_tags = repr.(collect(tags))
            println("  generated contact tags: ", join(s_tags,", ", " and "))
        end
    end

    return mesh
end


"""
    add_cohesive_elements(mesh, selector=nothing; layers=2, tag="",
                          midnodes_tag="", inter_regions=false,
                          auto_tag=false, quiet=true)

Inserts cohesive (joint) elements into `mesh` to connect bulk elements.

- By default, joits are added globally. If `selector` is specified (element tag or expression), joints are added only in the selected region and all receive the specified `tag`.
- When `inter_regions=true`, joints are created only between regions defined by bulk element tags. In this case, tags are automatically generated for joints based on the regions they connect.
- `layers=2` creates standard joint elements, while `layers=3` inserts additional mid-layer nodes (assigned `midnodes_tag`, if provided).
- If `auto_tag=true`, tags for joints are automatically generated based on the connected regions.
- Set `quiet=false` to print summary information.

# Arguments
- `mesh::Mesh`: The input mesh object.
- `selector::Union{Expr, Symbolic, String, Nothing}`: Region selector (optional).
- `layers::Int`: Number of layers in joint elements (2 or 3). Defaults to 2.
- `tag::String`: Tag assigned to generated joints. Defaults to `""`.
- `midnodes_tag::String`: Tag assigned to mid-layer nodes (only if `layers=3`).
- `inter_regions::Bool`: If `true`, creates joints only between distinct regions.
- `auto_tag::Bool`: Automatically generate joint tags based on connected regions.
- `quiet::Bool`: Suppress console output if `true`.

# Returns
- `Mesh`: The updated mesh including the new cohesive elements.
"""
function add_cohesive_elements(
    mesh         ::Mesh,
    selector     ::Union{Expr,Symbolic,Tuple,String,Nothing}=nothing;
    tag          ::String="",
    quiet        ::Bool=false,
)

    quiet || printstyled("Addition of cohesive elements:\n", bold=true, color=:cyan)

    # Target and locked cells: includes bulk, lines, etc.
    if selector===nothing
        target_cells = select(mesh.elems, :bulk)
        locked_cells = setdiff(mesh.elems, target_cells)
        # remove previous cohesive elements at locked region
        locked_cells = setdiff(locked_cells, select(locked_cells, :cohesive))
    else
        target_cells = select(mesh.elems, selector, :bulk)
        length(target_cells)==0 && error("add_cohesive_elements: no target cells found for selector $selector")
        locked_cells = setdiff(mesh.elems, target_cells)
        # remove previous cohesive elements at filtered region
        locked_cells = setdiff(locked_cells, select(locked_cells, selector, :cohesive))
    end

    # Get trial faces
    locked_outer_faces = get_outer_facets(locked_cells)
    trial_faces = [ face for cell in target_cells for face in get_facets(cell) ]
    trial_faces = [ trial_faces; locked_outer_faces ]

    # ❱❱❱ Get face pairs
    face_pairs  = Tuple{Cell, Cell}[]
    face_d      = Dict{UInt64, Cell}()

    for face in trial_faces
        hs = hash(face)
        f  = get(face_d, hs, nothing)
        if f===nothing
            face_d[hs] = face
        else
            push!(face_pairs, (face, f))
            delete!(face_d, hs)
        end
    end

    # ❱❱❱ Generate cohesive elements
    cohesive_cells = Cell[]
    for (f1, f2) in face_pairs
        n   = length(f1.nodes)
        con = Array{Node}(undef, 2*n)
        k = 0
        for (i,p1) in enumerate(f1.nodes)
            for p2 in f2.nodes
                if hash(p1)==hash(p2)
                    k += 1
                    con[i]   = p1
                    con[n+i] = p2
                    break
                end
            end
        end
        k==n || error("add_cohesive_elements: faces f1 and f2 are not coincident.")

        cell           = Cell(f1.shape, :cohesive, con, tag=tag)
        cell.couplings = [f1.owner, f2.owner]
        push!(cohesive_cells, cell)
    end

    # ❱❱❱ Remove cohesive elements that overlap with contact elements
    contact_cells            = select(locked_cells, :contact)
    # @show length(contact_cells)
    # @show length(cohesive_cells)
    contact_cells_d = Dict{UInt64, Cell}( hash(c) => c for c in contact_cells )
    cells_to_remove = Cell[]
    for cell in cohesive_cells
        hs = hash(cell)
        if haskey(contact_cells_d, hs)
            push!(cells_to_remove, cell)
        end
    end
    # @show length(cells_to_remove)
    setdiff!(cohesive_cells, cells_to_remove)
    
    # @show length(cohesive_cells)

    # All cells
    mesh.elems = [ mesh.elems; cohesive_cells ]

    # Update and reorder mesh
    synchronize(mesh, sort=true, cleandata=true)

    if !quiet
        @printf "  %5d new cohesive elements\n" length(cohesive_cells)
        @printf "  %5d total cells\n" length(mesh.elems)
        @printf "  %5d total nodes\n" length(mesh.nodes)
    end

    return mesh
end


"""
    split_cohesive_element(model::AbstractDomain, elem::AbstractCell)

Splits a cohesive element to simulate crack opening by duplicating nodes at the crack front.

This function is intended to be called during a simulation when a cohesive element `elem` is determined to be near failure. It modifies the `model`'s topology in-place.

The process involves:
1.  Iterating through the node pairs that define the cohesive element.
2.  For each pair that has not yet been split (i.e., both sides reference the same node ID), it duplicates the node.
3.  The neighboring bulk elements are reconnected to these new, distinct nodes, effectively creating a separation.
4.  It intelligently "re-joins" nodes of adjacent bulk elements that are not separated by a cohesive element, ensuring the crack only propagates along the intended path.
5.  The connectivity of the cohesive element itself and other special elements (like `:tip` or `:line_interface`) is updated to reference the new nodes.
6.  The new nodes are added to the `model`.

# Arguments
- `model::AbstractDomain`: The finite element model, which will be modified in-place. It must contain all nodes and elements of the simulation.
- `elem::AbstractCell`: The cohesive element that is to be split. Its `open_state` cache field is updated to track its status.

# Notes
- This is a topology-modifying operation and should be used within an analysis step that handles such changes (e.g., in a fracture mechanics simulation).
- The function relies on the `role` of elements (`:bulk`, `:cohesive`, etc.) and the `couplings` field to understand the model's structure.
- New nodes are assigned temporary negative IDs before being renumbered and added to the model to avoid conflicts.

"""
function split_cohesive_element(model::AbstractDomain, elem::AbstractCell)

    # A cohesive element has two faces; m is the number of nodes per face.
    m = div(length(elem.nodes), 2)
    ndim = model.ctx.ndim
    new_nodes = Node[]


    # Iterate through the nodes of the first face of the cohesive element.
    # Each node `node` has a corresponding pair `elem.nodes[m+k]` on the opposite face.
    for (k, node) in enumerate(elem.nodes[1:m])

        n_neg_id = 0
        # If the node IDs are different, this node pair has already been split.
        node.id != elem.nodes[m+k].id && continue
        
        # Find all bulk elements connected to the current node.
        bulks = [ e for e in node.elems if e.role == :bulk ]

        # Find all cohesive elements (that are not yet fully split) connected to the node.
        cohes = [ e for e in node.elems if e.role == :cohesive && e.cache.open_state != :split ]

        # Find special elements like interface tips or line interfaces connected to the node.
        special_elems = [ e for e in node.elems if e.role in (:line_interface, :tip) ]

        # ❱❱ 1. Split: Duplicate the node for each neighboring bulk element.
        # We start by assuming the crack propagates everywhere around the node.
        # Later, we will merge nodes where the there is no cohesive element separating them.
        elem_pos_d = Dict{Int,Int}()
        for (i, bulk) in enumerate(bulks)
            # Store the local index of the node within the bulk element's connectivity.
            pos = findfirst( n->n.id==node.id, bulk.nodes )
            elem_pos_d[bulk.id] = pos
            
            # The first bulk element keeps the original node.
            if i == 1
                node.elems = [] # Clear its element affiliations; they will be rebuilt.
                continue
            end
            
            # For all other bulk elements, create a new node and replace the original one.
            new_node = copy(node)
            n_neg_id -= 1 # Assign a temporary negative ID to avoid conflicts.
            new_node.id = n_neg_id
            bulk.nodes[pos] = new_node
            new_node.elems = [] # No element affiliations yet.

            # push!(new_nodes, new_node)
        end

        # ❱❱ 2. Join: remove duplicated nodes from bulks pairs that are not linked by a cohesive element
        paired = [ e for co in cohes for e in co.couplings ]
        paired2p = [ b for (b,cnt) in countmap(paired) if cnt>=2 ]
        if length(paired2p) != length(bulks)

            # Identify pairs of bulk elements that are connected by a cohesive element.
            # These are the interfaces where the crack is allowed to exist.
            with_cohesive_pairs = Tuple{Int,Int}[]    
            for cohe in cohes
                b1 = cohe.couplings[1]
                b2 = cohe.couplings[2]
                min_id = min(b1.id, b2.id)
                max_id = max(b1.id, b2.id)
                push!(with_cohesive_pairs, (min_id, max_id))
            end

            # Identify pairs of bulk elements that are adjacent but *not* separated by a cohesive element.
            # These are the interfaces where the nodes should be re-joined.
            bare_pairs = Tuple{Int,Int}[]

            for i in 1:length(bulks)
                for j in i+1:length(bulks)
                    b1 = bulks[i]
                    b2 = bulks[j]
                    min_id = min(b1.id, b2.id)
                    max_id = max(b1.id, b2.id)

                    # Skip if they are already known to be separated by a cohesive element.
                    (min_id, max_id) in with_cohesive_pairs && continue

                    # Check if they share a face by seeing if they share at least `ndim` nodes.
                    length(intersect(b1.nodes, b2.nodes)) < ndim && continue

                    # Check if they share a face by seeing if they share any facets.
                    if length(intersect(get_facets(b1), get_facets(b2))) > 0
                        push!(bare_pairs, (min_id, max_id))
                    end
                end
            end

            # Iteratively merge the duplicated nodes for the "bare pairs".
            # This loop handles cases where multiple elements need to be joined to the same node.
            changed = true
            while changed
                changed = false
                for (id1, id2) in bare_pairs
                    b1 = model.elems[id1]
                    b2 = model.elems[id2]
                    
                    # Get the nodes at the split location for both elements.
                    pos1 = elem_pos_d[id1]
                    pos2 = elem_pos_d[id2]
                    node1 = b1.nodes[pos1]
                    node2 = b2.nodes[pos2]

                    # Skip if they are already the same node.
                    node1.id == node2.id && continue

                    changed = true

                    # Update the connectivity in both elements to reference the same node.
                    if node1.id > node2.id # preserves the node with higher id
                        b2.nodes[pos2] = node1
                    else
                        b1.nodes[pos1] = node2
                    end
                end
            end
        end

        # ❱❱ 3. Rewire: Update the connectivity of cohesive and special elements.
        for cohe in cohes
            # Mark the cohesive element as being at the leading edge of the crack.
            cohe.cache.open_state = :leading # assume all as leading
            face1_nodes = cohe.nodes[1:m]
            face2_nodes = cohe.nodes[m+1:end]

            # For each face of the cohesive element...
            for (i, face_nodes) in enumerate((face1_nodes, face2_nodes))
                bulk = cohe.couplings[i]  # ...get the corresponding bulk element.
                # Find the (potentially new) node in the bulk element's connectivity.
                pos_b = findfirst( n->hash(n)==hash(node), bulk.nodes )
                bulk_node = bulk.nodes[pos_b]
                # Find the original node in the cohesive element's face.
                pos_c = findfirst(n->hash(n)==hash(node), face_nodes)
                push!(bulk_node.elems, cohe) # Add this cohesive element to the node's affiliations.
                face_nodes[pos_c] = bulk_node # Update the cohesive face to point to the new node.
            end
            # Reassemble the cohesive element's node list.
            cohe.nodes = [ face1_nodes; face2_nodes ]
        end

        # Collect all unique new nodes that were created
        vertex_nodes = [ model.elems[id].nodes[pos] for (id,pos) in elem_pos_d ]
        # unique!(n -> n.id, new_nodes)  # Check this!!! TODO
        
        # Remove any nodes that were not actually duplicated.
        filter!(n -> n.id < 0, vertex_nodes)
        # @show [ n.id for n in vertex_nodes ]
        append!(new_nodes, vertex_nodes)

        # Rebuild the element affiliations for the affected nodes.
        for (id, pos) in elem_pos_d
            node = model.elems[id].nodes[pos]
            push!(node.elems, model.elems[id])
        end

        # Rewire special elements (e.g., line_interfaces) to follow the connectivity of their host bulk element.
        for selem in special_elems
            host = selem.couplings[1]
            n_host_nodes = length(host.nodes)
            selem.nodes[1:n_host_nodes] .= host.nodes
        end

    end # node

    
    # @show length(model.nodes) 
    # @show length(new_nodes) 
    
    # update node numbers
    n_id = length(model.nodes)
    for node in new_nodes
        n_id += 1
        node.id = n_id
    end

    @assert all( n.id>0 for n in new_nodes ) 

    # Add the newly created and numbered nodes to the model.
    append!(model.nodes, new_nodes)

    elem.cache.open_state = :split

    return new_nodes
end


function cracksmesh(mesh::Mesh, opening::Real)

    # Get paired faces
    facedict = Dict{UInt64, Cell}()
    face_pairs = Tuple{Cell, Cell}[]
    for cell in mesh.elems
        for face in get_facets(cell)
            hs = hash(face)
            f  = get(facedict, hs, nothing)
            if f===nothing
                facedict[hs] = face
            else
                push!(face_pairs, (face, f))
                delete!(facedict, hs)
            end
        end
    end

    # Get normals and distances
    U = mesh.node_data["U"]
    crack_faces = Cell[]
    for pair in face_pairs
        face1, face2 = pair

        #error()
        X1 = face1.nodes[1].coord
        X2 = face1.nodes[2].coord
        X3 = face1.nodes[3].coord
        n = cross(X2-X1, X3-X1)
        normalize!(n)
        nnodes = length(face1.nodes)
        node_map = [node.id for node in face1.nodes]
        U1 = mean(U[node_map,:], dims=1)
        node_map = [node.id for node in face2.nodes]
        U2 = mean(U[node_map,:], dims=1)

        #dn = maximum((U2-U1)*n) # normal distance
        dn = dot(U2-U1,n) # normal distance
        #d  = norm(mean(U2-U1, dims=1)) # total distance

        #U2 = mean(U2, dims=1)
        #U1 = mean(U1, dims=1)
        d  = norm(U2-U1) # total distance
        #d = maximum(norm.(eachrow(U2-U1)))
        if dn>0
            #display(U2-U1)
            #error()
        end

        if dn>0 && d>=opening
            push!(crack_faces, face1)
        end
    end

    nodes = get_nodes(crack_faces)
    ids = [ node.id for node in nodes ]

    newsmesh = Mesh(crack_faces)

    for (k,v) in mesh.node_data
        k=="id" && continue
        newsmesh.node_data[k] = v[ids]
    end

    return newsmesh
end
