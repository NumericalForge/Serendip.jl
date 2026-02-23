# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl


"""
    add_cohesive_elements(mesh, selector=nothing; tag="", implicit=false, quiet=false)

Insert cohesive elements into `mesh` by pairing coincident faces.
Prevents insertion if the interface is already occupied by a contact element.

Algorithm:
1. Split elements into target (where interfaces are searched) and locked (kept as-is).
2. Optionally duplicate target-cell nodes (`implicit=false`) to make interfaces topologically open.
3. Build a geometric index of existing contact interfaces to avoid overlap.
4. Pair coincident faces from target and neighboring locked boundaries.
5. Create cohesive elements with aligned node ordering across both faces.
6. Remove old cohesive elements that match the same owners and interface geometry.

# Arguments
- `mesh::Mesh`: Mesh to modify.
- `selector::Union{Expr,Symbol,Symbolic,Tuple,String,Nothing}=nothing`: Optional selector restricting where cohesive elements are created.
- `tag::String=""`: Tag assigned to generated cohesive elements.
- `implicit::Bool=false`: If `true`, keeps the original shared topology
  (intrinsic/implicit interface activation). If `false`, duplicates target-region
  cell nodes before pairing faces so cohesive elements are created with distinct
  node sets on each side.
- `quiet::Bool=false`: Suppress output.

# Returns
- `Mesh`: The updated mesh including the new cohesive elements.
"""
function add_cohesive_elements(
    mesh         ::Mesh,
    selector     ::Union{Expr,Symbol,Symbolic,Tuple,String,Nothing}=nothing;
    tag          ::String="",
    implicit     ::Bool=false,
    quiet        ::Bool=false,
)

    quiet || printstyled("Addition of cohesive elements:\n", bold=true, color=:cyan)

    # 1. Select target bulk cells and keep the remaining cells as locked
    if selector === nothing
        target_cells = select(mesh.elems, :cont)
        locked_cells = setdiff(mesh.elems, target_cells)
    else
        target_cells = select(mesh.elems, selector, :cont)
        if isempty(target_cells)
            error("add_cohesive_elements: no target cells found for selector $selector")
        end
        locked_cells = setdiff(mesh.elems, target_cells)
    end
    
    # locked_cells may include contact/special elements used later.

    # 1.5 Optional topology split for non-implicit mode (default).
    # Duplicate all nodes of each target cell. This disconnects target cells
    # from neighbors before interface generation.
    if !implicit
        replacements = Dict{AbstractCell, Dict{Node, Node}}()

        for cell in target_cells
            node_map = Dict{Node, Node}()
            for (i, node) in enumerate(cell.nodes)
                if !haskey(node_map, node)
                    new_node = copy(node)
                    new_node.id = -1
                    node_map[node] = new_node
                end
                cell.nodes[i] = node_map[node]
            end
            replacements[cell] = node_map
        end

        # Keep special interface cells consistent with their updated host cell.
        target_set = Set{AbstractCell}(target_cells)
        for cell in locked_cells
            cell.role in (:line_interface, :tip) || continue
            isempty(cell.couplings) && continue
            host = cell.couplings[1]
            host in target_set || continue

            node_map = replacements[host]
            for (i, node) in enumerate(cell.nodes)
                haskey(node_map, node) && (cell.nodes[i] = node_map[node])
            end
        end
    end

    # 2. Index existing contact interfaces by geometry
    _pos_key(n) = (round(n.coord.x, digits=8)+0.0, round(n.coord.y, digits=8)+0.0, round(n.coord.z, digits=8)+0.0)
    
    contact_geo_keys = Set{Vector{Tuple{Float64,Float64,Float64}}}()
    
    # Contacts are non-continuum, so they are expected in locked_cells.
    for cell in locked_cells
        if cell.role == :contact
            # Use the first half (master side) as interface footprint.
            m = div(length(cell.nodes), 2)
            geo_key = sort([_pos_key(n) for n in cell.nodes[1:m]])
            push!(contact_geo_keys, geo_key)
        end
    end

    # 3. Pair coincident faces
    face_pairs = Tuple{Cell, Cell}[] 
    
    # Geometric key -> first unmatched face.
    pending_faces = Dict{Vector{Tuple{Float64,Float64,Float64}}, Cell}()
    
    # A) target internal faces
    for cell in target_cells
        for face in get_facets(cell)
            geo_key = sort([_pos_key(n) for n in face.nodes])
            if haskey(pending_faces, geo_key)
                other_face = pending_faces[geo_key]
                push!(face_pairs, (other_face, face))
                delete!(pending_faces, geo_key)
            else
                pending_faces[geo_key] = face
            end
        end
    end
    
    # B) locked outer faces (interfaces between target and non-target)
    locked_outer_faces = get_outer_facets(collect(locked_cells))
    for face in locked_outer_faces
        geo_key = sort([_pos_key(n) for n in face.nodes])
        if haskey(pending_faces, geo_key)
            other_face = pending_faces[geo_key]
            push!(face_pairs, (other_face, face))
            delete!(pending_faces, geo_key)
        else
            pending_faces[geo_key] = face
        end
    end

    # 4. Build cohesive elements from paired faces
    cohesive_cells = Cell[]
    
    # Signature (owner pair + geometry) used for duplicate cleanup.
    new_signatures = Set{Tuple{Set{Int}, Vector{Tuple{Float64,Float64,Float64}}}}()
    
    for (f1, f2) in face_pairs
        # Skip if a contact element already occupies this interface.
        geo_key = sort([_pos_key(n) for n in f1.nodes])
        if geo_key in contact_geo_keys
            continue
        end
        
        # Build connectivity with face-2 nodes reordered to face-1 sequence.
        n = length(f1.nodes)
        length(f2.nodes) == n || error("add_cohesive_elements: Mismatched face node counts.")

        f2_by_key = Dict{Tuple{Float64,Float64,Float64}, Node}()
        for node in f2.nodes
            f2_by_key[_pos_key(node)] = node
        end

        con = copy(f1.nodes)
        for node in f1.nodes
            node2 = get(f2_by_key, _pos_key(node), nothing)
            node2 === nothing && error("add_cohesive_elements: faces are not coincident.")
            push!(con, node2)
        end

        cell = Cell(f1.shape, :cohesive, con, tag=tag)
        cell.couplings = [f1.owner, f2.owner]
        push!(cohesive_cells, cell)
        
        # Record signature with owner order independence.
        owner_ids = Set([f1.owner.id, f2.owner.id])
        push!(new_signatures, (owner_ids, geo_key))
    end

    # 5. Remove old cohesive elements replaced by the new set
    removed_count = 0
    filter!(locked_cells) do cell
        if cell.role == :cohesive
            if length(cell.couplings) >= 2
                owner_ids = Set([c.id for c in cell.couplings])
                
                # Compare owner pair and first-face geometry.
                m = div(length(cell.nodes), 2)
                geo_key = sort([_pos_key(n) for n in cell.nodes[1:m]])
                
                if (owner_ids, geo_key) in new_signatures
                    removed_count += 1
                    return false
                end
            end
        end
        return true
    end

    # 6. Finalize mesh arrays and synchronize
    append!(mesh.elems, cohesive_cells)
    mesh.elems = [locked_cells; target_cells; cohesive_cells]

    if !implicit
        # Rebuild node list from final connectivity to avoid stale/orphan nodes.
        mesh.nodes = get_nodes(mesh.elems)
    end
    
    synchronize(mesh, sort=true, cleandata=true)

    if !quiet
        @printf "  %5d new cohesive elements\n" length(cohesive_cells)
        if removed_count > 0
            @printf "  %5d existing cohesive elements removed (update)\n" removed_count
        end
        @printf "  %5d total cells\n" length(mesh.elems)
    end

    return mesh
end


export load_cracked_mesh


"""
    load_cracked_mesh(filename::String; tag::String="") -> Mesh

Load a mesh containing crack marker elements and replace them with cohesive elements.

This function reads a mesh from `filename` in which cracks are indicated by
special elements placed along the crack path. These elements serve
only as markers and are removed after cohesive elements are inserted.

Marker elements definition by dimension:
- In 2D, crack paths are marked by line elements located on element edges.
- In 3D, crack surfaces are marked by flat (surface) elements.

# Arguments
- `filename::String`: Path to the input mesh file containing crack markers.

# Keyword Arguments
- `tag::String=""`: Optional tag assigned to the generated cohesive elements.

# Returns
- `Mesh`: A new mesh where crack markers have been replaced by cohesive elements.

# Typical Use Case
Pre-processing of fractured or pre-cracked meshes for cohesive-zone FEM
simulations, where cracks are prescribed geometrically rather than inserted
during analysis.
"""

function load_cracked_mesh(filename::String; tag::String="")
    mesh = Mesh(filename)
    ndim = mesh.ctx.ndim

    role = ndim==2 ? :line : :surface

    markers = select(mesh.elems, role)

    # Iterate over marker elements
    for marker in markers
        neighbors = filter(e -> e.role == :cont, marker.nodes[1].elems)

        # Look for neighboring elements
        for i in 2:length(marker.nodes)
            intersect!(neighbors, filter(e -> e.role == :cont, marker.nodes[i].elems))
        end

        # add cohesive element
        length(neighbors) == 2 || continue

        cohesive = Cell(marker.shape, :cohesive, [marker.nodes; marker.nodes], tag=tag)
        cohesive.couplings = neighbors
        push!(mesh.elems, cohesive)
    end

    # Remove all markers
    mesh.elems = filter(e -> e.role != role, mesh.elems)
    synchronize(mesh, sort=true, cleandata=true)

    return mesh
end




# function cracksmesh(mesh::Mesh, opening::Real)

#     # Get paired faces
#     facedict = Dict{UInt64, Cell}()
#     face_pairs = Tuple{Cell, Cell}[]
#     for cell in mesh.elems
#         for face in get_facets(cell)
#             hs = hash(face)
#             f  = get(facedict, hs, nothing)
#             if f===nothing
#                 facedict[hs] = face
#             else
#                 push!(face_pairs, (face, f))
#                 delete!(facedict, hs)
#             end
#         end
#     end

#     # Get normals and distances
#     U = mesh.node_fields["U"]
#     crack_faces = Cell[]
#     for pair in face_pairs
#         face1, face2 = pair

#         #error()
#         X1 = face1.nodes[1].coord
#         X2 = face1.nodes[2].coord
#         X3 = face1.nodes[3].coord
#         n = cross(X2-X1, X3-X1)
#         normalize!(n)
#         nnodes = length(face1.nodes)
#         node_map = [node.id for node in face1.nodes]
#         U1 = mean(U[node_map,:], dims=1)
#         node_map = [node.id for node in face2.nodes]
#         U2 = mean(U[node_map,:], dims=1)

#         #dn = maximum((U2-U1)*n) # normal distance
#         dn = dot(U2-U1,n) # normal distance
#         #d  = norm(mean(U2-U1, dims=1)) # total distance

#         #U2 = mean(U2, dims=1)
#         #U1 = mean(U1, dims=1)
#         d  = norm(U2-U1) # total distance
#         #d = maximum(norm.(eachrow(U2-U1)))
#         if dn>0
#             #display(U2-U1)
#             #error()
#         end

#         if dn>0 && d>=opening
#             push!(crack_faces, face1)
#         end
#     end

#     nodes = get_nodes(crack_faces)
#     ids = [ node.id for node in nodes ]

#     newsmesh = Mesh(crack_faces)

#     for (k,v) in mesh.node_fields
#         k=="id" && continue
#         newsmesh.node_fields[k] = v[ids]
#     end

#     return newsmesh
# end
