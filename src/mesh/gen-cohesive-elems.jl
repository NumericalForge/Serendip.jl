# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl


# const joint_shape_dict = Dict(
#                         ("LIN2" ,2) => JLIN2,
#                         ("LIN3" ,2) => JLIN3,
#                         ("LIN4" ,2) => JLIN4,
#                         ("TRI3" ,2) => JTRI3,
#                         ("TRI6" ,2) => JTRI6,
#                         ("QUAD4",2) => JQUAD4,
#                         ("QUAD8",2) => JQUAD8,
#                         ("LIN2" ,3) => J3LIN2,
#                         ("LIN3" ,3) => J3LIN3,
#                         ("LIN4" ,3) => J3LIN4,
#                         ("TRI3" ,3) => J3TRI3,
#                         ("TRI6" ,3) => J3TRI6,
#                         ("QUAD4",3) => J3QUAD4,
#                         ("QUAD8",3) => J3QUAD8,
#                        )

"""
    add_boundary_interface_elements(mesh, selector; tag="", support_tag="", quiet=false)

Adds interface elements to the boundary of `mesh`, useful for simulating boundary interactions such as elastic supports in the Winkler foundation model.

This function creates new interface elements at the faces selected by `selector`. Each face gets coupled with duplicated support nodes tagged by `support_tag`. The interface elements are assigned the specified `tag`.

# Arguments
- `mesh::Mesh`: The mesh where interface elements will be added.
- `selector::Union{Expr, Symbolic, String}`: Region selector for boundary faces where interface elements will be created.
- `tag::String=""`: Tag to assign to the generated interface elements.
- `support_tag::String`: Tag for the duplicated support nodes (required).
- `quiet::Bool=false`: Suppress console output if `true`.

# Returns
- `Mesh`: The updated mesh including the new boundary interface elements.

# Example
Use this to model elastic supports on the boundary:
```julia
add_boundary_interface_elements(mesh, selector="boundary_faces", tag="spring", support_tag="foundation")
```
"""
function add_boundary_interface_elements(
    mesh        :: Mesh,
    selector    :: Union{Expr,Symbolic,Tuple,String};
    tag         :: String="",
    support_tag :: String="",
    quiet       :: Bool=false,
)

    quiet || printstyled("Addition of boundary interface elements:\n", bold=true, color=:cyan)
    @check support_tag != "" "add_boundary_interface_elements: tag for support nodes 'support_tag' is required"

    # Get target cells
    faces = mesh.faces[selector]
    @check !isempty(faces) "add_boundary_interface_elements: no target cells found for selector $(repr(selector))"

    # Add interface elements
    joint_cells = Cell[]
    new_nodes_d = Dict{UInt64, Node}()
    for face in faces
        conn = copy(face.nodes)
        for (i, node) in enumerate(face.nodes)
            hs = hash(node)
            n  = get!(new_nodes_d, hs) do
                Node(node.coord, tag=support_tag)
            end
            push!(conn, n)
        end

        jshape = face.shape
        cell = Cell(jshape, conn, tag=tag)
        cell.couplings = [ face.owner ]
        push!(joint_cells, cell)
    end
    new_nodes = collect(values(new_nodes_d))

    # All cells
    append!(mesh.elems, joint_cells)
    append!(mesh.nodes, new_nodes)

    # Update and reorder mesh
    synchronize!(mesh, sort=true, cleandata=true)

    if !quiet
        @printf "  %5d new interface cells\n" length(joint_cells)
        @printf "  %5d total cells\n" length(mesh.elems)
    end

    return mesh
end


function add_boundary_shell_elements(
    mesh          :: Mesh,
    selector      :: Union{Expr,Symbolic,Tuple,String};
    tag           :: String="",
    interface_tag :: String="",
    quiet         :: Bool=false,
)

    quiet || printstyled("Addition of boundary shell elements:\n", bold=true, color=:cyan)
    
    faces = mesh.faces[selector]
    gen_interface_elems = interface_tag != ""

    # Add shell elements
    new_interface_cells   = Cell[]
    new_nodes_d = Dict{UInt64, Node}()
    for face in faces
        conn = copy(face.nodes)
        
        if gen_interface_elems
            conn2 = Node[]
            for (i, node) in enumerate(face.nodes)
                hs = hash(node)
                n  = get!(new_nodes_d, hs) do
                    Node(node.coord)
                end
                push!(conn2, n)
            end
            shell_cell = Cell(face.shape, conn2, tag=tag)
            shell_cell.role = :surface
            push!(new_shell_cells, shell_cell)

            interface_cell = Cell(face.shape, [conn;conn2], tag=interface_tag)
            interface_cell.couplings = [ face.owner, shell_cell ]
            interface_cell.role = :interface
            push!(new_interface_cells, interface_cell)
        else
            shell_cell = Cell(face.shape, conn, tag=tag)
            shell_cell.role = :surface
            push!(new_shell_cells, shell_cell)
        end

    end
    new_nodes = collect(values(new_nodes_d))

    # All cells
    append!(new_shell_cells, shell_cells)
    append!(mesh.elems, new_interface_cells)
    append!(mesh.nodes, new_nodes)

    # Update and reorder mesh
    synchronize!(mesh, sort=true, cleandata=true)

    if !quiet
        @printf "  %4d dimensions                           \n" mesh.ctx.ndim
        @printf "  %5d nodes\n" length(mesh.nodes)
        @printf "  %5d new shell cells\n" length(new_shell_cells)
        length(new_interface_cells)>0 && @printf("  %5d new interface cells\n", length(new_interface_cells))
        @printf "  %5d total cells\n" length(mesh.elems)
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
    inter_regions::Bool=false,
    auto_tag     ::Bool=false,
    layers       ::Int64=2,
    tag          ::String="",
    quiet        ::Bool=false,
    midnodes_tag ::String=""
)

    quiet || printstyled("Addition of cohesive elements:\n", bold=true, color=:cyan)

    layers in (2,3) || error("add_cohesive_elements: wrong number of layers ($layers).")

    inter_regions && tag=="" && (auto_tag=true)

    # Target and locked cells
    # locked cells include solids, lines, beams, etc.
    if selector===nothing
        targetcells = select(mesh.elems, :bulk)
        lockedcells = setdiff(mesh.elems, targetcells)
        # remove previous joints at locked region
        lockedcells = setdiff(lockedcells, select(lockedcells, :interface))
    else
        targetcells = select(mesh.elems, selector, :bulk)
        # targetcells = mesh.elems[selector].solids
        if length(targetcells)==0
            error("add_cohesive_elements: no targetcells found for selector $selector")
        end
        lockedcells = setdiff(mesh.elems, targetcells)
        # remove previous joints at filtered region
        lockedcells = setdiff(lockedcells, select(lockedcells, selector, :interface))
    end

    if !inter_regions
            # Splitting
            # generating new nodes at target cells
            for c in targetcells
                for (i,p) in enumerate(c.nodes)
                    newp = Node(p.coord, tag=p.tag)
                    newp.id = -1
                    c.nodes[i] = newp
                end
            end

            lockedouterfaces = get_outer_facets(lockedcells)

            # Get faces pairs
            facedict = Dict{UInt64, Cell}()
            face_pairs = Tuple{Cell, Cell}[]
            for cell in targetcells
                for face in getfacets(cell)
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

            # Add pairs using surface faces from locked cells
            for face in lockedouterfaces
                hs = hash(face)
                f  = get(facedict, hs, nothing)
                f===nothing && continue
                push!(face_pairs, (face, f))
                delete!(facedict, hs)
            end

    else # generate joints between tagged regions only
        # Get tags
        tag_set = Set{String}()
        for cell in targetcells
            # cell.tag != "" && push!(tag_set, cell.tag)
            push!(tag_set, cell.tag)
        end

        # Get joint faces
        trial_faces = CellFace[]
        for tag in tag_set
            for face in get_outer_facets(targetcells[tag])
                push!(trial_faces, face)
            end
        end

        # @show length(trial_faces)

        # Get nodes to duplicate
        nodes_to_dup = Set{Node}()
        facedict = Dict{UInt64, Cell}()
        ocells   = Cell[] # face owner cells
        for face in trial_faces
            hs = hash(face)
            f  = get(facedict, hs, nothing)
            if f===nothing
                facedict[hs] = face
            else # a matching face was found

                push!(ocells, face.owner)
                push!(ocells, f.owner)
                for node in face.nodes
                    push!(nodes_to_dup, node)
                end
                delete!(facedict, hs)
            end
        end

        # @show length(facedict)
        # @show length(nodes_to_dup)

        # Duplicate nodes
        for cell in targetcells
            for (i,node) in enumerate(cell.nodes)
                if node in nodes_to_dup
                    newnode = Node(node.coord)
                    newnode.id = -1
                    cell.nodes[i] = newnode
                end
            end
        end

        # Join nodes per tag
        for tag in tag_set
            nodedict = Dict{UInt64, Node}()
            for cell in targetcells[tag]
                for (i,node) in enumerate(cell.nodes)
                    hs = hash(node)
                    n  = get(nodedict, hs, nothing)
                    if n===nothing
                        nodedict[hs] = node
                    else
                        cell.nodes[i] = n
                    end
                end
            end
        end

        # Update joint faces (now with new nodes)
        trial_faces = CellFace[]
        for tag in tag_set
            for face in get_outer_facets(ocells[tag])
                push!(trial_faces, face)
            end
        end

        # Get paired faces
        facedict = Dict{UInt64, Cell}()
        face_pairs = Tuple{Cell, Cell}[]
        for face in trial_faces
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

    # Generate joint elements
    jointcells = Cell[]
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

        #jshape = joint_shape(f1.shape)
        # jshape = joint_shape_dict[(f1.shape.name,layers)]
        jshape = f1.shape

        # if tag=="" && auto_tag
        if auto_tag
            tagA, tagB = sort([f1.owner.tag, f2.owner.tag])
            if tagA==tagB
                tag = "joint-"*tagA
            else
                tag = "joint-"*tagA*"-"*tagB
            end
            # tag = join( sort([f1.owner.tag, f2.owner.tag]), "-" )
        end
        cell = Cell(jshape, :interface, con, tag=tag)
        cell.couplings = [f1.owner, f2.owner]
        push!(jointcells, cell)
    end

    # Generate inner nodes at joints (used in hydromechanical analyses)
    if layers==3

        auxnodedict = Dict{UInt64,Node}()

        for jcell in jointcells
            npts = jcell.shape.base_shape.npoints
            sample_pts = jcell.nodes[1:npts]
            for p in sample_pts
                hs = hash(p)
                if haskey(auxnodedict, hs)
                    newp = auxnodedict[hs]
                    push!(jcell.nodes, newp)
                else
                    newp = Node(p.coord.x, p.coord.y, p.coord.z, tag=midnodes_tag)
                    auxnodedict[hs] = newp
                    push!(jcell.nodes, newp)
                end
            end
        end
    end

    # Fix cells connectivities for special interface elements
    for c in lockedcells
        c.role in (:tip, :line_interface) || continue
        scell = c.couplings[1]
        nspts = length(scell.nodes)
        c.nodes[1:nspts] .= scell.nodes
    end

    if haskey(mesh.elem_data, "inset-data")
        idata = mesh.elem_data["inset-data"]
        mesh.elem_data["inset-data"] = [ idata; zeros(Int, length(jointcells), 3) ]
    end

    # All cells
    mesh.elems  = vcat(lockedcells, targetcells, jointcells)

    # Nodes dict
    nodesdict = Dict{Int,Node}()
    idx = length(mesh.nodes)
    # idx = 10000
    # idx = maximum( n.id for n in mesh.nodes )
    for cell in mesh.elems
        for node in cell.nodes
            if node.id<0
                idx += 1
                node.id = idx # new id
            end
            nodesdict[node.id] = node
        end
    end

    # All nodes
    mesh.nodes = collect(values(nodesdict))

    # Update and reorder mesh
    synchronize!(mesh, sort=true, cleandata=true)

    if !quiet
        @printf "  %4dd mesh                             \n" mesh.ctx.ndim
        @printf "  %5d nodes\n" length(mesh.nodes)
        @printf "  %5d total cells\n" length(mesh.elems)
        @printf "  %5d new joint cells\n" length(jointcells)
        nfaces = length(mesh.faces)
        nfaces>0 && @printf("  %5d faces\n", nfaces)
        nedges = length(mesh.edges)
        nedges>0 && @printf("  %5d edges\n", nedges)
    end

    return mesh

end


# function generate_joints_by_tag!(mesh::Mesh; layers::Int64=2, verbose::Bool=true)
#     # Get tags
#     tag_set = Set{String}()
#     for cell in mesh.cells
#         push!(tag_set, cell.tag)
#     end

#     # Get joint faces
#     joint_faces = CellFace[]
#     for tag in tag_set
#         for face in get_outer_facets(mesh.cells[tag])
#             push!(joint_faces, face)
#         end
#     end

#     # Get nodes to duplicate
#     nodes_to_dup = Set{Node}()
#     for face in joint_faces
#         for node in face
#             push!(nodes_to_dup, node)
#         end
#     end

#     # List of owner cells
#     ocells = Cell[ f.owner for f in joint_faces ]

#     # Duplicate nodes
#     for cell in ocells
#         for (i,node) in enumerate(cell.nodes)
#             if node in nodes_to_dup
#                 cell.nodes[i] = copy(node)
#             end
#         end
#     end

#     # Get joint faces (now with new nodes)
#     joint_faces = CellFace[]
#     for tag in tag_set
#         for face in get_outer_facets(ocells[tag])
#             push!(joint_faces, face)
#         end
#     end

#     # Get paired faces
#     facedict = Dict{UInt64, Cell}()
#     face_pairs = Tuple{Cell, Cell}[]
#     for face in joint_faces
#         hs = hash(face)
#         f  = get(facedict, hs, nothing)
#         if f===nothing
#             facedict[hs] = face
#         else
#             push!(face_pairs, (face, f))
#             delete!(facedict, hs)
#         end
#     end

# end


function cracksmesh(mesh::Mesh, opening::Real)

    # Get paired faces
    facedict = Dict{UInt64, Cell}()
    face_pairs = Tuple{Cell, Cell}[]
    for cell in mesh.elems
        for face in getfacets(cell)
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


