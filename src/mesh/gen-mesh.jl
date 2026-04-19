function _candidate_constraint_hosts(node::Node, host_dimtags; tol=1e-8)
    X = node.coord
    hosts = Tuple{Int,Int}[]

    for (dim, id) in host_dimtags
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(dim, id)
        if xmin - tol <= X.x <= xmax + tol &&
           ymin - tol <= X.y <= ymax + tol &&
           zmin - tol <= X.z <= zmax + tol
            push!(hosts, (Int(dim), Int(id)))
        end
    end

    return hosts
end


function _add_mesh_constraints!(meshes, target_dim::Integer)
    target_dim in (2, 3) || return 0

    dim_ids = gmsh.model.occ.getEntities(target_dim)
    host_dimtags = gmsh.model.getBoundary(dim_ids, false, false, false)
    point_hosts = Tuple{Int,Vector{Tuple{Int,Int}}}[]
    curve_node_keys = Dict{Int,Set{Tuple{Float64,Float64,Float64}}}()
    node_keys = Set{Tuple{Float64,Float64,Float64}}()
    nconstraints = 0

    for mesh in meshes
        facets = get_outer_facets(select(mesh.elems, :solid))
        for facet in facets
            facet.shape.ndim == target_dim - 1 || continue

            nvertices = facet.shape.base_shape.npoints
            for node in facet.nodes[1:nvertices]
                key = node_pos_key(node)
                key in node_keys && continue
                push!(node_keys, key)

                hosts = _candidate_constraint_hosts(node, host_dimtags)
                isempty(hosts) && continue

                for (dim, host) in hosts
                    dim == 1 || continue
                    keys = get!(curve_node_keys, host, Set{Tuple{Float64,Float64,Float64}}())
                    push!(keys, key)
                end

                tag = gmsh.model.occ.addPoint(node.coord.x, node.coord.y, node.coord.z)
                push!(point_hosts, (tag, hosts))
            end
        end
    end

    isempty(point_hosts) && return 0

    gmsh.model.occ.synchronize()

    # When stored mesh nodes constrain a boundary curve, keep Gmsh from adding
    # extra curve nodes that would later become hanging nodes at the interface.
    if target_dim == 2
        for (curve, keys) in curve_node_keys
            length(keys) >= 2 || continue
            gmsh.model.mesh.set_transfinite_curve(curve, length(keys))
        end
    end

    for (tag, hosts) in point_hosts
        for (dim, host) in hosts
            gmsh.model.mesh.embed(0, [tag], dim, host)
            nconstraints += 1
        end
    end

    return nconstraints
end


"""
    Mesh(geo::GeoModel; ndim=0, recombine=false, quadratic=false, algorithm=:delaunay,
         sort=true, quiet=false) -> Mesh

Generate a finite element mesh from a geometric model using Gmsh or structured
blocks. Supports both **unstructured** meshing from OCC entities and **structured**
meshing from `geo.blocks`. If both are present, block boundaries constrain the
OCC mesh and both meshes are joined.

# Arguments
- `geo::GeoModel`: Geometry model containing OCC entities and/or structured blocks.
- `ndim::Int=0`: Target mesh dimension (1, 2, or 3). Defaults to the maximum dimension
  present in the geometry.
- `recombine::Bool=false`: If `true`, recombines simplices (triangles, tetrahedra)
  into quads/hexas when possible.
- `quadratic::Bool=false`: If `true`, generates quadratic (second-order) elements.
- `algorithm::Symbol=:delaunay`: Meshing algorithm:
    - `:delaunay` → Delaunay (surface and 3D).
    - `:mesh_adapt` → Adaptive meshing.
    - `:best` → Frontal Delaunay (surface) + Delaunay (3D).
    - `:frontal` → Frontal Delaunay (surface and 3D).
- `sort::Bool=true`: If `true`, reorder mesh entities (nodes, elements)
  consistently after generation.
- `quiet::Bool=false`: If `true`, suppress console messages.

# Behavior
- Sets Gmsh meshing options according to algorithm and element order.
- Defines physical groups for all top-dimensional entities.
- Handles embedded points (`p.embedded == true`) by embedding them in host surfaces.
- Removes orphan vertices without adjacencies.
- Builds `Node` and `Cell` objects from Gmsh mesh data, including element connectivity.
- For structured meshes, uses block definitions from `geo`.
- Existing meshes added with `add_mesh(geo, mesh)` are included in the result.
- If structured blocks and/or stored meshes are combined with OCC entities, their
  boundary vertex nodes are used as Gmsh constraints before meshing the OCC region.
- After all parts are generated, meshes are joined with node welding and
  conformity checks enabled.
- If `geo.gpaths` are present, generates insets and re-synchronizes.

# Notes
- For quadratic stored or block meshes, only boundary vertex nodes are used as
  Gmsh constraints. High-order midside nodes are checked later during the final
  mesh join.

# Returns
- `mesh::Mesh`: A mesh object containing nodes, elements, faces, edges, and context information.

# Example
```julia
geo = GeoModel()
# ... build OCC surfaces/volumes or add blocks ...
mesh = Mesh(geo; ndim=3, quadratic=true, algorithm=:frontal, sort=true)

```
"""
function Mesh(geo::GeoModel;
    ndim      ::Int = 0,
    recombine ::Bool = false,
    quadratic ::Bool = false,
    algorithm ::Symbol = :delaunay,
    sort      ::Bool = true,
    quiet     ::Bool = false,
)

    source_meshes = AbstractDomain[geo.meshes...]
    has_blocks = length(geo.blocks)>0
    has_occ = length(gmsh.model.occ.getEntities(1))>0

    if has_blocks
        quiet || printstyled("Structured mesh generation:\n", bold=true, color=:cyan)
        block_mesh = mesh_structured(geo, ndim)
        push!(source_meshes, block_mesh)
    end

    if has_occ
        quiet || printstyled("Unstructured mesh generation:\n", bold=true, color=:cyan)
        quadratic && gmsh.option.setNumber("Mesh.ElementOrder", 2)

        if algorithm==:delaunay
            gmsh.option.setNumber("Mesh.Algorithm", 5) # delaunay
            gmsh.option.setNumber("Mesh.Algorithm3D", 1) # delaunay
        elseif algorithm==:mesh_adapt
            gmsh.option.setNumber("Mesh.Algorithm", 1) # mesh adapt
            gmsh.option.setNumber("Mesh.Algorithm3D", 1) # delaunay
        elseif algorithm==:best
            gmsh.option.setNumber("Mesh.Algorithm", 6) # frontal delaunay
            gmsh.option.setNumber("Mesh.Algorithm3D", 1) # delaunay
        elseif algorithm==:frontal
            gmsh.option.setNumber("Mesh.Algorithm", 6) # frontal delaunay
            gmsh.option.setNumber("Mesh.Algorithm3D", 4) # frontal
        end

        # mesh
        gmsh_ndim = gmsh.model.getDimension()
        dim_ids = gmsh.model.occ.getEntities(gmsh_ndim)

        # set physical groups by tags from synchronized OCC entity ids:
        # - tagged entities are grouped together by tag
        # - untagged entities keep one physical group per entity
        tag_to_ids = Dict{String, Vector{Int}}()
        ordered_tags = String[]
        untagged_ids = Int[]

        for (_, id) in dim_ids
            ent = get(geo.entities, (gmsh_ndim, id), nothing)
            if ent === nothing || ent.tag == ""
                push!(untagged_ids, id)
                continue
            end

            if !haskey(tag_to_ids, ent.tag)
                tag_to_ids[ent.tag] = Int[]
                push!(ordered_tags, ent.tag)
            end
            push!(tag_to_ids[ent.tag], id)
        end

        gidx = 1
        for tag in ordered_tags
            gmsh.model.addPhysicalGroup(gmsh_ndim, tag_to_ids[tag], gidx)
            gmsh.model.setPhysicalName(gmsh_ndim, gidx, tag)
            gidx += 1
        end

        for id in untagged_ids
            gmsh.model.addPhysicalGroup(gmsh_ndim, [id], gidx)
            gidx += 1
        end

        # make model coherent by fragmentation, removing duplicates and joining points
        # gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()

        # embedded points
        embedded_point = [ p for p in values(geo.entities) if p isa Point && p.embedded ]
        for p in embedded_point
            ent = _get_entity(2, p.coord)
            ent === nothing && error("No surface found at point $(p.coord)")
            gmsh.model.mesh.embed(0, [p.id], ent...)
        end

        gmsh.model.occ.synchronize()

        # refinement fields
        if length(geo.fields)>0
            # see https://gmsh.info/doc/texinfo/gmsh.html#Gmsh-mesh-size-fields
            field_ids = Int[]
            for field in geo.fields
                size1, size2 = field.size1, field.size2
                a, b, c = field.rx, field.ry, field.rz
                x, y, z = field.coord

                field_id = gmsh.model.mesh.field.add("MathEval")
                n = 2/max(field.roundness,0.01) # roundness 1 -> n=2 (ellipsoid), roundness 0 -> n=∞ (box)
                g = field.transition # 0 -> sharp, 1 -> linear

                an, bn, cn = a^n, b^n, c^n
                expr = "$size2 + ($size1 - $size2) * abs( 1 - (abs(x-$x)^$n/$an + abs(y-$y)^$n/$bn + abs(z-$z)^$n/$cn)^(1/$n))^$g * step(1 - ( abs(x-$x)^$n/$an + abs(y-$y)^$n/$bn + abs(z-$z)^$n/$cn ))"

                gmsh.model.mesh.field.setString(field_id, "F", expr)
                push!(field_ids, field_id)
            end

            fmin = gmsh.model.mesh.field.add("Min")
            gmsh.model.mesh.field.setNumbers(fmin,"FieldsList", field_ids)
            gmsh.model.mesh.field.setAsBackgroundMesh(fmin)
        end

        # remove orphan points
        for (_, id) in gmsh.model.getEntities(0)
            # id in embedded_point && continue
            upward, _ = gmsh.model.getAdjacencies(0, id)
            if isempty(upward)
                gmsh.model.occ.remove([(0, id)], true)
            end
        end
        gmsh.model.occ.synchronize()

        _add_mesh_constraints!(source_meshes, gmsh_ndim)

        # ❱❱❱ Mesh generation
        logfile = "_gmsh.log"
        try
            open(logfile, "w") do out
                redirect_stdout(out) do
                    # gmsh.write("_temp.geo_unrolled")
                    gmsh.model.mesh.generate(gmsh_ndim) # mesh generation
                    quadratic && gmsh.model.mesh.setOrder(2)
                    recombine && gmsh.model.mesh.recombine()
                    # gmsh.write("_temp.vtk")
                end
            end
        catch err
            error("Error generating unstructured mesh.")
        end

        node_ids, coord_list, _ = gmsh.model.mesh.getNodes()
        nnodes = length(node_ids)
        nodes = [ Node(coord_list[ 3*(i-1)+1 : 3*i ], id=i) for i in 1:nnodes ]
        node_by_tag = Dict{Int, Node}(Int(node_ids[i]) => nodes[i] for i in 1:nnodes)

        shape_dict = Dict( 1=>LIN2, 2=>TRI3, 3=>QUAD4, 4=>TET4, 5=>HEX8, 6=>WED6, 7=>PYR5, 8=>LIN3, 9=>TRI6, 10=>QUAD9, 11=>TET10, 12=>HEX27, 16=>QUAD8, 17=>HEX20 )
        nnodes_dict = Dict( 1=>2, 2=>3, 3=>4, 4=>4, 5=>8, 6=>6, 7=>5, 8=>3, 9=>6, 10=>9, 11=>10, 12=>27, 16=>8, 17=>20 )

        cells = Cell[]
        for (_, ent_id) in gmsh.model.getEntities(gmsh_ndim)
            ent = get(geo.entities, (gmsh_ndim, ent_id), nothing)
            ent_tag = ent === nothing ? "" : ent.tag

            elem_types, elem_ids, elem_conns = gmsh.model.mesh.getElements(gmsh_ndim, ent_id)
            for (i, ty) in enumerate(elem_types)
                shape = shape_dict[ty]
                nenodes = nnodes_dict[ty]
                role = ty in (1,8) ? :line : :solid
                for (j,id) in enumerate(elem_ids[i])
                    conn = elem_conns[i][ (j-1)*nenodes + 1 : j*nenodes ]
                    if shape == TET10 # Convert Gmsh TET10 ordering to Serendip's FE convention
                        conn[9], conn[10] = conn[10], conn[9] # swap last two nodes
                    end
                    enodes = Node[node_by_tag[Int(tag)] for tag in conn]
                    cell = Cell(shape, role, enodes, tag=ent_tag, id=Int(id))
                    push!(cells, cell)
                end
            end
        end

        mesh = Mesh(max(gmsh_ndim, ndim))
        mesh.nodes = isempty(source_meshes) ? nodes : get_nodes(cells)
        mesh.elems = cells
        synchronize(mesh, sort=sort)

        if !isempty(source_meshes)
            mesh = join_meshes(source_meshes..., mesh, check=true)
        end

    else
        if isempty(source_meshes)
            error("Mesh: No blocks, surfaces/volumes, or stored meshes found")
        end
        mesh = has_blocks && length(source_meshes)==1 ? source_meshes[1] :
               length(source_meshes)==1 ? copy(source_meshes[1]) :
               join_meshes(source_meshes..., check=true)
    end

    mesh.ctx.ndim = max(mesh.ctx.ndim, ndim)

    # fix role for surface elements
    if mesh.ctx.ndim == 3
        for elem in mesh.elems
            if elem.shape.ndim == 2
                elem.role = :surface
            end
        end
    end

    if length(geo.gpaths)>0
        gen_path_insets(mesh, geo.gpaths)
        synchronize(mesh, sort=sort)
    end

    if !quiet
        npoints = length(mesh.nodes)
        ncells  = length(mesh.elems)
        println("  $(mesh.ctx.ndim)d Mesh:")
        # @printf "  %4dd mesh\033[K\n" mesh.ctx.ndim
        @printf "  %5d nodes\n" npoints
        @printf "  %5d cells\n" ncells
        nfaces  = length(mesh.faces)
        @printf "  %5d faces\n" nfaces
        nedges  = length(mesh.edges)
        nedges>0 && @printf "  %5d surface edges\n" nedges
    end

    return mesh

end
