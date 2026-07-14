# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl


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

    dim_ids         = gmsh.model.occ.getEntities(target_dim)
    host_dimtags    = gmsh.model.getBoundary(dim_ids, false, false, false)
    point_hosts     = Tuple{Int,Vector{Tuple{Int,Int}}}[]
    curve_node_keys = Dict{Int,Set{Tuple{Float64,Float64,Float64}}}()
    node_keys       = Set{Tuple{Float64,Float64,Float64}}()
    nconstraints    = 0

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
    mesh_unstructured(geo, constraint_meshes; ndim=0, recombine=false,
                      quadratic=false, algorithm=:delaunay, sort=true,
                      quiet=false) -> Mesh

Generate the OCC portion of a `GeoModel` with Gmsh. Meshes in
`constraint_meshes` constrain matching OCC boundaries but are not joined into
the returned mesh.

This is an internal generator used by `Mesh(geo)`. Point tags are transferred
to mesh nodes through Gmsh's geometric-entity classification.
"""
function mesh_unstructured(
    geo               ::GeoModel,
    constraint_meshes ::AbstractVector{<:AbstractDomain};
    ndim      ::Int = 0,
    recombine ::Bool = false,
    quadratic ::Bool = false,
    algorithm ::Symbol = :delaunay,
    sort      ::Bool = true,
    quiet     ::Bool = false,
)
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

    gmsh_ndim = gmsh.model.getDimension()
    dim_ids = gmsh.model.occ.getEntities(gmsh_ndim)

    # Set physical groups by tags from synchronized OCC entity ids:
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

    gmsh.model.occ.synchronize()

    embedded_points = [p for p in values(geo.entities) if p isa Point && p.embedded]
    for point in embedded_points
        ent = _get_entity(2, point.coord)
        ent === nothing && error("No surface found at point $(point.coord)")
        gmsh.model.mesh.embed(0, [point.id], ent...)
    end

    gmsh.model.occ.synchronize()

    if !isempty(geo.fields)
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
        gmsh.model.mesh.field.setNumbers(fmin, "FieldsList", field_ids)
        gmsh.model.mesh.field.setAsBackgroundMesh(fmin)
    end

    # Orphan points do not participate in the requested mesh and must not
    # survive solely because they carry metadata.
    for (_, id) in gmsh.model.getEntities(0)
        upward, _ = gmsh.model.getAdjacencies(0, id)
        isempty(upward) && gmsh.model.occ.remove([(0, id)], true)
    end
    gmsh.model.occ.synchronize()

    _add_mesh_constraints!(constraint_meshes, gmsh_ndim)

    logfile = "_gmsh.log"
    try
        open(logfile, "w") do out
            redirect_stdout(out) do
                gmsh.model.mesh.generate(gmsh_ndim)
                quadratic && gmsh.model.mesh.setOrder(2)
                recombine && gmsh.model.mesh.recombine()
            end
        end
    catch
        error("Error generating unstructured mesh.")
    end

    node_ids, coord_list, _ = gmsh.model.mesh.getNodes()
    nnodes = length(node_ids)
    nodes = [Node(coord_list[3*(i-1)+1 : 3*i], id=i) for i in 1:nnodes]
    node_by_gmsh_id = Dict{Int,Node}(Int(node_ids[i]) => nodes[i] for i in 1:nnodes)

    active_point_ids = Set(Int(id) for (_, id) in gmsh.model.getEntities(0))
    for ent in values(geo.entities)
        ent isa Point || continue
        isempty(ent.tag) && continue
        ent.id in active_point_ids || continue

        point_node_ids, _, _ = gmsh.model.mesh.getNodes(0, ent.id)
        for node_id in point_node_ids
            node = get(node_by_gmsh_id, Int(node_id), nothing)
            node === nothing && continue
            isempty(node.tag) && (node.tag = ent.tag)
        end
    end

    shape_dict = Dict(1=>LIN2, 2=>TRI3, 3=>QUAD4, 4=>TET4, 5=>HEX8, 6=>WED6, 7=>PYR5, 8=>LIN3, 9=>TRI6, 10=>QUAD9, 11=>TET10, 12=>HEX27, 16=>QUAD8, 17=>HEX20)
    nnodes_dict = Dict(1=>2, 2=>3, 3=>4, 4=>4, 5=>8, 6=>6, 7=>5, 8=>3, 9=>6, 10=>9, 11=>10, 12=>27, 16=>8, 17=>20)

    cells = Cell[]
    for (_, ent_id) in gmsh.model.getEntities(gmsh_ndim)
        ent = get(geo.entities, (gmsh_ndim, ent_id), nothing)
        ent_tag = ent === nothing ? "" : ent.tag

        elem_types, elem_ids, elem_conns = gmsh.model.mesh.getElements(gmsh_ndim, ent_id)
        for (i, ty) in enumerate(elem_types)
            shape = shape_dict[ty]
            nenodes = nnodes_dict[ty]
            role = ty in (1,8) ? :line : :solid
            for (j, id) in enumerate(elem_ids[i])
                conn = elem_conns[i][(j-1)*nenodes + 1 : j*nenodes]
                if shape == TET10 # Convert Gmsh TET10 ordering to Serendip's FE convention
                    conn[9], conn[10] = conn[10], conn[9] # swap last two nodes
                end
                enodes = Node[node_by_gmsh_id[Int(node_id)] for node_id in conn]
                cell = Cell(shape, role, enodes, tag=ent_tag, id=Int(id))
                push!(cells, cell)
            end
        end
    end

    mesh = Mesh(max(gmsh_ndim, ndim))
    mesh.nodes = isempty(constraint_meshes) ? nodes : get_nodes(cells)
    mesh.elems = cells
    synchronize(mesh, sort=sort)

    return mesh
end
