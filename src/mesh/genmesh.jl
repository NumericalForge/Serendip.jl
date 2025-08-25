

function Mesh(geo::GeoModel;
    ndim      ::Int = 0,
    recombine ::Bool = false,
    quadratic ::Bool = false,
    algorithm ::Symbol = :delaunay,
    sort_nodes::Bool = true,
    quiet     ::Bool = false,
)
    if length(gmsh.model.occ.getEntities(1))>0
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

        length(geo.blocks)>0 && warn("Mesh: Blocks are being ignored")

        # mesh
        gmsh_ndim = gmsh.model.getDimension()
        dim_ids = gmsh.model.occ.getEntities(gmsh_ndim)

        # set physical groups
        for (i, (_, idx)) in enumerate(dim_ids)
            gmsh.model.addPhysicalGroup(gmsh_ndim, [idx], i)
        end

        # make model coherent by fragmentation, removing duplicates and joining points
        # gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()

        # embedded points
        embedded_point = [ p for p in values(geo.added_entities) if p isa Point && p.embedded ]
        for p in embedded_point
            ent = _get_entity(2, p.coord)
            ent === nothing && error("No surface found at point $(p.coord)")
            gmsh.model.mesh.embed(0, [p.id], ent...)
        end

        gmsh.model.occ.synchronize()

        # remove orphan points
        for (_, id) in gmsh.model.getEntities(0)
            # id in embedded_point && continue
            upward, _ = gmsh.model.getAdjacencies(0, id)
            if isempty(upward)
                gmsh.model.occ.remove([(0, id)], true)
            end
        end
        gmsh.model.occ.synchronize()

        tempfile = "_temp.vtk"
        logfile = "_gmsh.log"
        try
            open(logfile, "w") do out
                redirect_stdout(out) do
                    # gmsh.write("_temp.geo_unrolled")
                    gmsh.model.mesh.generate(gmsh_ndim)
                    quadratic && gmsh.model.mesh.setOrder(2)
                    recombine && gmsh.model.mesh.recombine()
                    gmsh.write(tempfile)
                end
            end
        catch err
            error("Error generating unstructured mesh.")
        end

        node_ids, coord_list, _ = gmsh.model.mesh.getNodes()
        nnodes = length(node_ids)
        nodes = [ Node(coord_list[ 3*(i-1)+1 : 3*i ], id=i) for i in 1:nnodes ]

        shape_dict = Dict( 1=>LIN2, 2=>TRI3, 3=>QUAD4, 4=>TET4, 7=>PYR5, 8=>LIN3, 9=>TRI6, 11=>TET10 )
        nnodes_dict = Dict( 1=>2, 2=>3, 3=>4, 4=>4, 7=>5, 8=>3, 9=>6, 11=>10 )

        elem_types, elem_ids, elem_conns = gmsh.model.mesh.getElements(gmsh_ndim)
        cells = Cell[]
        for (i,ty) in enumerate(elem_types)
            shape = shape_dict[ty]
            nenodes = nnodes_dict[ty]
            role = ty in (1,8) ? :line : :bulk
            for (j,id) in enumerate(elem_ids[i])
                conn = elem_conns[i][ (j-1)*nenodes + 1 : j*nenodes ]
                enodes = nodes[conn]
                cell = Cell(shape, role, enodes, id=Int(id))
                push!(cells, cell)
            end
        end

        mesh = Mesh(max(gmsh_ndim, ndim))
        mesh.nodes = nodes
        mesh.elems = cells
        synchronize!(mesh, sortnodes=sort_nodes)

    elseif length(geo.blocks)>0
        quiet || printstyled("Structured mesh generation:\n", bold=true, color=:cyan)
        mesh = mesh_structured(geo, ndim)
    else
        error("Mesh: No blocks or surfaces/volumes found")
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

    # @show mesh.ctx.ndim

    if length(geo.gpaths)>0
        gen_insets!(mesh, geo.gpaths)
        synchronize!(mesh, sortnodes=sort_nodes)
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

