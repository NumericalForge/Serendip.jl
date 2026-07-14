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
- Transfers tags from retained geometric points to their mesh nodes.
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
        occ_mesh = mesh_unstructured(
            geo,
            source_meshes;
            ndim=ndim,
            recombine=recombine,
            quadratic=quadratic,
            algorithm=algorithm,
            sort=sort,
            quiet=quiet,
        )
        mesh = isempty(source_meshes) ? occ_mesh :
               join_meshes(source_meshes..., occ_mesh, check=true)
    else
        if isempty(source_meshes)
            if isempty(geo.gpaths)
                error("Mesh: No blocks, surfaces/volumes, stored meshes, or paths found")
            end

            all(gpath.mode == :free for gpath in geo.gpaths) ||
                error("Mesh: Non-free paths require a base mesh or OCC geometry")

            mesh = Mesh(ndim)
        else
            mesh = has_blocks && length(source_meshes)==1 ? source_meshes[1] :
                   length(source_meshes)==1 ? copy(source_meshes[1]) :
                   join_meshes(source_meshes..., check=true)
        end
        synchronize(mesh, sort=sort)
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
