mutable struct Stage
    id          ::Int
    name        ::String
    bcs         ::Vector{BoundaryCondition}
    constraints ::Vector{Constraint}
    nincs       ::Int
    nouts       ::Int
    tspan       ::Float64
    activate    ::Vector{Element}
    deactivate  ::Vector{Element}
    status      ::Symbol # :idle, :pending, :solving, :done, :error
    analysis    ::Analysis

    function Stage(name::String="";
                   nincs       ::Int     = 1,
                   nouts       ::Int     = 0,
                   tspan       ::Real  = 0.0,
                   activate  ::Vector{<:Element}=Element[],
                   deactivate::Vector{<:Element}=Element[],
    )
        @check nincs>0
        @check nouts>=0
        bcs = BoundaryCondition[]
        constraints = Constraint[]
        return new(-1, name, bcs, constraints, nincs, nouts, tspan, activate, deactivate, :idle)
    end
end


function _resolve_stage_target(stage::Stage, kind::Symbol, selector; prefix::String)
    kind in (:node, :face, :edge, :body) || error("$prefix: Invalid condition kind: $kind. Use :node, :face, :edge or :body.")

    model = stage.analysis.model
    target_type = kind == :node ? Node : kind == :body ? AbstractCell : kind == :face ? CellFace : CellEdge
    item_name = kind == :node ? :node : kind == :body ? :element : kind == :face ? :face : :edge
    items = kind == :node ? model.nodes : kind == :body ? model.elems : kind == :face ? model.faces : model.edges

    if kind == :node && selector isa AbstractArray
        X = Vec3(selector)
        x, y, z = X
        selector = :(x==$x && y==$y && z==$z)
    end

    target = select(items, selector)
    length(target) == 0 && alert("$prefix: No $(item_name)s found for selector: ", selector)

    return target_type, selector, target
end


"""
    add_bc(stage::Stage, kind::Symbol, selector; conds...)

Add a boundary condition (BC) to the given `stage` in the analysis.

This function attaches a boundary condition to a specific set of entities in the model, identified by `kind` and filtered using a spatial expression or coordinates.

# Arguments
- `stage::Stage`: The analysis stage where the boundary condition will be applied.
- `kind::Symbol`: The type of boundary entity to apply the condition on. Options:
    - `:node` – apply on nodes.
    - `:face` – apply on faces (in 3D) or edges (in 2D).
    - `:edge` – apply on edges explicitly.
- `selector`: A filtering expression or array of coordinates to select entities for the BC. If a coordinate array is provided for `:node`, it is converted into an exact coordinate-equality selector.
- `conds...`: Named keyword arguments specifying the boundary conditions to apply (e.g., `ux=0`, `uy=0`, `tz=-5`).

# Behavior
- Resolves the target entities (`nodes`, `faces`, or `edges`) in the finite element model using the provided `selector`.
- For `kind == :node`, an array selector like `[x, y, z]` is treated as an exact point lookup, not a nearest search.
- If `kind == :face` in a 2D model, surface BCs are automatically mapped to edges with a notification.
- Adds the resulting `BoundaryCondition` to `stage.bcs`.

# Returns
- `bc::BoundaryCondition`: The created boundary condition object.

# Examples
```julia
add_bc(stage, :node, x==0, ux=0, uy=0)      # Fix displacement on nodes at x==0
add_bc(stage, :face, z==1, tz=-10)       # Apply surface traction at z==1
add_bc(stage, :edge, (y==0,z==1), qx=10)    # Apply linear traction at edges where y==0 && z==1
add_bc(stage, :node, [0.0, 0.0, 0.0], ux=0) # Constrain node at origin
```
"""
function add_bc(
    stage::Stage,
    kind::Symbol,
    selector;
    conds...
    )

    target_type, selector, target = _resolve_stage_target(stage, kind, selector; prefix="add_bc")

    bc = BoundaryCondition{target_type}(kind, selector, conds, target)
    push!(stage.bcs, bc)

    return bc
end


function add_constraint(stage::Stage, kind::Symbol, selector, expr::Union{Expr,Symbolic})
    kind in (:node, :face, :edge) || error("add_constraint: Invalid constraint kind: $kind. Use :node, :face or :edge.")

    target_type, selector, target = _resolve_stage_target(stage, kind, selector; prefix="add_constraint")
    cexpr = getexpr(expr)
    cexpr isa Expr || error("add_constraint: Constraint expression must be an Expr equality.")
    terms, rhs = get_affine_terms(cexpr)

    constraint = Constraint{target_type}(kind, selector, cexpr, target, terms, rhs)
    push!(stage.constraints, constraint)

    return constraint
end


# function add_body_load(stage::Stage, selector; conds...)  # Body load treated as a special case of boundary condition
#     model = stage.analysis.model
#     target = model.elems.active[selector]
#     length(target) == 0 && notify("add_body_load: No elements found for selector: ", selector)

#     return BoundaryCondition(:body, selector, conds, target)
# end
