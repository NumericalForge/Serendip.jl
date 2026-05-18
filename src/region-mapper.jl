struct RegionMapping
    selector::Any
    etype::Type{<:ElementFormulation}
    cmodel::Type{<:Constitutive}
    quadrature::Union{Int, Tuple}
    params::NamedTuple
    state::NamedTuple
end


function validate_quadrature(quadrature)
    if quadrature isa Int
        quadrature >= 0 || error("add_mapping: quadrature must be >= 0")
        return
    end

    if quadrature isa Tuple
        length(quadrature) in (1, 2, 3) || error("add_mapping: tuple quadrature must have 1, 2, or 3 entries")
        all(x -> x isa Int, quadrature) || error("add_mapping: tuple quadrature entries must be integers")
        return
    end

    error("add_mapping: quadrature must be an integer or a tuple with 1, 2, or 3 entries")
end


"""
    RegionMapper()

Creates an empty `RegionMapper` object.

A `RegionMapper` holds a list of region mappings that associate parts of the mesh (defined by filters) with element formulations, constitutive models, and their parameters.
Mappings can later be added using [`add_mapping`](@ref).

Each mapping associates a filtered region of the mesh with:
- An element formulation (`etype`),
- A constitutive model (`cmodel`),
- A list of parameter values (`params`).
"""
mutable struct RegionMapper
    mappings::Vector{RegionMapping}

    function RegionMapper()
        return new(RegionMapping[])
    end
end


"""
    add_mapping(mapper::RegionMapper, selector, etype, cmodel; quadrature=0, params...)

Adds a new region mapping to the given `RegionMapper`.

Each mapping associates a filtered region of the mesh with:
- An element formulation (`etype`),
- A constitutive model (`cmodel`),
- A list of parameter values (`params`).

# Arguments
- `mapper::RegionMapper`: The mapper to add the region mapping to.
- `selector`: A filtering expression defining the mesh region (e.g., `x==0`, `:all`).
- `etype::Type`: The element formulation type (e.g., `MechSolid`).
- `cmodel::Type`: The constitutive model type (e.g., `LinearElastic`).
- `quadrature=0`: Quadrature request forwarded to each selected element's `set_quadrature` method.
- `params...`: Named parameters for the constitutive model (e.g., `rho=10.0, E=30.0e6`).

# Example
```julia
add_mapping(mapper, x>=0, MechSolid, LinearElastic; rho=10.0, E=30.0e6, nu=0.3)
```

# Trows
An error if a mapping with the same `selector` already exists in the mapper.
"""
function add_mapping(mapper::RegionMapper, selector, etype::Type{S}, cmodel::Type{T}; state::NamedTuple=(;), quadrature::Union{Int,Tuple}=0, params...) where S<:ElementFormulation where T<:Constitutive
    validate_quadrature(quadrature)
    mapping = RegionMapping(selector, etype, cmodel, quadrature, NamedTuple(params), state)
    for m in mapper.mappings
        m.selector == selector && error("Mapping already exists for selector: $selector")
    end
    push!(mapper.mappings, mapping)
end

const add_map = add_mapping  # alias


"""
    RegionModel(etype, cmodel; quadrature=0, params...)

Creates a `RegionMapper` for simple cases where the **same element formulation** and **constitutive model** are applied to the entire mesh.

This is a convenience shortcut equivalent to manually creating a `RegionMapper` and adding a mapping with `selector=:all`.

# Arguments
- `etype::Type`: The element formulation type (e.g., `MechSolid`).
- `cmodel::Type`: The constitutive model type (e.g., `LinearElastic`).
- `quadrature=0`: Quadrature request forwarded to the mapped elements.
- `params...`: Named parameters for the constitutive model.

# Example
```julia
model = RegionModel(MechSolid, LinearElastic; rho=10, E=1.0, nu=0.3)
```
"""
function RegionModel(etype::Type{S}, cmodel::Type{T}; quadrature::Union{Int,Tuple}=0, params...) where S<:ElementFormulation where T<:Constitutive
    mapper = RegionMapper()
    add_mapping(mapper, :all, etype, cmodel; quadrature=quadrature, params...)
    return mapper
end
