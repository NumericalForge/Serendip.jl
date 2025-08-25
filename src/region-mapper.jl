
struct RegionMapping
    selector::Any
    eform::Type{<:ElementFormulation}
    pmodel::Type{<:PhysicsModel}
    params::NamedTuple
end

"""
    RegionMapper()

Creates an empty `RegionMapper` object.

A `RegionMapper` holds a list of region mappings that associate parts of the mesh (defined by filters) with element formulations, physics models, and their parameters.
Mappings can later be added using [`add_mapping`](@ref).
"""
mutable struct RegionMapper
    mappings::Vector{RegionMapping}

    function RegionMapper()
        return new(RegionMapping[])
    end
end


"""
    add_mapping(mapper::RegionMapper, selector, eform, pmodel; params...)

Adds a new region mapping to the given `RegionMapper`.

Each mapping associates a filtered region of the mesh with:
- An element formulation (`eform`),
- A physics model (`pmodel`),
- And user-defined parameter values (`params`).

# Arguments
- `mapper::RegionMapper`: The mapper to add the region mapping to.
- `selector`: A selector expression defining the mesh region (e.g., `x==0`, `:all`).
- `eform::Type`: The element formulation type (e.g., `MechBulk`).
- `pmodel::Type`: The physics model type (e.g., `LinearElastic`).
- `params...`: Named parameters for the physics model (e.g., `rho=10.0, E=30.0e6`).

# Example
```julia
add_mapping(mapper, x>=0, MechBulk, LinearElastic; rho=10.0, E=30.0e6, nu=0.3)
```

# Trows
An error if a mapping with the same `selector` already exists in the mapper.
"""
function add_mapping(mapper::RegionMapper, selector, eform::Type{S}, pmodel::Type{T}; params...) where S<:ElementFormulation where T<:PhysicsModel
    mapping = RegionMapping(selector, eform, pmodel, NamedTuple(params))
    for m in mapper.mappings
        # m.selector == selector && error("Mapping already exists for selector: $selector")
    end
    push!(mapper.mappings, mapping)
end


"""
    RegionModel(eform, pmodel; params...)

Creates a `RegionMapper` for simple cases where the **same element formulation** and **physics model** are applied to the entire mesh.

This is a convenience shortcut equivalent to manually creating a `RegionMapper` and adding a mapping with `selector=:all`.

# Arguments
- `eform::Type`: The element formulation type (e.g., `MechBulk`).
- `pmodel::Type`: The physics model type (e.g., `LinearElastic`).
- `params...`: Named parameters for the physics model.

# Example
```julia
model = RegionModel(MechBulk, LinearElastic; rho=10, E=1.0, nu=0.3)
```
"""
function RegionModel(eform::Type{S}, pmodel::Type{T}; params...) where S<:ElementFormulation where T<:PhysicsModel
    mapper = RegionMapper()
    add_mapping(mapper, :all, eform, pmodel; params...)
    return mapper
end

