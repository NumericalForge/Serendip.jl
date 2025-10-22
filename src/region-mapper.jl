
struct RegionMapping
    selector::Any
    etype::Type{<:ElementFormulation}
    cmodel::Type{<:Constitutive}
    params::NamedTuple
    state::NamedTuple
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
    add_mapping(mapper::RegionMapper, selector, etype, cmodel; params...)

Adds a new region mapping to the given `RegionMapper`.

Each mapping associates a filtered region of the mesh with:
- An element formulation (`etype`),
- A constitutive model (`cmodel`),
- A list of parameter values (`params`).

# Arguments
- `mapper::RegionMapper`: The mapper to add the region mapping to.
- `selector`: A filtering expression defining the mesh region (e.g., `x==0`, `:all`).
- `etype::Type`: The element formulation type (e.g., `MechBulk`).
- `cmodel::Type`: The constitutive model type (e.g., `LinearElastic`).
- `params...`: Named parameters for the constitutive model (e.g., `rho=10.0, E=30.0e6`).

# Example
```julia
add_mapping(mapper, x>=0, MechBulk, LinearElastic; rho=10.0, E=30.0e6, nu=0.3)
```

# Trows
An error if a mapping with the same `selector` already exists in the mapper.
"""
function add_mapping(mapper::RegionMapper, selector, etype::Type{S}, cmodel::Type{T}; state::NamedTuple=(;), params...) where S<:ElementFormulation where T<:Constitutive
    mapping = RegionMapping(selector, etype, cmodel, NamedTuple(params), state)
    for m in mapper.mappings
        m.selector == selector && error("Mapping already exists for selector: $selector")
    end
    push!(mapper.mappings, mapping)
end


"""
    RegionModel(etype, cmodel; params...)

Creates a `RegionMapper` for simple cases where the **same element formulation** and **constitutive model** are applied to the entire mesh.

This is a convenience shortcut equivalent to manually creating a `RegionMapper` and adding a mapping with `selector=:all`.

# Arguments
- `etype::Type`: The element formulation type (e.g., `MechBulk`).
- `cmodel::Type`: The constitutive model type (e.g., `LinearElastic`).
- `params...`: Named parameters for the constitutive model.

# Example
```julia
model = RegionModel(MechBulk, LinearElastic; rho=10, E=1.0, nu=0.3)
```
"""
function RegionModel(etype::Type{S}, cmodel::Type{T}; params...) where S<:ElementFormulation where T<:Constitutive
    mapper = RegionMapper()
    add_mapping(mapper, :all, etype, cmodel; params...)
    return mapper
end

