mutable struct Stage
    id          ::Int
    name        ::String
    bcs         ::AbstractArray
    nincs       ::Int
    nouts       ::Int
    tspan       ::Float64
    activate  ::Vector{Element}
    deactivate::Vector{Element}
    status      ::Symbol # :idle, :pending, :solving, :done, :error
    analysis    ::Analysis

    function Stage(name::String="";
                    # bcs       ::Vector{BoundaryCondition}=BoundaryCondition[];
                   nincs       ::Int     = 1,
                   nouts       ::Int     = 0,
                   tspan       ::Real  = 0.0,
                   activate  ::Vector{<:Element}=Element[],
                   deactivate::Vector{<:Element}=Element[],
    )
        @check nincs>0
        @check nouts>=0
        bcs = BoundaryCondition[]
        return new(-1, name, bcs, nincs, nouts, tspan, activate, deactivate, :idle)
    end
end