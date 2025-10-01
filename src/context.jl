

"""
    Context(; ndim=3, stress_state=:auto, transient=nothing, thickness=1.0, g=0.0, T0=0.0)

Defines global analysis metadata and configuration parameters for the finite element model.

This struct encapsulates the physical and numerical context in which the FE model is defined, including dimensionality, stress state, analysis type, and environmental parameters such as gravity and reference temperature.

# Keyword Arguments
- `ndim::Int=3`:
  Number of spatial dimensions.
  - `3` for 3D analysis.
  - `2` for 2D analysis.

- `stress_state::Symbol=:auto`:
  Stress state for 2D problems.
  Accepted values:
  - `:auto` — Automatically determined from geometry and settings.
  - `:plane_stress` — Plane stress assumption (e.g., thin plates).
  - `:plane_strain` — Plane strain assumption (e.g., long tunnels).
  - `:axisymmetric` — Axisymmetric assumption for rotationally symmetric models.

- `transient::Flag=nothing`:
  Flag for transient analysis.
  Accepted values:
  - `true` — Transient (time-dependent) analysis.
  - `false` — Steady-state analysis.
  - `nothing` — Auto-detect based on problem setup.

- `thickness::Float64=1.0`:
  Thickness for 2D analyses. Ignored for 3D problems.

- `g::Float64=0.0`:
  Gravity acceleration.

- `T0::Float64=0.0`:
  Reference temperature for thermal effects (if applicable).

# Example
```julia
ctx = Context(ndim=2, stress_state=:plane_strain, thickness=0.1, g=9.81)
```
"""
mutable struct Context
    ndim::Int
    stress_state::Symbol
    transient::Flag
    thickness::Float64
    g::Float64
    T0::Float64

    function Context(;
        ndim::Int            = 3,
        stress_state::Symbol = :auto,
        transient::Flag      = nothing,
        thickness::Float64   = 1.0,
        g::Float64           = 0.0,
        T0::Float64          = 0.0,
    )
        @check ndim == 2 || ndim == 3 "Context: ndim must be 2 or 3."
        @check stress_state in (:auto, :plane_stress, :plane_strain, :axisymmetric)
        return new(ndim, stress_state, transient, thickness, g, T0)
    end
end
