# Workflow

This page captures a recommended end-to-end workflow for mechanical analyses in Serendip.

## Typical pipeline

1. Build or import geometry and mesh (`GeoModel`, `Mesh`, mesh transforms).
2. Define mappings from regions to element formulation + constitutive model (`RegionMapper`, `add_mapping`).
3. Construct the finite element model (`FEModel`).
4. Create an analysis (`MechAnalysis` or other analysis type).
5. Add stages, boundary conditions, monitors, and loggers.
6. Run analysis and inspect output files/results.

## Minimal sketch

```julia
using Serendip

# Mesh and mapping
mesh = Mesh(rand_mesh(2, 2))
mapper = RegionMapper()
add_mapping(mapper, :all, MechBulk, LinearElastic; E=30e6, nu=0.2)

# Model + analysis
model = FEModel(mesh, mapper)
ana = MechAnalysis(model)

# Stage and BCs
stage = add_stage(ana)
add_bc(stage, :node, x==0, ux=0, uy=0, uz=0)

# Solve
run(ana)
```

Notes:

1. Use selectors by tag/expr/geometry for robust BC assignment.
2. Prefer `RegionMapper` APIs used in current examples/tests over legacy constructors.
3. Keep scripts under `examples/` as executable references, and mirror them in docs narrative pages.
