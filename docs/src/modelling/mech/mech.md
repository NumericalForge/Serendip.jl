
# Mechanical Analysis

## Elements

```@docs
MechCont
MechBar
MechBeam
MechShell
MechMembrane
MechLSJoint
MechContact
```

## Constitutive models

```@docs
LinearElastic
DruckerPrager
VonMises
TensionFieldSoftening
WillamWarnke
```

## Membranes

`MechMembrane` uses three translations per node and a fixed reference
surface in a 3D analysis. It accepts 2D area shapes with implemented
interpolation, derivatives, and surface quadrature, with positive constant
thickness. Patch, rotation, load, and mass tests cover `Tri3`, `Tri6`,
`Quad4`, `Quad8`, and `Quad9`. There is no bending or geometric stiffness;
a flat membrane needs transverse restraint from supports, solids, or interfaces.

```julia
# LinearElastic retains in-plane compressive resistance.
elastic = RegionModel(MechMembrane, LinearElastic,
    E=200.0, nu=0.3, thickness=0.01)

# Omit GF for perfect plasticity; provide GF for power-law softening.
tension = RegionModel(MechMembrane, TensionFieldSoftening,
    E=200.0, nu=0.3, fy=2.0, GF=0.02, beta=2.0, thickness=0.01)
```

`TensionFieldSoftening` relaxes compressive stresses and preserves plastic history
through reversible wrinkling. `GF` is plastic dissipation per rupture area;
the critical accumulated plastic strain is calculated using the element's
reference-area length `h = sqrt(Ae)`. This length is approximate for elongated
elements. All constitutive vectors use local in-plane Mandel components.

A ruptured material point returns zero stress and tangent permanently. Failed
local iterations instead request global cutback (`autoinc=true`); they do not
mark rupture. Steep softening can prevent a continuous strain-controlled
solution even with cutbacks. The initial tangent uses the elastic tensile-side
predictor, while accepted slack and ruptured states have exactly zero tangent.

Boundary cells created with `add_boundary_shell_elements` have role `:surface`
and can be mapped to `MechMembrane`; the mesh helper does not prescribe the
mechanical formulation.

## Surface stress recovery

Membranes and shells provide the nodal fields `σ1′`, `σ2′`, and `σvm′`.
The first two are the ordered in-plane principal stresses; the von Mises
stress includes the full constitutive stress state. Materials that store
accumulated equivalent plastic strain also provide `εp′`. The prime `′`
identifies a surface result, so these fields remain separate from solid
fields at shared nodes. For example:

```julia
mplot = DomainPlot()
add_plot(mplot, model; field="σ1′")
```

Surface patches use a best-fit plane from the midsurface nodes. Samples and
nodes are projected into centered, scaled coordinates for linear
regression using the complete basis `(1, u, v)`. A constant fit is used when
there are fewer than three samples, the sample layout has insufficient rank,
or the plane is ambiguous. Predictions from overlapping patches are averaged into one
value per node; edges and corners therefore show a blended estimate.
Regression is unweighted and does not clip recovered values.

At each shell surface quadrature position, all recovered fields come from
the thickness layer with greatest von Mises stress. Exact ties retain the
first layer in quadrature order. The sample is located on the midsurface
for recovery, independently of the selected layer's thickness offset.

Continuous constitutive fields are no longer exported as cell averages or
extremes for these formulations. The original integration-point fields
remain available for inspection and loggers, including the existing `εp`
output of `VonMises`. For membrane materials with a rupture flag, the cell
field `ruptured` remains the area-weighted fraction of ruptured integration
points; the Boolean history flag itself is not regressed. Principal strains
are not included in surface recovery.
