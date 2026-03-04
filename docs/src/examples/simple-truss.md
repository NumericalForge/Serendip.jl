```@meta
EditURL = "../../../examples/docs/simple-truss.jl"
```

# Simple truss

This page is generated from `examples/docs/simple-truss.jl` using `Literate.jl`.

## Loading Serendip

````@example simple-truss
using Serendip
````

## Mesh

Coordinates and connectivity

````@example simple-truss
coord = [ 0 0; 1 0; 1 1; 0 1]
conn  = [[1, 2], [2, 3], [3, 4], [4, 1], [1, 3]]

mesh = Mesh(coord, conn, tag="bars")
````

## FEM analysis

Material mapping

````@example simple-truss
mapper = RegionMapper()
add_mapping(mapper, "bars", MechBar, LinearElastic, E=10000, A=0.01)
````

FE model and analysis

````@example simple-truss
model = FEModel(mesh, mapper)
ana = MechAnalysis(model)
````

Analysis stage

````@example simple-truss
stage = add_stage(ana)
````

Boundary conditions

````@example simple-truss
add_bc(stage, :node, (x==0, y==0), ux=0, uy=0)
add_bc(stage, :node, (x==0, y==1), ux=0)
add_bc(stage, :node, (x==1, y==1), fy=-1)
````

Run analysis

````@example simple-truss
run(ana)
````

## Post-processing

````@example simple-truss
plot = DomainPlot(model,
    field = "σx´",
    colormap = :coolwarm,
    label = t"$σ_x$ [kN]",
    warp = 10
)
save(plot, "simple-truss.svg")
````

## Generated figure

```@raw html
<img src="simple-truss.svg" alt="simple-truss plot" width="100%">
<p><a href="simple-truss.svg">Open/download SVG</a></p>
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

