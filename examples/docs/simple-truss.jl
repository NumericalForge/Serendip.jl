
# # Simple truss
#
# This page is generated from `examples/docs/simple-truss.jl` using `Literate.jl`.
#
# ## Loading Serendip

using Serendip

# ## Mesh

# Coordinates and connectivity
coord = [ 0 0; 1 0; 1 1; 0 1]
conn  = [[1, 2], [2, 3], [3, 4], [4, 1], [1, 3]]

mesh = Mesh(coord, conn, tag="bars")

# ## FEM analysis

# Material mapping
mapper = RegionMapper()
add_mapping(mapper, "bars", MechBar, LinearElastic, E=10000, A=0.01)

# FE model and analysis
model = FEModel(mesh, mapper)
ana = MechAnalysis(model)

# Analysis stage
stage = add_stage(ana)

# Boundary conditions
add_bc(stage, :node, (x==0, y==0), ux=0, uy=0)
add_bc(stage, :node, (x==0, y==1), ux=0)
add_bc(stage, :node, (x==1, y==1), fy=-1)

# Run analysis
run(ana)

# ## Post-processing

plot = DomainPlot(model,
    field = "σx´",
    colormap = :coolwarm,
    label = t"$σ_x$ [kN]",
    warp = 10
)
save(plot, "simple-truss.svg")

# ## Generated figure
#
# ```@raw html
# <img src="simple-truss.svg" alt="simple-truss plot" width="100%">
# <p><a href="simple-truss.svg">Open/download SVG</a></p>
# ```
