#❱❱❱ Simple 2D truss example ❰❰❰#

using Serendip

#❱❱❱ Mesh

# Coordinates and connectivity
coord = [ 0 0; 9 0; 18 0; 0 9; 9 9; 18 9.]
conn  = [[1, 2], [1, 5], [2, 3], [2, 6], [2, 5], [2, 4], [3, 6], [3, 5], [4, 5], [5, 6]]

mesh = Mesh(coord, conn, tag="bars")

#❱❱❱ FEM analysis

# Material mapping
mapper = RegionMapper()
add_mapping(mapper, "bars", MechBar, LinearElastic, E=6.894757e7, A=0.043)

# FE model
model = FEModel(mesh, mapper)
ana = MechAnalysis(model, outkey="truss-2d")

stage = add_stage(ana)

add_bc(stage, :node, (x==0, y==0), ux=0, uy=0)
add_bc(stage, :node, (x==0, y==9), ux=0, uy=0)
add_bc(stage, :node, (x==9, y==0), fy=-450.0)
add_bc(stage, :node, (x==18, y==0), fy=-450.0)
add_bc(stage, :body, (x>=0), wy=-10.0)

run(ana)

# ❱❱❱ Post-processing

plot = DomainPlot(model,
    field = "σx´",
    colormap = :coolwarm,
    label = t"$σ_x$ [MPa]",
    warp = 100,
    node_labels = true
)
save(plot, "truss-2d.pdf")