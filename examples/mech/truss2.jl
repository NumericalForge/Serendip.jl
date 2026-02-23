using Serendip

# 2D Truss
coord = [ 0 0; 2 0; 3 0; 0 -2; 2 -2; 3 -2]
conn  = [[1, 2], [2, 3], [4, 5], [5, 6], [1, 4], [2, 5], [3, 6], [2, 4], [3, 5]]

mesh = Mesh(coord, conn, tag="bars")

mapper = RegionMapper()
add_mapping(mapper, "bars", MechBar, LinearElastic, E=1e5, A=1.0)

model = FEModel(mesh, mapper)
ana = MechAnalysis(model)

stage = add_stage(ana)

add_bc(stage, :node, (x==0, y==0), ux=0, uy=0)
add_bc(stage, :node, (x==0, y==-2), ux=0)
add_bc(stage, :node, (x==2, y==0), fy=-1e4)
add_bc(stage, :node, (x==3, y==0), fx=5e3)
add_bc(stage, :node, (x==3, y==-2), fy=-1e4)

run(ana)

# save(model, "truss2.vtu")

plot = DomainPlot(model,
    field = "σx´",
    colormap = :coolwarm,
    label = t"$σ_x$ [kN]",
    warp = 0.2
)
save(plot, "truss.pdf")