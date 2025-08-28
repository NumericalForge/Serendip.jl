using Serendip

# # Mesh generation
geo = GeoModel()
add_rectangle(geo, [0.0, 0.0], 3.0, 0.4)
mesh = Mesh(geo)
select(mesh, :element, tag="solids")


# Finite element modeling
mapper = RegionMapper()
add_mapping(mapper, "solids", MechBulk, LinearElastic, E=200e6, nu=0.2)

ctx = Context(stress_state=:plane_stress)            
model = FEModel(mesh, mapper, ctx)
ana = MechAnalysis(model)

add_logger(ana, :node, (x==0, y==0), "one-node.dat")
add_logger(ana, :ip, (y<0.025), "ip-list.dat")
add_monitor(ana, :node, (x==3, y==0.4), :uy)

stage = add_stage(ana)
add_bc(stage, :node, (x==0, y==0), ux=0, uy=0)
add_bc(stage, :node, (x==3, y==0), uy=0)
add_bc(stage, :face, (y==0.4), ty=:(-0.1*x)) # triangular load

run(ana)

plot = DomainPlot(model,
    field = "σxx",
    colormap = :coolwarm,
    label = L"\sigma_x",
    warp = 20
)
save(plot, "2d-static.pdf")