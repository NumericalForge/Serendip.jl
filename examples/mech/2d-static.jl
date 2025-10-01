using Serendip

# ❱❱❱ Geometry and mesh generation

geo = GeoModel()
add_rectangle(geo, [0.0, 0.0], 3.0, 0.4) # geo, corner, width, height

mesh = Mesh(geo)
select(mesh, :element, tag="solids") # tag all elements as "solids"

# ❱❱❱ Finite element modeling

mapper = RegionMapper() 
add_mapping(mapper, "solids", MechBulk, LinearElastic, E=200e6, nu=0.2)

model = FEModel(mesh, mapper, stress_state=:plane_stress)
ana = MechAnalysis(model)

add_logger(ana, :node, (x==0, y==0), "one-node.dat") # analysis, kind, filter, filename
add_logger(ana, :ip, (y<0.025), "ip-list.dat")

add_monitor(ana, :node, (x==3, y==0.4), :uy) # analysis, kind, filter, variable

stage = add_stage(ana)
add_bc(stage, :node, (x==0, y==0), ux=0, uy=0) # stage, kind, filter, dof=value, ...
add_bc(stage, :node, (x==3, y==0), uy=0)
add_bc(stage, :face, (y==0.4), ty=:(-0.1*x)) # triangular load

run(ana)

# ❱❱❱ Post-processing

plot = DomainPlot(model,
    field = "σxx",
    colormap = :coolwarm,
    label = L"\sigma_x",
    warp = 20
)
save(plot, "2d-static.pdf")