using Serendip
using Test

#❱❱❱ Mesh
geo = GeoModel()
add_block(geo, [0,0,0], 2, 0, 0, tag="beam", nx=2, shape=:lin3)

mesh = Mesh(geo, ndim=2)
#❱❱❱ FEM analysis
mapper = RegionMapper()
add_mapping(mapper, "beam", MechBeam, LinearElastic, E=300e3, b=0.05, h=0.2)

# FE model
model = FEModel(mesh, mapper)
select(model, :edge, x>=1, tag="right")

ana = MechAnalysis(model, outdir="beam")

log1 = add_logger(ana, :node, x==1)
add_monitor(ana, :node, x==1, :uy)

stage = add_stage(ana)

add_bc(stage, :node, (x==0, y==0), ux=0, uy=0)
add_bc(stage, :node, (x==2, y==0), ux=0, uy=0)
add_bc(stage, :edge, "right", qy=-12.0)
add_bc(stage, :node, (x==1, y==0), fy=-1)
add_bc(stage, :node, (x==1, y==0), mz=1)

run(ana) 

@test log1.table[:uy][end] ≈ -0.14447 atol=1e-5

# ❱❱❱ Post-processing
plot = DomainPlot(model,
    size            = (280,100),
    warp            = 2,
    field           = "σx´",
    colormap        = :coolwarm,
    line_elem_width = 3,
    colorbar_ratio  = 1,
    font_size       = 8,
    outerpad        = 15.0,
    node_labels     = true
)

save(plot, "beam.pdf")
