# # Static 2D analysis
#
# This page is generated from `examples/docs/static-2d.jl` using `Literate.jl`.
#
# ## Loading Serendip

using Serendip

# ## Geometry and mesh generation

geo = GeoModel(size=0.1)
add_rectangle(geo, [0.0, 0.0], 3.0, 0.4) # geo, corner, width, height

mesh = Mesh(geo)
select(mesh, :element, tag="solids"); # tag all elements as "solids"

# ## Finite element model

mapper = RegionMapper()
add_mapping(mapper, "solids", MechBulk, LinearElastic, E=200e6, nu=0.2)
model = FEModel(mesh, mapper, stress_state=:plane_stress)

ana = MechAnalysis(model, outkey="static-2d")

add_logger(ana, :node, (x==0, y==0), "one-node.dat")
add_logger(ana, :ip, (y<0.025), "ip-list.dat")

add_monitor(ana, :node, (x==3, y==0.4), :uy)

stage = add_stage(ana)
add_bc(stage, :node, (x==0, y==0), ux=0, uy=0)
add_bc(stage, :node, (x==3, y==0), uy=0)
add_bc(stage, :face, (y==0.4), ty=-10.0*x)

run(ana)

# ## Post-processing

plot = DomainPlot(model,
    field = "σxx",
    colormap = :coolwarm,
    label = t"σ_x",
    warp = 10000
)
save(plot, "static-2d.svg")

# ## Generated figure
#
# ```@raw html
# <img src="static-2d.svg" alt="static-2d plot" width="100%">
# <p><a href="static-2d.svg">Open/download SVG</a></p>
# ```
