#❱❱❱ Simple beam example ❰❰❰#

using Serendip

# ❱❱❱ Mesh generation

l = 4.0   # x direction
b = 0.3   # y direction
h = 0.6   # z direction

geo  = GeoModel(size=0.1)
box  = add_box(geo, [0.0, 0.0, 0.0], l, b, h)
mesh = Mesh(geo)
save(mesh, "simple-beam-mesh.vtu")

plot = DomainPlot(mesh)
save(plot, "simple-beam-mesh.pdf")

# ❱❱❱ Finite element analysis

mapper = RegionModel(MechBulk, LinearElastic, E=22.0e6, nu=0.25)
model = FEModel(mesh, mapper)

ana   = MechAnalysis(model, outkey="simple-beam")
stage = add_stage(ana, nincs=10, nouts=5)

add_bc(stage, :node, (x==0, z==0), ux=0, uy=0, uz=0) # fixed support
add_bc(stage, :node, (x==l, z==0), uy=0, uz=0) # roller support
add_bc(stage, :face, (z==h), tz=-10) # surface load

run(ana)

# ❱❱❱ Post-processing

plot = DomainPlot(model, field="uy", colormap=:spectral)
save(plot, "simple-beam.pdf")

