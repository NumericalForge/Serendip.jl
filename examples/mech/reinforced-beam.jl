#❱❱❱ Reinforced beam example ❰❰❰#

using Serendip

# ❱❱❱ Mesh generation

# geo = GeoModel(size=0.08)
geo = GeoModel(size=0.15)

ℓ = 6.0   # x direction
b = 0.3   # y direction
h = 0.6   # z direction
box = add_box(geo, [0.0, 0.0, 0.0], ℓ, b, h)

p1 = add_point(geo, [ℓ/2, 0, h])
p2 = add_point(geo, [ℓ/2, b, h])
l1 = add_line(geo, p1, p2)
fragment(geo, l1, box)

c     = 0.05 # cover
pa    = add_point(geo, [c, c, c])
pb    = add_point(geo, [ℓ-c, c, c])
lin   = add_line(geo, pa, pb)
path1 = add_path(geo, [lin], tag="bottom-steel")
# add_array(geo, path1; nx=1, ny=5, nz=2, dx=0.0, dy=0.05, dz=0.05)
add_array(geo, path1; ny=3, dy=0.1)

p3 = add_point(geo, [c, c, h-c])
p4 = add_point(geo, [ℓ-c, c, h-c])
l2 = add_line(geo, p3, p4)
path2 = add_path(geo, [l2], tag="bars", interface_tag="bond")
add_array(geo, path2, ny=3, dy=0.1)

save(geo, "medium-geo.step")
save(geo, "medium-geo.geo_unrolled")

# Mesh generation
mesh = Mesh(geo)
save(mesh, "reinforced-beam-mesh.vtu")
mplot = DomainPlot(mesh, view_mode = :wireframe, elevation=10, azimuth=30)
save(mplot, "reinforced-beam-mesh.pdf")

select(mesh, :element, :bulk, tag="concrete")
select(mesh, :element, :line, tag="bars")
select(mesh, :element, :line_interface, tag="bond")

# ❱❱❱ Finite element modeling

Ec = 24e6 # kPa
fc = -30.00e3
ft = 3.0e3
Es = 210e6 # kPa
fy = 240e3 # kPa
φ  = 0.01 # m
p  = π*φ
A  = π*φ^2/4
H  = 0.0
nu = 0.3

mapper = RegionMapper()
add_mapping(mapper, "concrete", MechBulk, LinearElastic, E=E, nu=nu)
# add_mapping(mapper, "concrete", MechBulk, UCP, E=Ec, nu=0.25, alpha=0.66, beta=1.15, fc=fc, epsc=-0.002, eta=2.2, ft=ft, wc=0.0001) 
# UCP is not fully working in tension yet
add_mapping(mapper, "bars", MechBar, VonMises, E=Es, A=A, fy=fy, H=H)
add_mapping(mapper, "bond", MechBondSlip, LinearBondSlip, ks=1e7, kn=1e9, p=p)

model = FEModel(mesh, mapper)

# Mechanical analysis
ana   = MechAnalysis(model, outkey="reinforced")
stage = add_stage(ana, nouts=50)

add_bc(stage, :node, (x==0, z==0), ux=0, uy=0, uz=0) # fixed support
add_bc(stage, :node, (x==ℓ, z==0), uy=0, uz=0) # roller support
# add_bc(stage, :face, (z==h), tz=-100) # surface load
add_bc(stage, :edge, (x==ℓ/2, z==h), uz=-0.01)

run(ana, tol=1.0, rspan=0.05, autoinc=true)

# ❱❱❱ Post-processing

mplot = DomainPlot(model, field="σx´", colormap=:spectral, view_mode=:wireframe)
save(mplot, "reinforced-beam.pdf")
