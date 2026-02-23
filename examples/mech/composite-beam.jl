#❱❱❱ Composite concrete-steel beam ❰❰❰#

using Serendip

# ❱❱❱ Geometry and mesh generation
ℓ  = 4.5   # length
hf = 0.06  # flange thickness
hw = 0.12  # web thickness
h  = hf+hw # total height
b  = 0.35  # flange width

geo = GeoModel()
sz = 0.2

# profile (Defined in Y-Z plane now, X=0)
p1 = add_point(geo, [0, 0.143, 0], size=sz)
p2 = add_point(geo, [0, 0.207, 0], size=sz)
p3 = add_point(geo, [0, 0.207, 0.0063], size=sz)
p4 = add_point(geo, [0, 0.177, 0.0063], size=sz)
p5 = add_point(geo, [0, 0.177, 0.1137], size=sz)
p6 = add_point(geo, [0, 0.207, 0.1137], size=sz)
p7 = add_point(geo, [0, 0.207, 0.120], size=sz)
p8 = add_point(geo, [0, 0.143, 0.120], size=sz)
p9 = add_point(geo, [0, 0.143, 0.1137], size=sz)
p10 = add_point(geo, [0, 0.173, 0.1137], size=sz)
p11 = add_point(geo, [0, 0.173, 0.0063], size=sz)
p12 = add_point(geo, [0, 0.143, 0.0063], size=sz)

l1 = add_line(geo, p1, p2)
l2 = add_line(geo, p2, p3)
l3 = add_line(geo, p3, p4)
l4 = add_line(geo, p4, p5)
l5 = add_line(geo, p5, p6)
l6 = add_line(geo, p6, p7)
l7 = add_line(geo, p7, p8)
l8 = add_line(geo, p8, p9)
l9 = add_line(geo, p9, p10)
l10 = add_line(geo, p10, p11)
l11 = add_line(geo, p11, p12)
l12 = add_line(geo, p12, p1)

loop = add_loop(geo, [l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12])
section_profile = add_plane_surface(geo, loop)

# Extrude along X (Length)
profile = extrude(geo, section_profile, [ℓ, 0, 0])

# Flange (Concrete slab): Dimensions are now [Length, Width, Height]
flange = add_box(geo, [0, 0, hw], ℓ, b, hf) 
fragment(geo, profile, flange)

# reinforcement paths (Longitudinal along X)
# Coordinates changed from [width_pos, length_pos, z] to [length_pos, width_pos, z]
pa = add_point(geo, [0, 0.05, 0.170], size=sz)      # Start
pb = add_point(geo, [2.25, 0.05, 0.140], size=sz)   # Mid (control point)
pc = add_point(geo, [ℓ, 0.05, 0.170], size=sz)      # End

bez  = add_bezier(geo, [pa, pb, pc])
path = add_path(geo, [bez], tag="rebar", interface_tag="rebar-interface")

# Array copies along Y (Width) instead of X
add_array(geo, path; ny=5, dy=0.05)

# stud paths
# Position changed: [Length_pos, Width_pos, Z]
p1 = add_point(geo, [0.3, b/2, hw - 0.02 ])
p2 = add_point(geo, [0.3, b/2, hw + 0.04 ])
stud = add_line(geo, p1, p2)
path = add_path(geo, [stud], tag="stud", tips=:both, interface_tag="stud-interface", tip_tag="stud-tip")

# Array copies along X (Length) instead of Y
add_array(geo, path, nx=14, dx=0.3)

# embedded point
add_point(geo, [ℓ/2, b/2, 0.0], size=sz*0.5, embedded=true)

# setting mesh size
select(geo, :surface, [0.1,0.1,0.18], tag="surf")
set_size(geo, :surface, "surf", sz*0.8)

save(geo, "test.step")

mesh = Mesh(geo, quadratic=true)
save(mesh, "medium-mesh.vtu")

# applying tags to regions
select(mesh, :element, :bulk, z<=0.120, tag="steel-beam")
select(mesh, :element, :bulk, z>=0.120, tag="concrete")
# select(mesh, :element, :line, tag="rebar")
# select(mesh, :element, :line_interface, tag="interface")

add_cohesive_elements(mesh, z>=0.120, tag="cohesive")
add_contact_elements(mesh, "steel-beam", "concrete", tag="contact")

# ❱❱❱ Finite element modeling

mapper = RegionMapper()
add_mapping(mapper, "steel-beam", MechBulk, LinearElastic, E=200e6, nu=0.2)
add_mapping(mapper, "concrete", MechBulk, LinearElastic, E=22.1e6, nu=0.25)
add_mapping(mapper, "cohesive", MechCohesive, MohrCoulombCohesive, E=22.1e6, ft=3.0e3, mu=1.4, GF=0.1)
add_mapping(mapper, "contact", MechContact, LinearInterface, ks=1e7, kn=1e10)
add_mapping(mapper, "rebar", MechBar, LinearElastic, E=200e6, A=0.005)
add_mapping(mapper, "stud", MechBeam, LinearElastic, E=200e6, A=0.005)
add_mapping(mapper, "rebar-interface", MechBondSlip, LinearBondSlip, ks=1e6, kn=1e10, p=0.01) 
add_mapping(mapper, "stud-interface", MechBondSlip, LinearBondSlip, ks=1e6, kn=1e10, p=0.01) 
add_mapping(mapper, "stud-tip", MechBondTip, LinearTip, k=1e10) 

model = FEModel(mesh, mapper)

ana = MechAnalysis(model, outkey="composite", outdir="composite")
stage = add_stage(ana, nouts=50)

# Boundary Conditions adjusted for X-axis length
# Supports at x=0 and x=ℓ
add_bc(stage, :node, (x==0, z==0), ux=0, uy=0, uz=0) # fixed support
add_bc(stage, :node, (x==ℓ, z==0), uy=0, uz=0)       # roller support
add_bc(stage, :face, (z==0.180), tz=-20)             # surface load (Z normal remains same)

run(ana, autoinc=true, tol=5.0, tangent_scheme=:ralston, rspan=0.03)
save(model, "medium-model.vtu")

# ❱❱❱ Post-processing

plot = DomainPlot(model,
    field = "σx´", # stress in the bars (local axis likely aligns with bar, so this might stay)
    field_mult = 1e-3,
    colormap = :rainbow,
    edge_color = :gray,
    edge_width = 0.1,
    view_mode = :wireframe,
    colorbar = :bottom,
    label = t"σ_(x')",
    warp = 20,
)
save(plot, "composite-beam.pdf")