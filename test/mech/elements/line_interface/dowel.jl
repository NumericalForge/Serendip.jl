using Serendip

geo = GeoModel(size=0.1)

# Bulk
add_rectangle(geo, [0.0, 0.0, 0.0], 1.0, 1.0, tag="bulks")

# Bar
p1 = add_point(geo, [0.2, 0.2, 0.0])
p2 = add_point(geo, [0.5, 0.5, 0.0])
p3 = add_point(geo, [0.6, 0.5, 0.0])
p4 = add_point(geo, [0.9, 0.5, 0.0])
p5 = add_point(geo, [1.2, 0.5, 0.0])

c1 = add_bezier(geo, [p1, p2, p3, p4])
l1 = add_line(geo, p4, p5)
add_path(geo, [c1, l1], tag="lines", interface_tag="linejoints")

mesh = Mesh(geo)
save(mesh, "dowel.vtu")

# FEM
mapper = RegionMapper()
add_mapping(mapper, :bulk, MechBulk, LinearElastic, E=24e2, nu=0.2)
add_mapping(mapper, "lines", MechBeam, LinearElastic, E=200e6, A=0.00011)
add_mapping(mapper, "linejoints", MechBondSlip, LinearBondSlip, kn=5000, ks=6000, p=0.25)

model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=0.1)

ana = MechAnalysis(model)

stage = add_stage(ana, nincs=10, nouts=2)
add_bc(stage, :node, y==0, ux=0, uy=0)
add_bc(stage, :node, (x==1.2, y==0.5), ux=0.01, uy=-0.01)

run(ana, tol=0.01, maxits=3, autoinc=true)

save(model, "dowel.vtu")
