using Serendip

# Geometry
geo = GeoModel()

r = add_rectangle(geo, [0, 0, 0], 1, 1, tag="base" )
c = add_disk(geo, [0.5, 0.5, 0], [0, 0, 1], 0.2, tag="base" )
holed = cut(geo, r, c)

extrude(geo, holed, [0, 0, 0.5])
cir = get_curve(geo, [0.5,0.5,0.5])
set_transfinite_curve(geo, cir, 20)
cir = get_curve(geo, [0.5,0.5,0])
set_transfinite_curve(geo, cir, 20)
save(geo, "test.step")

mesh = Mesh(geo, quiet=true)
save(mesh, "test.vtu")

# Finite elements
mapper = RegionMapper()
add_mapping(mapper, x>=0, MechBulk, LinearElastic, rho=10, E=1.0, nu=0.3, alpha_s=5/6)

select(mesh, :element, x==0, tag="base")
select(mesh, :node, x==0, tag="base")

model = FEModel(mesh, mapper, thickness=1.0)
ana      = MechAnalysis(model)

add_logger(ana, :node, x==0, "a.table")
add_logger(ana, :nodegroup, x==0, "b.table")
add_logger(ana, :nodalreduce, x==0, "c.table")
add_logger(ana, :ip, x==0, "d.table")
add_logger(ana, :ipgroup, x==0, "ea.table")

add_monitor(ana, :node, x==0, :ux)
add_monitor(ana, :ip, x==0, :ux)

stage = add_stage(ana, nincs=10)

add_bc(stage, :node, x==0, ux=0, uy=0)
add_bc(stage, :node, x==1, uy=0)
add_bc(stage, :face, z==0, tz=0)
add_bc(stage, :body, z==1, fz=1)


run(ana, tol=1e-6)

