using Serendip

# Geometry
geo = GeoModel(size=1.0)

r = add_rectangle(geo, [0, 0, 0], 1, 1, tag="base" )
c = add_disk(geo, [0.5, 0.5, 0], [0, 0, 1], 0.2, tag="base" )
holed = cut(geo, r, c)

extrude(geo, holed, [0, 0, 0.5])
cir = get_curve(geo, [0.5,0.5,0.5])
set_transfinite_curve(geo, cir, 10)
cir = get_curve(geo, [0.5,0.5,0])
set_transfinite_curve(geo, cir, 10)
# save(geo, "test.step")

mesh = Mesh(geo, quiet=true)
# save(mesh, "test.vtu")

# Finite elements
mapper = RegionMapper()
add_mapping(mapper, x>=0, MechBulk, LinearElastic, rho=10, E=1.0, nu=0.3)

select(mesh, :element, x==0, tag="base")
select(mesh, :node, x==0, tag="base")

model = FEModel(mesh, mapper, thickness=1.0)
ana   = MechAnalysis(model)

add_logger(ana, :node, (x==0, y==0, z==0.5), "a.dat")
add_logger(ana, :nodegroup, x==0, "b.dat")
add_logger(ana, :nodalreduce, x==0, "c.dat")
add_logger(ana, :ip, [0,0,0], "d.dat")
add_logger(ana, :ipgroup, x<0.1, "ea.dat")

add_monitor(ana, :node, x==0, :ux)
add_monitor(ana, :ip, [0,0,0], :σyy)

stage = add_stage(ana, nincs=1)

add_bc(stage, :node, z==0, ux=0, uy=0, uz=0)
add_bc(stage, :node, x==1, uy=0)
add_bc(stage, :face, z==0, tz=0)
add_bc(stage, :body, z>=0, wz=-1)

run(ana, tol=1e-6)

