using Serendip
using Test

L  = 6.0
R  = 3.0
nx = 20
n  = 20

# Generate the geometry model
geo = GeoModel()
p1  = add_point(geo, [-R, 0, 0])
p2  = add_point(geo, [-R, L, 0])
l   = add_line(geo, p1, p2)
set_transfinite_curve(geo, l, nx)

# Generate the mesh
mesh = Mesh(geo, quadratic=true)
mesh = revolve(mesh, angle=180, n=n, base=[0,0,0], axis=[0,1,0])

# Finite element model
mapper = RegionModel(MechShell, LinearElastic, E=3e4, nu=0.3, thickness=0.03)

model = FEModel(mesh, mapper, ndim=3)
ana   = MechAnalysis(model)

# Add logger and monitor
add_logger(ana, :nodegroup, (y==L/2, z==R))
add_monitor(ana, :node, (y==L/2, z==R), :uz)

# Add a stage and boundary conditions
stage = add_stage(ana)
add_bc(stage, :node, y==0, ux=0, uy=0, rz=0)
add_bc(stage, :node, y==L, ux=0, uz=0, rz=0)
add_bc(stage, :node, z==0, uz=0, ry=0)
add_bc(stage, :node, (y==L/2, z==R), fz=-0.01)
add_bc(stage, :face, z>=0, tz=-0.001)

# Run the analysis
run(ana)

@test ana.data.monitors[1].table.uz[end] ≈ -0.001823 atol=1e-5