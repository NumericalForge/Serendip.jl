# Example of a solid beam using Von-Mises model

A finite element analysis of a cantilever beam fixed at the left end and a prescribed displacement at the right end.

```@example
using Serendip

# Data
th = 0.05  # thickness (m)
L  = 1.0   # beam length (m)
h  = 0.1   # beam height (m)
E  = 210e6 # Young modulus (kPa)
nu = 0.3   # Poisson ratio
fy = 240e3 # Yied strength (kPa)
H  = 0.0   # Hardening modulus (kPa)

# Geometry model
geo = GeoModel()
sz = 0.05 # mesh size

# List of points
p1 = add_point(geo, [0, 0, 0], size=sz)
p2 = add_point(geo, [L, 0, 0], size=sz)
p3 = add_point(geo, [L, 0, h/2], size=sz)
p4 = add_point(geo, [L, 0, h], size=sz)
p5 = add_point(geo, [0, 0, h], size=sz)

# Define a closed loop to create a surface
l1 = add_line(geo, p1, p2)
l2 = add_line(geo, p2, p3)
l3 = add_line(geo, p3, p4)
l4 = add_line(geo, p4, p5)
l5 = add_line(geo, p5, p1)
loop = add_loop(geo, [l1, l2, l3, l4, l5])
surf = add_plane_surface(geo, loop)

# Extrude the surface to create a solid
extrude(geo, surf, [0, 1, 0], heights=[th])

# Finite element mesh
mesh= Mesh(geo)

# List of element types and constitutive model
mapper = RegionModel(MechBulk, VonMises; E=E, nu=nu, fy=fy, H=H)

# A mechanical analysis context
ctx = Context()

# A finite element model
model = FEModel(mesh, mapper)

# A finite element analysis object
ana = MechAnalysis(model)

add_logger(ana, :nodalreduce, (x==L, z==h/2), "file.dat")

# List of monitors
add_monitor(ana, :node, (x==L, y==0, z==h/2), :fz)

# List of boundary conditions
add_bc(ana.current_stage, :node, x==0, ux=0, uy=0, uz=0)
add_bc(ana.current_stage, :node, (x==L, z==h/2), uz=-0.08)

# Adds a load stage
addstage!(ana, bcs, nincs=20, nouts=1)

# Run the analysis
solve!(ana, autoinc=true)
```
