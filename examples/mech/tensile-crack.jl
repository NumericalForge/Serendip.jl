using Serendip

# ❱❱❱ Geometry and mesh

hplate    = 0.005
hsample   = 0.12
htop      = hsample + hplate
thickness = 0.01

sz = 0.004

geo = GeoModel(size=sz)
p1  = add_point(geo, [0.0, 0.0, 0.0])
p2  = add_point(geo, [0.06, 0.0, 0.0])
p3  = add_point(geo, [0.06, 0.0665, 0.0])
p4  = add_point(geo, [0.05, 0.0665, 0.0])
p5  = add_point(geo, [0.05, 0.0685, 0.0])
p6  = add_point(geo, [0.06, 0.0685, 0.0])
p7  = add_point(geo, [0.06, 0.12, 0.0])
p8  = add_point(geo, [0.0, 0.12, 0.0])
p9  = add_point(geo, [0.0, 0.0535, 0.0])
p10 = add_point(geo, [0.01, 0.0535, 0.0])
p11 = add_point(geo, [0.01, 0.0515, 0.0])
p12 = add_point(geo, [0.0, 0.0515, 0.0])

conc  = add_polygon(geo, [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12])
steel = add_rectangle(geo, [0.0, 0.12, 0.0], 0.06, 0.005)
fragment(geo, conc, steel)

q1a = add_point(geo, [0.005, 0.0375, 0.0], embedded=true)
q2b = add_point(geo, [0.005, 0.0675, 0.0], embedded=true)
q2a = add_point(geo, [0.055, 0.0525, 0.0], embedded=true)
q1b = add_point(geo, [0.055, 0.0825, 0.0], embedded=true)

set_refinement(geo, [0.03, 0.06, 0], 0.03, 0.02, 1, 0.4*sz, sz, gradient=0.7, roundness=0.5)

mesh = Mesh(geo)

select(mesh, :element, y<=0.12, tag="concrete")
select(mesh, :element, y>=0.12, tag="steel")

# ❱❱❱ Finite element analysis

# Material properties
  
Ec = 31000e3 # modulo de elasticidade do concreto (kPa)
ft = 3e3	 # resistencia a tração concreto (kPa)
GF = 0.05	 # energia de fratura do concreto (kN/m)
mu = 1.4
Es = 200e6  # modulo de elasticidade do aço (kPa)

w = 0.0002
w = 0.0001

select(mesh, :element, tag="concrete")
add_cohesive_elements(mesh, (y>=0.04, y<=0.08), tag="cohesive")

# Finite element model
mapper = RegionMapper()
add_mapping(mapper, "concrete", MechCont, LinearElastic, E=Ec, nu=0.2)
add_mapping(mapper, "steel", MechCont, LinearElastic, E=Es, nu=0.3)
add_mapping(mapper, "cohesive", MechCohesive, MohrCoulombCohesive, E=Ec, nu=0.2, ft=ft, mu=mu, GF=GF, ft_law=:hordijk, zeta=10)

model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=thickness)
ana   = MechAnalysis(model, outkey="tensile-crack", outdir="tensile-crack")

# Loggers
topL = add_logger(ana, :nodalreduce, (y==htop), "top.table")
q1aL = add_logger(ana, :node, q1a.coord, "q1a.dat")
q1bL = add_logger(ana, :node, q1b.coord, "q1b.dat")
q2aL = add_logger(ana, :node, q2a.coord, "q2a.dat")
q2bL = add_logger(ana, :node, q2b.coord, "q2b.dat")

# Boundary conditions and solution
stage = add_stage(ana, nouts=50)
add_bc(stage, :node, (y==0), ux=0, uy=0)
add_bc(stage, :node, (y==htop), uy=w)

run(ana, autoinc=true, tol=0.01, rspan=0.03, dTmax=0.005, tangent_scheme=:ralston, quiet=false)

# ❱❱❱ Post-processing

chart = Chart(xlabel="CMOD", ylabel="Load [kN]")

δ = 0.5*(q1bL.table[:uy] - q1aL.table[:uy] + q2bL.table[:uy] - q2aL.table[:uy])
F = topL.table[:fy]

add_series(chart, δ*1e3, F, mark=:circle)
save(chart, "tensile-crack.pdf")
