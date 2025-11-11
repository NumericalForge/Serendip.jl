using Serendip
using Test

h  = 0.1
th = 0.05
L  = 1.0
E  = 210e6 # kPa
fy = 240e3 # kPa
H  = 0.0
nu = 0.3

geo = GeoModel()
add_block(geo, [0, 0, 0], L, 0, 0, nx=50, shape=LIN3, tag="beam")
mesh = Mesh(geo, ndim=3)

mapper = RegionMapper()
add_mapping(mapper, "beam", MechBeam, VonMises, E=E, nu=nu, fy=fy, H=H, b=th, h=h)
model = FEModel(mesh, mapper)

ana = MechAnalysis(model)
log = add_logger(ana, :node, (x==L))
mon = add_monitor(ana, :node, (x==L), :fz)

stage = add_stage(ana, nincs=30, nouts=1)
add_bc(stage, :node, (x==0), ux=0, uy=0, uz=0, rx=0, ry=0, rz=0)
add_bc(stage, :node, (x==L), uz = -0.03)
run(ana, autoinc=true)
@test log.table["fz"][end]≈-30 atol=5.0