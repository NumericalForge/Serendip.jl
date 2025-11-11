using Serendip
using Test


th = 0.05
E  = 210e6 # kPa
fy = 240e3 # kPa
H  = 0
nu = 0.3

geo = GeoModel()
add_block(geo, [0, 0, -th], th, 1.0, 2*th, nx=1, ny=15, nz=2, shape=HEX20, tag="beam")
mesh = Mesh(geo)


mapper = RegionMapper()
add_mapping(mapper, "beam", MechBulk, VonMises, E=E, nu=nu, fy=fy, H=H)
model = FEModel(mesh, mapper)

ana = MechAnalysis(model)
log = add_logger(ana, :node, (x==th/2, y==1, z==0))
add_monitor(ana, :node, (x==th/2, y==1, z==0), :fz)

stage = add_stage(ana, nincs=20, nouts=1)
add_bc(stage, :node, (y==0), uy=0)
add_bc(stage, :node, (y==0, z==0), uz=0)
add_bc(stage, :node, (x==th/2, y==0, z==0), ux=0)
add_bc(stage, :node, (x==th/2, y==1, z==0), uz=-0.03)

run(ana, autoinc=true, tol=0.001)
@test log.table.fz[end]≈-30 atol=1.1

makeplots = false
if @isdefined(makeplots) && makeplots
    tab = log.table
    chart = Chart(;
        xlabel = "Displacement uz [m]",
        ylabel = "Force fz [kN]",
    )
    add_series(chart, -tab["uz"], -tab["fz"], mark=:circle)
    save(chart, "vm-3d.pdf")
end
