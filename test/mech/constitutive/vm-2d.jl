using Serendip
using Test

h  = 0.1
th = 0.05
L  = 1.0
E  = 210e6 # kPa
fy = 240e3 # kPa
H  = 0
nu = 0.3

geo = GeoModel()
add_block(geo, [0.0, 0.0], [L, h], nx=50, ny=2, shape=QUAD8, tag="beam")
mesh = Mesh(geo)

mapper = RegionMapper()
add_mapping(mapper, "beam", MechBulk, VonMises, E=E, nu=nu, fy=fy, H=H)

model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=th)

ana   = MechAnalysis(model)
log = add_logger(ana, :node, (y==h/2, x==1))
mon = add_monitor(ana, :node, (y==h/2, x==1), :fy)

stage = add_stage(ana, nincs=30, nouts=1)
add_bc(stage, :node, (x==0), ux=0)
add_bc(stage, :node, (x==0, y==h/2), uy=0)
add_bc(stage, :node, (x==1, y==h/2), uy = -0.03)

run(ana, autoinc=true)

@test log.table["fy"][end] ≈ -30 atol=0.4

makeplots = false
if @isdefined(makeplots) && makeplots
    tab = log.table
    chart = Chart(;
        xlabel = "Displacement uy [m]",
        ylabel = "Force fy [kN]",
    )
    add_series(chart, -tab["uy"], -tab["fy"], mark=:circle)
    save(chart, "vm-2d.pdf")
end

