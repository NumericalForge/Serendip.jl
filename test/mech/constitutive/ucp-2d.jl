using Serendip
using Test

h  = 0.1
th = 0.05
L  = 1.0

geo = GeoModel()
add_block(geo, [0.0, 0.0], L, h, 0, nx=50, ny=6, quadratic=true, tag="beam")
mesh = Mesh(geo)

mapper = RegionMapper()
add_mapping(mapper, "beam", MechSolid, ECP, 
    E=30e6, nu=0.25, beta=1.15, 
    fc=-30e3, epsc=-0.002, ft=3e3, wc=0.0001)

model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=th)

ana   = MechAnalysis(model)
log = add_logger(ana, :node, (x==L/2, y==h))
mon = add_monitor(ana, :node, (x==L/2, y==h), (:uy, :fy))

stage = add_stage(ana, nincs=30, nouts=20)
add_bc(stage, :node, (x==0, y==0), ux=0, uy=0)
add_bc(stage, :node, (x==L), uy=0)
add_bc(stage, :node, (x==L/2), uy=-0.004)

try run(ana, autoinc=true, tol=0.1, dTmax=0.001, rspan=0.02, quiet=false)
catch err
end

# @test log.table["fy"][end] ≈ -30 atol=0.4

tab = log.table
chart = Chart(
    xlabel = "Displacement uy [m]",
    ylabel = "Force fy [kN]",
)
add_line(chart, -tab["uy"], -tab["fy"], mark=:circle)
save(chart, "ucp-2d.pdf")

