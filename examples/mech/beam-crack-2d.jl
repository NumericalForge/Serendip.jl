# ❱❱❱ Cracking in a 2D cantilever beam ❰❰❰#

using Serendip

ℓ = 0.4   # x direction
b = 0.025 # y direction
h = 0.1   # z direction

geo = GeoModel(size=0.01)
add_rectangle(geo, [0,0,0], ℓ, h; tag="solid")

p1 = add_point(geo, [0.00*ℓ, 0.95*h, 0])
p2 = add_point(geo, [0.95*ℓ, 0.95*h, 0])
bar = add_line(geo, p1, p2)
add_path(geo, [bar], tag="bar", interface_tag="bar-interface")

mesh = Mesh(geo)

select(mesh.elems, :solid, tag="solid")
add_cohesive_elements(mesh, x<ℓ/2, tag="cohesive", implicit=false)

# finite element analysis
E     = 37.e6
ft    = 2.4e3
fc    = -24e3
alpha = 1.5
GF    = 0.07
zeta  = 5

mapper = RegionMapper()
add_mapping(mapper, "solid", MechSolid, LinearElastic, E=E, nu=0.2)
add_mapping(mapper, "cohesive", MechCohesive, MohrCoulombCohesive, E=E, nu=0.2, ft=ft, GF=GF, mu=1.5, psi=2.5, zeta=zeta)
add_mapping(mapper, "bar", MechBar, LinearElastic, E=2e5, A=0.005)
add_mapping(mapper, "bar-interface", MechBondSlip, LinearBondSlip, ks=1e10, kn=1e9, p=0.01*3)

model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=b)
ana   = MechAnalysis(model, outkey="crack-2d", outdir="crack-2d")

log1 = add_logger(ana, :nodalreduce, x==ℓ)
add_monitor(ana, :node, (x==ℓ, y==0), :uy)

stage = add_stage(ana, nincs=10, nouts=50)
add_bc(stage, :node, x==0, ux=0, uy=0)
add_bc(stage, :node, x==ℓ, uy=-0.002)

run(ana, autoinc=true, maxits=3, tol=0.1, rspan=0.03, quiet=false)

save(model, "crack-2d.vtk")

# ❱❱❱ Post-processing

plot_kwargs = (
    field      = "σxx",
    field_mult = 1e-3,
    warp       = 25,
    colormap   = :spectral,
    diverging  = true,
    line_color = :gray,
    line_width = 0.1,
    colorbar   = :bottom,
    label      = "`σ_(x x)` [MPa]",
)

# Charts

plot = DomainPlot()
add_plot(plot, model; plot_kwargs...)
save(plot, "beam-crack-2d.pdf")

chart = Chart(
    xlabel = "`u_y` [mm]",
    ylabel = "`f_y` [kN]",
)
add_line(chart, -log1.table[:uy]*1e3, -log1.table[:fy], mark=:circle)
save(chart, "beam-crack-2d-chart.pdf")

# Crack video
nouts  = stage.nouts

# nouts  = 50
# outdir = "crack-2d"
# outkey = "crack-2d"

video = VideoBuilder(freeze_scale=true, bounds_factor=1.05)

for i in 0:nouts
    file = joinpath(ana.data.outdir, "$(ana.data.outkey)-$i.vtu")
    # file = joinpath(outdir, "$outkey-$i.vtu")
    # print("frame $i\r")
    frame = DomainPlot()
    add_plot(frame, Mesh(file; quiet=true); plot_kwargs...)
    add_frame(video, frame)
end
print("\e[K")

save(video, "beam-crack-2d.mp4")
