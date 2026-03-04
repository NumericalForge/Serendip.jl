using Serendip
using Test

# Shared pullout setup and approximately aligned bond-slip parameters
τmax   = 12.0
τres   = 3.0
s_peak = 0.001
s_res  = 0.004
α      = 0.5
β      = 0.4
ks     = (τmax / s_peak) * 10.0
kn     = 5.0e3
p      = 0.25
outer_ip_coord = [0.5, 5.95, 0.5] # near loaded end (outermost bond-slip IP)

chart = Chart(
    xlabel="slip s",
    ylabel="bond stress τ",
    ylimits=[0.0, 1.4*τmax],
    xlimits=[0.0, 1.5*s_res],
    legend=:top_right,
)

cases = [
    (
        name = "LinearBondSlip",
        cmodel = LinearBondSlip,
        kwargs = (; ks=ks, kn=kn),
    ),
    (
        name = "CebBondSlip",
        cmodel = CebBondSlip,
        kwargs = (; taumax=τmax, taures=τres, s1=s_peak, s2=1.1*s_peak, s3=s_res, alpha=α, ks=ks, kn=kn),
    ),
    (
        name = "PowerExpBondSlip",
        cmodel = PowerExpBondSlip,
        kwargs = (; taumax=τmax, taures=τres, speak=s_peak, sc=s_res, alpha=α, beta=10.0, ks=ks, kn=kn),
    ),
    (
        name = "CyclicBondSlip",
        cmodel = CyclicBondSlip,
        kwargs = (; taumax=τmax, taures=τres, speak=s_peak, sres=s_res, alpha=α, beta=β, ks=ks, kn=kn, p=p),
    ),
]

for case in cases
    printstyled("\n", case.name, "\n", color=:yellow, bold=true)
    geo = GeoModel(quiet=true)
    add_block(geo, [0.0, 0.0, 0.0], 1.0, 6.0, 1.0, nx=1, ny=5, nz=1, tag="solids")
    p1 = add_point(geo, [0.5, 3.0, 0.5])
    p2 = add_point(geo, [0.5, 6.0, 0.5])
    edge = add_line(geo, p1, p2)
    add_path(geo, [edge]; tag="bars", interface_tag="joints")

    mesh = Mesh(geo)

    bar_points = get_nodes(select(mesh, :element, "bars"))
    select(bar_points, y==6, tag="tip")
    select(get_nodes(select(mesh, :element, "solids")), tag="fixed")

    mapper = RegionMapper()
    add_mapping(mapper, "solids", MechBulk, LinearElastic, E=24e6, nu=0.2)
    add_mapping(mapper, "bars", MechBar, LinearElastic, E=200e6, A=0.00011)
    add_mapping(mapper, "joints", MechBondSlip, case.cmodel; case.kwargs..., p=p)

    model = FEModel(mesh, mapper)
    ana = MechAnalysis(model)

    joint_log = add_logger(ana, :ip, outer_ip_coord)

    stage = add_stage(ana)
    add_bc(stage, :node, "fixed", ux=0, uy=0, uz=0)
    add_bc(stage, :node, "tip", uy=+0.007)

    status = run(ana, autoinc=true, tol=0.01, dTmax=0.05)

    table = joint_log.table

    add_series(chart, :line, table["s"], table["τ"], label=case.name, mark=:circle)

    case.name == "LinearBondSlip" && continue # skip assertions for linear case

    @test maximum(abs.(table["τ"])) ≈ τmax atol=0.01*τmax
    @test abs(table["τ"][end]) ≈ τres atol=0.01*τres
end

save(chart, "bondslip-monotonic-comparison.pdf")

nothing