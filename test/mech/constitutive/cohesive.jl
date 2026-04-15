using Serendip
using Test

# ❱❱❱ Geometry and mesh

geo = GeoModel()
bl1  = add_block(geo, [0, 0], 0.1, 0.1, 0.0; nx=1, ny=1, shape=:quad4, tag="bulk")
bl2  = add_block(geo, [0.1, 0], 0.1, 0.1, 0.0; nx=1, ny=1, shape=:quad4, tag="bulk")
mesh = Mesh(geo)
add_cohesive_elements(mesh, tag="interface")

left_elem = select(mesh, :element, :bulk, x<=0.1)
select(get_nodes(left_elem), tag="left")

right_elem = select(mesh, :element, :bulk, x>=0.1)
select(get_nodes(right_elem), tag="right")

E    = 27.e6
nu   = 0.2
ft   = 2.4e3
fc   = -24e3
wc   = 1.7e-4
zeta = 5.0

# ❱❱❱ Finite element analyses

trajectories = [ 
    "pure extension",
    "extension with shear",
    "pure shear",
    "compression with shear" 
]

models = (
    (cmodel=MohrCoulombCohesive, props=(E=E, nu=nu, ft=ft, mu=1.4, zeta=zeta, wc=wc), mark=:triangle),
    # (cmodel=PowerYieldCohesive, props=(E=E, nu=nu, fc=fc, ft=ft, zeta=zeta, wc=wc, alpha=1.5, gamma=0.05, theta=1.5), mark=:circle),
    (cmodel=AsinhPowerCohesive, props=(E=E, nu=nu, fc=fc, ft=ft, zeta=zeta, wc=wc, alpha=0.5, beta=1.0, theta=1.0, psi=1.4, B=ft), mark=:circle),
    (cmodel=AsinhYieldCohesive, props=(E=E, nu=nu, fc=fc, ft=ft, zeta=zeta, wc=wc, alpha=0.33, beta=0.2, theta=1.0, psi=1.4), mark=:circle),
)


grid = ChartGrid(
    # size=(16cm, 6cm), 
    size=(16cm, 24cm),
    background=:old_paper, 
    row_headers=trajectories,
    column_headers=[ "`σ_n` vs `τ`", "`w` vs `σ_n`" ],
    quiet=true,
)


for (i,trajectory) in enumerate(trajectories)
    printstyled("\nTrajectory: $(trajectory)\n\n", color=:cyan, bold=true)
    
    chart1 = Chart(legend=:top_left, xlabel="`σ_n`", ylabel="`τ`")
    chart2 = Chart(legend=:top_right, xlabel="`w`", ylabel="`σ_n`")

    for model in models
        
        cmodel = model.cmodel
        props = model.props
        mark = model.mark

        printstyled("\n$(string(cmodel)):\n\n", color=:yellow, bold=true)

        mapper = RegionMapper()
        add_mapping(mapper, "bulk", MechSolid, LinearElastic, E=E, nu=nu)
        add_mapping(mapper, "interface", MechCohesive, cmodel; props...)
        
        model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=1.0)
        ana = MechAnalysis(model)
        
        joint = select(model.elems, "interface")
        select( get_ips(joint), tag="jips" )
        log1 = add_logger(ana, :ip, "jips", string(cmodel)*".dat")
        add_monitor(ana, :ip, "jips", (:σn, :τ))
        
        stage = add_stage(ana, nincs=80, nouts=80)

        add_bc(stage, :node, "left", ux=0, uy=0)

        if trajectory == "pure extension"
            add_bc(stage, :node, "right", ux=0.0002)
        elseif trajectory == "extension with shear"
            add_bc(stage, :node, "right", ux=0.00001, uy=0.0001)
        elseif trajectory == "pure shear"
            add_bc(stage, :node, "right", ux=0.0, uy=0.001)
        elseif trajectory == "compression with shear"
            add_bc(stage, :node, "right", ux=-0.000001, uy=0.0002)
        end

        status = run(ana, autoinc=true, maxits=10, tol=0.1, rspan=0.03, dTmax=0.1, quiet=false)
        
        add_line(chart1, log1.table["σn"], log1.table["τ"], label=string(cmodel), mark=mark)
        add_line(chart2, log1.table["w"], log1.table["σn"], label=string(cmodel), mark=mark)

    end
    add_chart(grid, chart1, (i, 1))
    add_chart(grid, chart2, (i, 2))

end

save(grid, "cohesive-models.pdf")
