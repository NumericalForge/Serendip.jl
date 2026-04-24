using Serendip
using Test

# ❱❱❱ Geometry and mesh

geo = GeoModel()
bl1  = add_block(geo, [0, 0], 0.1, 0.1, 0.0; nx=1, ny=1, shape=:quad4, tag="left")
bl2  = add_block(geo, [0.1, 0], 0.1, 0.1, 0.0; nx=1, ny=1, shape=:quad4, tag="right")
mesh = Mesh(geo)

add_contact_elements(mesh, tag="interface")

select(mesh, :element, "left", :node, tag="left")
select(mesh, :element, "right", :node, tag="right")

# left_elem = select(mesh, :element, :bulk, x<=0.1)
# select(get_nodes(left_elem), tag="left")

# right_elem = select(mesh, :element, :bulk, x>=0.1)
# select(get_nodes(right_elem), tag="right")

E    = 27.e6
nu   = 0.2
ft   = 2.4e3
wc   = 1.7e-4
kn   = 2e10
ks   = 1.8e7
mu   = 1.5
zeta = 5.0

# ❱❱❱ Finite element analyses

trajectories = [ 
    "pure extension",
    "extension with shear",
    "pure shear",
    "compression with shear" 
]

models = (
    (cmodel=MohrCoulombContact, props=(ks=ks, kn=kn, ft=ft, wc=wc, mu=1.4), mark=:triangle),
    (cmodel=CoulombContact, props=(kn=kn, ks=ks, mu=mu), mark=:circle),
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
    printstyled("\n\nTrajectory: $(trajectory)\n\n", color=:yellow, bold=true)
    
    chart1 = Chart(legend=:top_left, xlabel="`σ_n`", ylabel="`τ`")
    chart2 = Chart(legend=:top_right, xlabel="`w`", ylabel="`σ_n`")

    for model in models
        cmodel = model.cmodel
        props  = model.props
        mark   = model.mark

        printstyled("\n$(string(cmodel))\n\n", color=:yellow, bold=true)

        mapper = RegionMapper()
        add_mapping(mapper, :solid, MechSolid, LinearElastic, E=E, nu=nu)
        add_mapping(mapper, "interface", MechContact, cmodel; props...)
        
        model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=1.0)
        ana = MechAnalysis(model)
        # global mm = model
        
        select(model, :element, "interface", :ip, tag="c-ips" )
        log1 = add_logger(ana, :ip, "c-ips", string(cmodel)*".dat")
        add_monitor(ana, :ip, "c-ips", (:σn, :τ))
        
        stage = add_stage(ana, nincs=80, nouts=20)

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

        status = run(ana, autoinc=true, tol=0.1, rspan=0.03, dTmax=0.1, quiet=false)
        # @show log1.table
        
        add_line(chart1, log1.table["σn"], log1.table["τ"], label=string(cmodel), mark=mark)
        add_line(chart2, log1.table["w"], log1.table["σn"], label=string(cmodel), mark=mark)

    end
    add_chart(grid, chart1, (i, 1))
    add_chart(grid, chart2, (i, 2))

end

save(grid, "contact-models.pdf")
