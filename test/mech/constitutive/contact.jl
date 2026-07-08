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
    (cmodel=MohrCoulombContact, props=(ks=ks, kn=kn, ft=ft, wc=wc, mu=1.4)),
    (cmodel=CoulombContact, props=(kn=kn, ks=ks, mu=mu)),
)

for trajectory in trajectories
    for model in models
        @announced_testset "$(trajectory) / $(model.cmodel)" begin
            cmodel = model.cmodel
            props = model.props

            mapper = RegionMapper()
            add_mapping(mapper, :solid, MechSolid, LinearElastic, E=E, nu=nu)
            add_mapping(mapper, "interface", MechContact, cmodel; props...)

            fe_model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=1.0)
            ana = MechAnalysis(fe_model)
            select(fe_model, :element, "interface", :ip, tag="c-ips")
            log1 = add_logger(ana, :ip, "c-ips", string(cmodel) * ".dat")
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

            status = run(ana, autoinc=true, tol=0.1, rspan=0.03, dTmax=0.1, quiet=true)
            @test status.successful
            @test size(log1.table, 1) > 0
        end
    end
end
