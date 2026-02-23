using Serendip
using Test

# ❱❱❱ Geometry and mesh

geo = GeoModel()
bl1  = add_block(geo, [0, 0], 0.1, 0.1, 0.0; nx=1, ny=1, shape=QUAD4, tag="bulk")
bl2  = add_block(geo, [0.1, 0], 0.1, 0.1, 0.0; nx=1, ny=1, shape=QUAD4, tag="bulk")
mesh = Mesh(geo)
select(mesh, :element, x<=0.1, tag="left")

add_contact_elements(mesh, tag="interface")

E    = 27.e6
nu   = 0.2
ft   = 2.4e3
fc   = -24e3
zeta = 5.0
wc   = 1.7e-4
kn   = 2e10
ks   = 1.8e7

# ❱❱❱ Finite element analyses

models_props_d = Dict(
    MohrCoulombContact => (ks=ks, kn=kn, ft=ft, wc=wc, mu=0.9)
)

model_elem_d = Dict(
    MohrCoulombContact => MechContact
)

# ❱❱❱ Finite element analysis using LinearCohesive
chart = Chart(legend=:top_left, xlabel="σn", ylabel="τ")

for (cmodel, props) in models_props_d
    printstyled("\n$(string(cmodel)):\n\n", color=:yellow, bold=true)
    elem_type = model_elem_d[cmodel]
    
    mapper = RegionMapper()
    add_mapping(mapper, :bulk, MechBulk, LinearElastic, E=E, nu=nu)
    add_mapping(mapper, "interface", elem_type, cmodel; props...)
    
    model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=1.0)
    ana = MechAnalysis(model)
    
    contact = select(model.elems, "interface")
    select( get_ips(contact), tag="jips" )
    log1 = add_logger(ana, :ip, "jips")
    add_monitor(ana, :ip, "jips", (:σn, :τ))
    
    stage = add_stage(ana, nincs=80, nouts=20)
    add_bc(stage, :node, (x==0,y==0), uy=0)
    add_bc(stage, :node, x==0, ux=0)
    add_bc(stage, :node, x==0.2, ux=0.0001)
    
    status = run(ana, autoinc=true, tol=0.05, rspan=0.02, dTmax=0.03, tangent_scheme=:ralston, quiet=false)
    add_series(chart, log1.table["σn"], log1.table["τ"], label=string(cmodel), mark=:circle)

    # break
    @test status.successful
end

save(chart, "contact-models-comparison.pdf")