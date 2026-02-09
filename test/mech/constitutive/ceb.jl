using Serendip
using Test

# ❱❱❱ Geometry
geo = GeoModel()
add_block(geo, [0.0, 0.0, 0.0], 1.0, 6.0, 1.0, nx=1, ny=5, nz=1, tag="solids")
p1 = add_point(geo, [0.5, 3.0, 0.5])
p2 = add_point(geo, [0.5, 6.0, 0.5])
edge = add_line(geo, p1, p2)
path = add_path(geo, [edge]; tag="bars", interface_tag="joints")


# ❱❱❱ Mesh
mesh = Mesh(geo)
save(mesh, "mesh.vtk")

bar_points = get_nodes(select(mesh, :element, "bars"))
select(bar_points, (y==6), tag="tip")
select(get_nodes(select(mesh, :element, "solids")), tag="fixed")

# ❱❱❱ Finite elements
mapper = RegionMapper()
add_mapping(mapper, "solids", MechBulk, LinearElastic, E=24e3, nu=0.2)
add_mapping(mapper, "bars", MechBar, LinearElastic, E=200e6, A=0.00011)
add_mapping(mapper, "joints", MechBondSlip, CebBondSlip,
            taumax=12.0, taures=3.0, s1=0.001, s2=0.0011, s3=0.004, alpha=0.5,
            ks=(12/0.001)*5, kn=5e3, p=0.25)

model = FEModel(mesh, mapper)
ana = MechAnalysis(model)

select(get_ips(select(model, :element, "joints")), tag="joint_ips")

# add_logger(ana, :node, "tip")
# add_logger(ana, :ipgroup, "joint_ips")
logg = add_logger(ana, :ip, "joint_ips")
add_monitor(ana, :node, "tip", :uy)
add_monitor(ana, :node, "tip", :fy)

stage = add_stage(ana)

add_bc(stage, :node, "fixed", ux=0, uy=0, uz=0)
add_bc(stage, :node, "tip", uy=0.0003)

stage = add_stage(ana)
add_bc(stage, :node, "fixed", ux=0, uy=0, uz=0)
add_bc(stage, :node, "tip", uy=-0.0001)

stage = add_stage(ana)
add_bc(stage, :node, "fixed", ux=0, uy=0, uz=0)
add_bc(stage, :node, "tip", uy=+0.0006)

stage = add_stage(ana)
add_bc(stage, :node, "fixed", ux=0, uy=0, uz=0)
add_bc(stage, :node, "tip", uy=-0.0005)

stage = add_stage(ana)
add_bc(stage, :node, "fixed", ux=0, uy=0, uz=0)
add_bc(stage, :node, "tip", uy=+0.005)

run(ana, autoinc=true, tol=0.01, maxits=3)

# ❱❱❱ Post-processing

table = logg.table
chart = Chart(xlabel=L"$u_r$", ylabel=L"$\tau$", legend=:bottom_right)
add_series(chart, :line, table["s"], table["τ"], mark=:circle, color=:red)
save(chart, "chart.pdf")
