using Serendip

geo = GeoModel()
add_block(geo, [0.0, 0.0], [1.0, 10.0], nx=1, ny=10, tag="bulks")

mesh = Mesh(geo)
select(mesh, :element, tag="bulks")

# mplot(mesh, "mesh.pdf")

mapper = RegionMapper()
add_mapping(mapper, "bulks", MechBulk, LinearElastic, E=2e6, nu=0.2, rho=15.0)

model = FEModel(mesh, mapper)
ana = MechModalAnalysis(model)

stage = add_stage(ana)
add_bc(stage, :node, (y==0.0), ux=0, uy=0)

run(ana, nmodes=5)

# mats = [ :bulks => MechBulk => LinearElastic => (E=2e6, nu=0.2, rho=15.0) ]

# model = FEModel(mesh, mats, )
# ana = MechModalAnalysis(model)

# bcs = [
#     :(y==0.0) => NodeBC(ux=0, uy=0)
# ]

# addstage!(ana, bcs)

# solve!(ana, quiet=false, nmodes=5)

