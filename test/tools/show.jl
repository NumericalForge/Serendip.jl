using Serendip

# Finite element entities
geo = GeoModel()
add_block(geo, [0,0], [1,1], nx=4, ny=4, shape=QUAD9, tag="solids")
mesh = Mesh(geo)

mapper = RegionMapper()
add_mapping(mapper, "solids", MechBulk, LinearElastic, E=100.0, nu=0.2)
model = FEModel(mesh, mapper)
ana = MechAnalysis(model)

log1 = add_logger(ana, :node, (x==1, y==1), "node.table")
log2 = add_logger(ana, :nodegroup, (y==1), "nodes.table")

stage = add_stage(ana, nincs=3)
add_bc(stage, :node, (y==0), ux=0, uy=0)
add_bc(stage, :face, (y==1), ty=2)

run(ana)


table1 = DataTable("output/node.table")
table2  = DataTable("output/nodes.table")

println(mesh)
println(log1)
println(stage.bcs)
println(model.nodes[1].dofs[1])
println(model.nodes[1].dofs)
println(model.nodes[1])
println(model.nodes)
println(model.elems[1])
println(model.elems)
println(model)

println(table1)
println(table2)
