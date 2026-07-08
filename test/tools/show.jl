using Serendip
using Test

# Finite element entities
geo = GeoModel()
add_block(geo, [0,0,0], 1,1,0, nx=4, ny=4, shape=:quad9, tag="solids")
mesh = Mesh(geo)

mapper = RegionMapper()
add_mapping(mapper, "solids", MechSolid, LinearElastic, E=100.0, nu=0.2)
model = FEModel(mesh, mapper)
ana = MechAnalysis(model)

log1 = add_logger(ana, :node, (x==1, y==1), "node.dat")
log2 = add_logger(ana, :nodegroup, (y==1), "nodes.dat")

stage = add_stage(ana, nincs=3)
add_bc(stage, :node, (y==0), ux=0, uy=0)
add_bc(stage, :face, (y==1), ty=2)

status = run(ana, quiet=true)


table1 = DataTable("output/node.dat")
table2  = DataTable("output/nodes.dat")

@test status.successful
@test isfile("output/node.dat")
@test isfile("output/nodes.dat")
@test size(table1, 1) > 0
@test size(table2, 1) > 0

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

for obj in (mesh, log1, stage.bcs, model.nodes[1].dofs[1], model.nodes[1].dofs, model.nodes[1], model.nodes, model.elems[1], model.elems, model, table1, table2)
    @test !isempty(sprint(show, obj))
end
