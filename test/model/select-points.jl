using Serendip
using Test

mutable struct DummyState <: Serendip.ConstState
    value::Float64
end

function make_point_test_model()
    geo = GeoModel(quiet=true)
    add_block(geo, [0,0,0], 1,1,1, nx=2, ny=1, nz=1, tag="solids")
    mesh = Mesh(geo, quiet=true)
    select(mesh, :element, x <= 0.5, tag="left")
    select(mesh, :element, x >= 0.5, tag="right")

    mapper = RegionMapper()
    add_mapping(mapper, "left", MechSolid, LinearElastic, E=100.0, nu=0.2)
    add_mapping(mapper, "right", MechSolid, LinearElastic, E=100.0, nu=0.2)

    model = FEModel(mesh, mapper, quiet=true)
    ana = MechAnalysis(model)
    return model, ana
end

function capture_stdout(f)
    pipe = Pipe()
    redirect_stdout(pipe) do
        f()
    end
    close(pipe.in)
    return read(pipe, String)
end


printstyled("\nPoint Selection Through select\n", color=:yellow, bold=true)

nodes = Node[
    Node(0.0, 0.0, 0.0, id=1),
    Node(0.0, 0.0, 0.0, id=2),
    Node(1.0, 0.0, 0.0, id=3),
]

@test length(select(nodes, [0.0, 0.0, 0.0])) == 2
@test length(select(nodes, [0.0, 0.0, 0.0]; nearest=true)) == 2
@test isempty(select(nodes, [0.25, 0.0, 0.0]))
@test select(nodes, [0.25, 0.0, 0.0]; nearest=true)[1].id == 1

out = capture_stdout() do
    select(nodes, [0.25, 0.0, 0.0])
end
@test occursin("select: No node found", out)

out = capture_stdout() do
    select(nodes, [0.25, 0.0, 0.0]; quiet=true)
end
@test isempty(out)

owner = Cell(Serendip.LIN2, :line, [Node(0.0, id=1), Node(1.0, id=2)])
ips = Ip[
    Ip([0.0], 1.0, owner, DummyState(0.0)),
    Ip([0.0], 1.0, owner, DummyState(0.0)),
    Ip([0.0], 1.0, owner, DummyState(0.0)),
]
ips[1].coord = Serendip.Vec3(0.5, 0.0, 0.0)
ips[2].coord = Serendip.Vec3(0.5, 0.0, 0.0)
ips[3].coord = Serendip.Vec3(1.0, 0.0, 0.0)
ips[1].id = 1
ips[2].id = 2
ips[3].id = 3

@test length(select(ips, [0.5, 0.0, 0.0])) == 2
@test length(select(ips, [0.5, 0.0, 0.0]; nearest=true)) == 2
@test isempty(select(ips, [0.75, 0.0, 0.0]))
@test select(ips, [0.75, 0.0, 0.0]; nearest=true)[1].id == 1

out = capture_stdout() do
    select(ips, [0.75, 0.0, 0.0]; prefix="custom")
end
@test occursin("custom: No ip found", out)


printstyled("Chained Point Selection", color=:yellow, bold=true)

model, _ = make_point_test_model()

left_node = select(model, :element, "left", :node, [0.0, 0.0, 0.0])
@test length(left_node) == 1
@test left_node[1].coord == Serendip.Vec3(0.0, 0.0, 0.0)

left_ip_exact = select(model, :element, "left", :ip, collect(select(model, :element, "left", :ip, :all)[1].coord))
@test length(left_ip_exact) >= 1
@test all(ip.tag != "left_ips" for ip in select(model, :ip, :all))

selected_ips = select(model, :element, "left", :ip, [0.75, 0.5, 0.5]; nearest=true, tag="left_ips")
@test length(selected_ips) == 1
@test selected_ips[1].tag == "left_ips"

left_candidates = select(model, :element, "left", :ip, :all)
right_candidates = select(model, :element, "right", :ip, :all)
@test selected_ips[1] in left_candidates
@test !(selected_ips[1] in right_candidates)


printstyled("Logger and Monitor Point Selectors", color=:yellow, bold=true)

model, ana = make_point_test_model()
select(model, :element, "left", :ip, tag="left_ips")

logger = add_logger(ana, :ip, ("left_ips", [0.75, 0.5, 0.5]))
@test length(logger.target) == 1
@test logger.target[1] in select(model, :element, "left", :ip, :all)

monitor = add_monitor(ana, :ip, ("left_ips", [0.75, 0.5, 0.5]), :σyy)
@test length(monitor.target) == 1
@test monitor.target[1] in select(model, :element, "left", :ip, :all)

out = capture_stdout() do
    add_logger(ana, :ip, ("left_ips", [0.75, 0.5, 0.5]))
end
@test occursin("add_logger: No ip found", out)

out = capture_stdout() do
    add_monitor(ana, :ip, ("left_ips", [0.75, 0.5, 0.5]), :σyy)
end
@test occursin("add_monitor: No ip found", out)
