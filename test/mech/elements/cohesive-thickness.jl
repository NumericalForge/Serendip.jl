using Serendip
using Test
using LinearAlgebra

const ECOH = 27.0e6
const NUCOH = 0.2
const ZETACOH = 5.0

function build_cohesive_model(; left_thickness, right_thickness, model_thickness)
    geo = GeoModel(quiet=true)
    add_block(geo, [0.0, 0.0], 1.0, 1.0, 0.0; nx=1, ny=1, shape=:quad4, tag="left")
    add_block(geo, [1.0, 0.0], 1.0, 1.0, 0.0; nx=1, ny=1, shape=:quad4, tag="right")

    mesh = Mesh(geo, quiet=true)
    add_cohesive_elements(mesh, tag="interface", quiet=true)

    mapper = RegionMapper()
    add_mapping(mapper, "left", MechSolid, LinearElastic, E=ECOH, nu=NUCOH, thickness=left_thickness)
    add_mapping(mapper, "right", MechSolid, LinearElastic, E=ECOH, nu=NUCOH, thickness=right_thickness)
    add_mapping(mapper, "interface", MechCohesive, LinearCohesive, E=ECOH, nu=NUCOH, zeta=ZETACOH)

    model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=model_thickness, quiet=true)
    joint = only(select(model, :element, "interface"))
    return model, joint
end


function build_invalid_cohesive_model(; left_thickness, right_thickness, model_thickness)
    geo = GeoModel(quiet=true)
    add_block(geo, [0.0, 0.0], 1.0, 1.0, 0.0; nx=1, ny=1, shape=:quad4, tag="left")
    add_block(geo, [1.0, 0.0], 1.0, 1.0, 0.0; nx=1, ny=1, shape=:quad4, tag="right")

    mesh = Mesh(geo, quiet=true)
    add_cohesive_elements(mesh, tag="interface", quiet=true)

    cell = only(select(mesh, :element, "interface"))
    cell.couplings = [cell.couplings[1]]

    mapper = RegionMapper()
    add_mapping(mapper, "left", MechSolid, LinearElastic, E=ECOH, nu=NUCOH, thickness=left_thickness)
    add_mapping(mapper, "right", MechSolid, LinearElastic, E=ECOH, nu=NUCOH, thickness=right_thickness)
    add_mapping(mapper, "interface", MechCohesive, LinearCohesive, E=ECOH, nu=NUCOH, zeta=ZETACOH)

    return FEModel(mesh, mapper, stress_state=:plane_stress, thickness=model_thickness, quiet=true)
end


@testset "Cohesive thickness from linked solids" begin
    _, joint = build_cohesive_model(left_thickness=0.4, right_thickness=0.9, model_thickness=3.0)

    @test all(th ≈ 0.4 for th in joint.cache.depth_list)
    @test joint.ips[1].state.σ isa Serendip.Vec3
    @test joint.ips[1].state.w isa Serendip.Vec3
end


@testset "Cohesive thickness with expression owner" begin
    right_thickness = :(0.2 + 0.3*y)
    _, joint = build_cohesive_model(left_thickness=0.5, right_thickness=right_thickness, model_thickness=7.0)

    expected = [ evaluate(right_thickness, x=ip.coord.x, y=ip.coord.y, z=ip.coord.z) for ip in joint.ips ]

    @test joint.cache.depth_list ≈ expected
    @test length(unique(round.(joint.cache.depth_list, digits=10))) > 1
end


@testset "Cohesive h matches geometric value for equal thickness" begin
    _, joint = build_cohesive_model(left_thickness=0.4, right_thickness=0.4, model_thickness=3.0)

    @test all(ip.state.h ≈ 1.0 for ip in joint.ips)
end


@testset "Cohesive h uses physical 2d volume area ratio" begin
    _, joint = build_cohesive_model(left_thickness=0.25, right_thickness=0.75, model_thickness=9.0)

    @test all(ip.state.h ≈ 2.0 for ip in joint.ips)
end


@testset "Cohesive h stays element constant with variable thickness" begin
    right_thickness = :(0.2 + 0.3*y)
    _, joint = build_cohesive_model(left_thickness=0.5, right_thickness=right_thickness, model_thickness=7.0)

    expected_h = 0.425/0.35

    @test all(ip.state.h ≈ expected_h for ip in joint.ips)
end


@testset "Cohesive thickness requires two linked solids" begin
    err = try
        build_invalid_cohesive_model(left_thickness=0.4, right_thickness=0.9, model_thickness=3.0)
        nothing
    catch err
        err
    end

    @test err isa ErrorException
    @test occursin("requires exactly two linked owners", sprint(showerror, err))
end


@testset "Cohesive stiffness uses resolved thickness" begin
    _, joint_a = build_cohesive_model(left_thickness=0.25, right_thickness=0.75, model_thickness=9.0)
    _, joint_b = build_cohesive_model(left_thickness=0.25, right_thickness=0.75, model_thickness=2.0)
    _, joint_c = build_cohesive_model(left_thickness=0.50, right_thickness=0.75, model_thickness=1.0)

    Ka, _, _ = Serendip.elem_stiffness(joint_a)
    Kb, _, _ = Serendip.elem_stiffness(joint_b)
    Kc, _, _ = Serendip.elem_stiffness(joint_c)

    ha = joint_a.ips[1].state.h
    hc = joint_c.ips[1].state.h
    expected_ratio = (joint_c.cache.depth_list[1]/hc)/(joint_a.cache.depth_list[1]/ha)

    @test joint_a.cache.depth_list ≈ fill(0.25, length(joint_a.ips))
    @test joint_c.cache.depth_list ≈ fill(0.50, length(joint_c.ips))
    @test Ka ≈ Kb
    @test norm(Kc)/norm(Ka) ≈ expected_ratio
end
