using Serendip
using Test
using LinearAlgebra

const ECONTACT = 27.0e6
const NUCONTACT = 0.2
const KNCONTACT = 2.0e8
const KSCONTACT = 1.8e8

function build_contact_model_2d(; left_thickness, right_thickness, model_thickness)
    geo = GeoModel(quiet=true)
    add_block(geo, [0.0, 0.0], 1.0, 1.0, 0.0; nx=1, ny=1, shape=:quad4, tag="left")
    add_block(geo, [1.0, 0.0], 1.0, 1.0, 0.0; nx=1, ny=1, shape=:quad4, tag="right")

    mesh = Mesh(geo, quiet=true)
    add_contact_elements(mesh, "left", "right", tag="interface", quiet=true)

    mapper = RegionMapper()
    add_mapping(mapper, "left", MechSolid, LinearElastic, E=ECONTACT, nu=NUCONTACT, thickness=left_thickness)
    add_mapping(mapper, "right", MechSolid, LinearElastic, E=ECONTACT, nu=NUCONTACT, thickness=right_thickness)
    add_mapping(mapper, "interface", MechContact, LinearContact, kn=KNCONTACT, ks=KSCONTACT)

    model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=model_thickness, quiet=true)
    contact = only(select(model, :element, "interface"))
    return model, contact
end


function build_boundary_contact_model_2d(; owner_thickness, model_thickness, clear_couplings=false)
    geo = GeoModel(quiet=true)
    add_block(geo, [0.0, 0.0], 1.0, 1.0, 0.0; nx=1, ny=1, shape=:quad4, tag="bulk")

    mesh = Mesh(geo, quiet=true)
    add_boundary_contact_elements(mesh, x==1.0, tag="interface", quiet=true)

    if clear_couplings
        cell = only(select(mesh, :element, "interface"))
        cell.couplings = typeof(cell.couplings)()
    end

    mapper = RegionMapper()
    add_mapping(mapper, "bulk", MechSolid, LinearElastic, E=ECONTACT, nu=NUCONTACT, thickness=owner_thickness)
    add_mapping(mapper, "interface", MechContact, LinearContact, kn=KNCONTACT, ks=KSCONTACT)

    model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=model_thickness, quiet=true)
    contact = only(select(model, :element, "interface"))
    return model, contact
end


function build_contact_model_3d(; model_thickness)
    geo = GeoModel(quiet=true)
    add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 1.0; nx=1, ny=1, nz=1, shape=:hex8, tag="left")
    add_block(geo, [1.0, 0.0, 0.0], 1.0, 1.0, 1.0; nx=1, ny=1, nz=1, shape=:hex8, tag="right")

    mesh = Mesh(geo, quiet=true)
    add_contact_elements(mesh, "left", "right", tag="interface", quiet=true)

    mapper = RegionMapper()
    add_mapping(mapper, "left", MechSolid, LinearElastic, E=ECONTACT, nu=NUCONTACT)
    add_mapping(mapper, "right", MechSolid, LinearElastic, E=ECONTACT, nu=NUCONTACT)
    add_mapping(mapper, "interface", MechContact, LinearContact, kn=KNCONTACT, ks=KSCONTACT)

    model = FEModel(mesh, mapper, thickness=model_thickness, quiet=true)
    contact = only(select(model, :element, "interface"))
    return model, contact
end


@testset "Contact thickness from linked solids" begin
    _, contact = build_contact_model_2d(left_thickness=0.4, right_thickness=0.9, model_thickness=3.0)

    @test contact.cache.depth_list ≈ fill(0.4, length(contact.ips))
end


@testset "Boundary contact thickness follows owner" begin
    owner_thickness = :(0.2 + 0.3*y)
    _, contact = build_boundary_contact_model_2d(owner_thickness=owner_thickness, model_thickness=7.0)

    expected = [evaluate(owner_thickness, x=ip.coord.x, y=ip.coord.y, z=ip.coord.z) for ip in contact.ips]

    @test contact.cache.depth_list ≈ expected
    @test length(unique(round.(contact.cache.depth_list, digits=10))) > 1
end


@testset "Boundary contact falls back to model thickness" begin
    _, contact = build_boundary_contact_model_2d(owner_thickness=0.5, model_thickness=1.7, clear_couplings=true)

    @test contact.cache.depth_list ≈ fill(1.7, length(contact.ips))
end


@testset "Contact stiffness uses resolved thickness" begin
    _, contact_a = build_contact_model_2d(left_thickness=0.25, right_thickness=0.75, model_thickness=9.0)
    _, contact_b = build_contact_model_2d(left_thickness=0.25, right_thickness=0.75, model_thickness=2.0)
    _, contact_c = build_contact_model_2d(left_thickness=0.50, right_thickness=0.75, model_thickness=1.0)

    Ka, _, _ = Serendip.elem_stiffness(contact_a)
    Kb, _, _ = Serendip.elem_stiffness(contact_b)
    Kc, _, _ = Serendip.elem_stiffness(contact_c)

    @test contact_a.cache.depth_list ≈ fill(0.25, length(contact_a.ips))
    @test contact_c.cache.depth_list ≈ fill(0.50, length(contact_c.ips))
    @test Ka ≈ Kb
    @test norm(Kc)/norm(Ka) ≈ 2.0
end


@testset "3d contact ignores model thickness" begin
    _, contact_a = build_contact_model_3d(model_thickness=1.0)
    _, contact_b = build_contact_model_3d(model_thickness=7.0)

    Ka, _, _ = Serendip.elem_stiffness(contact_a)
    Kb, _, _ = Serendip.elem_stiffness(contact_b)

    @test contact_a.cache.depth_list ≈ fill(1.0, length(contact_a.ips))
    @test contact_b.cache.depth_list ≈ fill(1.0, length(contact_b.ips))
    @test Ka ≈ Kb
end
