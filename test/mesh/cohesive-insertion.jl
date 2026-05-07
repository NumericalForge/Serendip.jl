using Serendip
using Test

printstyled("\nCohesive insertion (implicit/explicit, whole/partial)\n", color=:blue, bold=true)

function initial_mesh()
    geo = GeoModel()
    add_block(geo, [0.0, 0.0, 0.0], 1.0, 1.0, 0.0, nx=4, ny=4, shape=:quad4)
    return Mesh(geo)
end

function cohesive_pairwise_identity_flags(mesh)
    cohesives = select(mesh.elems, :cohesive)
    flags = Bool[]
    for c in cohesives
        m = div(length(c.nodes), 2)
        same = all(c.nodes[i] === c.nodes[m+i] for i in 1:m)
        push!(flags, same)
    end
    return flags
end

unique_node_objects(nodes) = length(unique(objectid.(nodes)))

@testset "Whole mesh insertion" begin
    mesh = initial_mesh()
    add_cohesive_elements(mesh, tag="coh", implicit=true)

    cohesives = select(mesh.elems, :cohesive)
    @test length(cohesives) == 24
    @test all(cohesive_pairwise_identity_flags(mesh))  # implicit: paired nodes are identical
    @test length(mesh.nodes) == 25

    mesh_exp = initial_mesh()
    add_cohesive_elements(mesh_exp, tag="coh", implicit=false)

    cohesives_exp = select(mesh_exp.elems, :cohesive)
    @test length(cohesives_exp) == 24
    @test all(.!cohesive_pairwise_identity_flags(mesh_exp))  # explicit: paired nodes differ
    @test length(mesh_exp.nodes) == 64
end

@testset "Partial region insertion" begin
    mesh = initial_mesh()
    add_cohesive_elements(mesh, x <= 0.5, tag="coh", implicit=true)

    cohesives = select(mesh.elems, :cohesive)
    @test length(cohesives) == 14
    @test all(cohesive_pairwise_identity_flags(mesh))  # implicit: paired nodes are identical
    @test length(mesh.nodes) == 25

    mesh_exp = initial_mesh()
    add_cohesive_elements(mesh_exp, x <= 0.5, tag="coh", implicit=false)

    cohesives_exp = select(mesh_exp.elems, :cohesive)
    @test length(cohesives_exp) == 14
    @test all(.!cohesive_pairwise_identity_flags(mesh_exp))  # explicit: paired nodes differ
    @test length(mesh_exp.nodes) == 47
end

@testset "Boundary shell generation stays connected after explicit split" begin
    mesh = initial_mesh()
    add_cohesive_elements(mesh, tag="coh", implicit=false, quiet=true)
    add_boundary_shell_elements(mesh, x == 0.0, tag="shell", contact_tag="shell-contact", quiet=true)

    shells = select(mesh.elems, :surface, "shell")
    contacts = select(mesh.elems, :contact, "shell-contact")

    @test length(shells) == 4
    @test length(contacts) == 4
    @test unique_node_objects(get_nodes(shells)) == 5

    for contact in contacts
        shell = contact.couplings[2]
        m = div(length(contact.nodes), 2)
        @test all(contact.nodes[m+i] === shell.nodes[i] for i in 1:m)
    end
end

@testset "Boundary shells without contact stay attached to boundary faces" begin
    mesh = initial_mesh()
    add_boundary_shell_elements(mesh, x == 0.0, tag="shell", quiet=true)

    shells = select(mesh.elems, :surface, "shell")

    @test length(shells) == 4
    for shell in shells
        owner = shell.couplings[1]
        matches_owner_facet = any(owner.shape.facet_idxs) do facet_idxs
            length(shell.nodes) == length(facet_idxs) || return false
            all(shell.nodes[i] === owner.nodes[facet_idxs[i]] for i in eachindex(facet_idxs))
        end
        @test matches_owner_facet
    end
end
