using Serendip
using Test


function build_single_node_constraint_model()
    geo = GeoModel()
    add_block(geo, [0, 0], 1, 1, 0, nx=1, ny=1, shape=:quad4)
    mesh = Mesh(geo, quiet=true)

    mapper = RegionModel(MechSolid, LinearElastic, E=10.0, nu=0.25)
    model = FEModel(mesh, mapper, stress_state=:plane_strain)
    ana = MechAnalysis(model)

    stage = add_stage(ana, nincs=1)
    add_bc(stage, :face, x==0, ux=0)
    add_bc(stage, :face, y==0, uy=0)
    add_constraint(stage, :node, [1.0, 1.0, 0.0], :(ux-uy=0))
    add_bc(stage, :node, [1.0, 1.0, 0.0], fx=1.0)

    return model, ana
end


function build_inclined_support_constraint_model()
    r = 0.3
    ℓ = 1.0 - r
    α = 45.0

    geo = GeoModel(quiet=true)
    add_block(geo, [
        r                         0.0
        1.0                       0.0
        cosd(α)                   sind(α)
        r*cosd(α)                 r*sind(α)
        (r + ℓ/2)                 0.0
        cosd(α/2)                 sind(α/2)
        ((r + ℓ/2)*cosd(α))       ((r + ℓ/2)*sind(α))
        r*cosd(α/2)               r*sind(α/2)
    ], nx=5, ny=5, shape=:quad8)
    mesh = Mesh(geo, quiet=true)

    mapper = RegionModel(MechSolid, LinearElastic, E=10.0, nu=0.25)
    model = FEModel(mesh, mapper, stress_state=:plane_strain)
    ana = MechAnalysis(model)

    stage = add_stage(ana, nincs=1)
    add_bc(stage, :edge, y==0, uy=0)
    add_constraint(stage, :edge, x==y, :(ux-uy=0))
    add_bc(stage, :edge, x*x + y*y > 0.95, qn=0.12)

    return model, ana, α
end


model1, ana1 = build_single_node_constraint_model()

run(ana1, tol=1e-8, maxits=10, quiet=true)

vals1 = get_values(select(model1, :node, [1.0, 1.0, 0.0])[1])

@test vals1[:ux] ≈ vals1[:uy] atol=1e-8
@test abs(vals1[:ux]) > 1e-8

model2, ana2, α = build_inclined_support_constraint_model()

run(ana2, tol=1e-8, utol=1e-8, maxits=10, quiet=true)

bottom_nodes = select(model2, :node, y==0)
top_nodes = select(model2, :node, x==y)
outer_mid = select(model2, :node, [cosd(α/2), sind(α/2), 0.0])[1]
outer_vals = get_values(outer_mid)

@test !isempty(bottom_nodes)
@test !isempty(top_nodes)
@test all(isapprox(get_values(node)[:uy], 0.0; atol=1e-8) for node in bottom_nodes)
@test all(isapprox(get_values(node)[:ux], get_values(node)[:uy]; atol=1e-8) for node in top_nodes)
@test hypot(outer_vals[:ux], outer_vals[:uy]) > 1e-8

