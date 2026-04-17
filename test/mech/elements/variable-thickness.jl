using Serendip
using Test

E      = 200e6
nu     = 0.2
L      = 4.0
H      = 0.2
th_bot = 0.1
th_top = th_bot/2
qy     = -7.0


function midspan_uy(model::FEModel)
    if model.ctx.ndim == 3
        node = select(model, :node, (x==L/2, y==H, z==th_bot/2))[1]
    else
        node = select(model, :node, (x==L/2, y==H))[1]
    end
    return get_dof(node, :uy).vals[:uy]
end


function run_layered_beam_3d()
    geo = GeoModel(quiet=true)
    add_block(geo, [0.0, 0.0, 0.0], L, H/2, th_bot, nx=20, ny=2, nz=4, shape=:hex20, tag="bottom")
    add_block(geo, [0.0, H/2, th_bot/4], L, H/2, th_top, nx=20, ny=2, nz=2, shape=:hex20, tag="top")
    mesh = Mesh(geo, quiet=true)

    mapper = RegionMapper()
    add_mapping(mapper, "bottom", MechSolid, LinearElastic, E=E, nu=nu)
    add_mapping(mapper, "top", MechSolid, LinearElastic, E=E, nu=nu)

    model = FEModel(mesh, mapper, quiet=true)
    ana = MechAnalysis(model)

    stage = add_stage(ana, nincs=1, nouts=1)
    add_bc(stage, :node, (x==0, y==0), uy=0)
    add_bc(stage, :node, (x==L, y==0), uy=0)
    add_bc(stage, :node, (x==0, y==0, z==th_bot/2), ux=0, uz=0)
    add_bc(stage, :node, (x==L, y==0, z==th_bot/2), uz=0)
    add_bc(stage, :body, x>=0, wy=qy)

    run(ana)

    return model, midspan_uy(model)
end


function run_layered_beam_2d()
    geo = GeoModel(quiet=true)
    add_block(geo, [0.0, 0.0], L, H/2, 0.0, nx=20, ny=2, shape=:quad8, tag="bottom")
    add_block(geo, [0.0, H/2], L, H/2, 0.0, nx=20, ny=2, shape=:quad8, tag="top")
    mesh = Mesh(geo, quiet=true)

    mapper = RegionMapper()
    add_mapping(mapper, "bottom", MechSolid, LinearElastic, E=E, nu=nu, thickness=th_bot)
    add_mapping(mapper, "top", MechSolid, LinearElastic, E=E, nu=nu, thickness=th_top)

    model = FEModel(mesh, mapper, stress_state=:plane_stress, thickness=1.0, quiet=true)
    ana = MechAnalysis(model)

    stage = add_stage(ana, nincs=1, nouts=1)
    add_bc(stage, :node, (x==0, y==0), ux=0, uy=0)
    add_bc(stage, :node, (x==L, y==0), uy=0)
    add_bc(stage, :body, x>=0, wy=qy)

    run(ana)

    return model, midspan_uy(model)
end


model3d, uy3d = run_layered_beam_3d()
model2d, uy2d = run_layered_beam_2d()

@test all(elem.etype.thickness ≈ th_bot for elem in select(model2d, :element, "bottom"))
@test all(elem.etype.thickness ≈ th_top for elem in select(model2d, :element, "top"))

@show uy2d
@show uy3d
@show uy2d/uy3d

@test uy2d ≈ uy3d atol=1e-1

