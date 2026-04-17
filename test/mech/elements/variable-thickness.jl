using Serendip
using Test

E  = 20e6
nu = 0.2
L  = 4.0
qy = -7.0

y1 = 0.20
y2 = 0.40
y3 = 0.60
y4 = 1.00
y5 = 1.20

t1 = 0.40
t2 = 0.20
t3 = 0.40

thickness = :(
    y < $y1 ? $t1 :
    y < $y2 ? $t1 - ($t1 - $t2)/($y2 - $y1)*(y - $y1) :
    y < $y3 ? $t2 :
    y < $y4 ? $t2 + ($t3 - $t2)/($y4 - $y3)*(y - $y3) :
    $t3
)


function midspan_uy(model::FEModel)
    if model.ctx.ndim == 3
        node = select(model, :node, (x==L/2, y==y5, z==0))[1]
    else
        node = select(model, :node, (x==L/2, y==y5))[1]
    end
    return get_dof(node, :uy).vals[:uy]
end


function run_variable_beam_3d()
    geo = GeoModel(size=0.5, quiet=true)

    pts = [
        add_point(geo, [0.0, 0.0, -t1/2]),
        add_point(geo, [0.0, 0.0,    0.0]),
        add_point(geo, [0.0, 0.0,  t1/2]),
        add_point(geo, [0.0, y1,   t1/2]),
        add_point(geo, [0.0, y2,   t2/2]),
        add_point(geo, [0.0, y3,   t2/2]),
        add_point(geo, [0.0, y4,   t3/2]),
        add_point(geo, [0.0, y5,   t3/2]),
        add_point(geo, [0.0, y5,    0.0]),
        add_point(geo, [0.0, y5,  -t3/2]),
        add_point(geo, [0.0, y4,  -t3/2]),
        add_point(geo, [0.0, y3,  -t2/2]),
        add_point(geo, [0.0, y2,  -t2/2]),
        add_point(geo, [0.0, y1,  -t1/2]),
    ]

    surf = add_polygon(geo, pts; tag="beam")
    extrude(geo, surf, [L, 0.0, 0.0], num_elements=[24])
    mesh = Mesh(geo, quadratic=true, quiet=true)

    mapper = RegionModel(MechSolid, LinearElastic, E=E, nu=nu)
    model = FEModel(mesh, mapper, quiet=true)
    ana = MechAnalysis(model)

    stage = add_stage(ana, nincs=1, nouts=1)
    add_bc(stage, :node, (x==0, y==0), uy=0)
    add_bc(stage, :node, (x==L, y==0), uy=0)
    add_bc(stage, :node, (x==0, y==0, z==0), ux=0, uz=0)
    add_bc(stage, :node, (x==L, y==0, z==0), uz=0)
    add_bc(stage, :body, x>=0, wy=qy)

    run(ana)

    return model, midspan_uy(model)
end


function run_variable_beam_2d()
    geo = GeoModel(size=0.2, quiet=true)
    pts = [
        add_point(geo, [0.0, 0.0, 0.0]),
        add_point(geo, [L,   0.0, 0.0]),
        add_point(geo, [L,    y5, 0.0]),
        add_point(geo, [L/2,  y5, 0.0]),
        add_point(geo, [0.0,  y5, 0.0]),
    ]
    add_polygon(geo, pts; tag="beam")
    mesh = Mesh(geo, quadratic=true, quiet=true)

    mapper = RegionModel(MechSolid, LinearElastic, E=E, nu=nu, thickness=thickness)
    model = FEModel(mesh, mapper, stress_state=:plane_stress, quiet=true)
    ana = MechAnalysis(model)

    stage = add_stage(ana, nincs=1, nouts=1)
    add_bc(stage, :node, (x==0, y==0), ux=0, uy=0)
    add_bc(stage, :node, (x==L, y==0), uy=0)
    add_bc(stage, :body, x>=0, wy=qy)

    run(ana)

    return model, midspan_uy(model)
end


_, uy3d = run_variable_beam_3d()
_, uy2d = run_variable_beam_2d()

@show uy2d
@show uy3d
@show uy2d/uy3d

@test uy2d ≈ uy3d rtol=0.02
