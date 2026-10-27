using Serendip
using Test
using LinearAlgebra

function surface_recovery_test_model(formulation=MechMembrane; material=LinearElastic,
        nx=2, ny=2, shape=:quad9, thickness=0.1, quadrature=0)
    geo = GeoModel()
    add_block(geo, [0.,0.,0.], 2., 1., 0.; nx=nx, ny=ny, shape=shape)
    mesh = Mesh(geo, ndim=3, quiet=true, sort=false)
    plastic = material == LinearElastic ? (;) : (; fy=2.)
    mapper = RegionModel(formulation, material; E=200., nu=0.3, thickness=thickness,
        quadrature=quadrature, plastic...)
    return FEModel(mesh, mapper, quiet=true)
end

function transform_surface_recovery_test!(model, transform)
    for node in model.nodes
        node.coord = transform(node.coord)
    end
    for elem in model.elems
        Serendip.set_quadrature(elem)
        Serendip.elem_init(elem)
    end
end

@testset "Projected surface regression fallback" begin
    for nx in (1,2)
        model = surface_recovery_test_model(nx=nx,ny=1,shape=:quad4,quadrature=1)
        for (i,elem) in enumerate(model.elems), ip in elem.ips
            ip.state.σ = Serendip.Vec3(2.0i,0.,0.)
        end
        Serendip.update_output_data!(model)
        expected = nx == 1 ? 2.0 : 3.0
        @test model.node_fields["σ1′"] ≈ fill(expected,length(model.nodes))
    end

    model = surface_recovery_test_model(nx=1,ny=1,shape=:quad4)
    elem = only(model.elems)
    # Collinear samples cannot support a plane polynomial even when the
    # element nodes span a plane.
    for (i,ip) in enumerate(elem.ips)
        ip.coord = Serendip.Vec3(0.2i,0.,0.)
        ip.state.σ = Serendip.Vec3(3.0i-2,0.,0.)
    end
    Serendip.update_output_data!(model)
    @test model.node_fields["σ1′"] ≈ fill(5.5,length(model.nodes))

    # Exercise recovery of degenerate geometry without invoking element
    # initialization or a mechanical solve on that geometry.
    for (i,node) in enumerate(model.nodes)
        node.coord = Serendip.Vec3(Float64(i),0.,0.)
    end
    Serendip.update_output_data!(model)
    @test model.node_fields["σ1′"] ≈ fill(5.5,length(model.nodes))

    tetrahedron = ([1.,1.,1.],[1.,-1.,-1.],[-1.,1.,-1.],[-1.,-1.,1.])
    for (node,coord) in zip(model.nodes,tetrahedron)
        node.coord = Serendip.Vec3(coord)
    end
    Serendip.update_output_data!(model)
    @test model.node_fields["σ1′"] ≈ fill(5.5,length(model.nodes))
end

@testset "Surface recovery invariance and independent fields" begin
    Q = Matrix(qr([1. 2. 1.; -1. 1. 2.; 2. -1. 1.]).Q)
    for (R, scale, offset) in (
        (Matrix{Float64}(LinearAlgebra.I,3,3),1.,zeros(3)),
        (Q,1.,[10.,-20.,30.]),
        (Q,1e-6,zeros(3)),
        (Q,1e3,[1e6,-2e6,3e6]),
    )
        model = surface_recovery_test_model(material=TensionFieldSoftening)
        transform_surface_recovery_test!(model, p -> R*(scale*p)+offset)
        f(p) = 3+0.5p[1]+0.25p[2]
        plastic(p) = 0.01+0.001p[1]+0.002p[2]
        for elem in model.elems, ip in elem.ips
            localcoord = R'*(ip.coord-offset)/scale
            ip.state.σ = Serendip.Vec3(f(localcoord),0.,0.)
            ip.state.εpa = plastic(localcoord)
        end
        Serendip.update_output_data!(model)
        expected = [f(R'*(node.coord-offset)/scale) for node in model.nodes]
        expectedp = [plastic(R'*(node.coord-offset)/scale) for node in model.nodes]
        @test model.node_fields["σ1′"] ≈ expected atol=1e-8
        @test model.node_fields["σvm′"] ≈ expected atol=1e-8
        @test model.node_fields["σ2′"] ≈ zeros(length(model.nodes)) atol=1e-12
        @test model.node_fields["εp′"] ≈ expectedp atol=1e-10
        @test all(!haskey(model.node_fields,k) for k in ("σvm","ε1′","ε2′","εpa′"))
        @test collect(keys(model.elem_fields)) == ["ruptured"]
        @test all(iszero,model.elem_fields["ruptured"])
    end
end

@testset "Nearly planar, curved and folded surface patches" begin
    for transform in (
        p -> Serendip.Vec3(p[1],p[2],1e-6*p[1]*p[2]),
        p -> Serendip.Vec3(p[1],sin(p[2]),1-cos(p[2])),
        p -> Serendip.Vec3(p[1],min(p[2],0.5),max(0.,p[2]-0.5)),
    )
        model = surface_recovery_test_model()
        transform_surface_recovery_test!(model,transform)
        for elem in model.elems, ip in elem.ips
            ip.state.σ = Serendip.Vec3(4.,1.,0.)
        end
        Serendip.update_output_data!(model)
        @test model.node_fields["σ1′"] ≈ fill(4.,length(model.nodes)) atol=1e-10
        @test model.node_fields["σ2′"] ≈ ones(length(model.nodes)) atol=1e-10
        @test all(isfinite,model.node_fields["σvm′"])
        @test !haskey(model.node_fields,"εp′")
    end
    nearly_planar = surface_recovery_test_model()
    transform_surface_recovery_test!(nearly_planar,p -> Serendip.Vec3(p[1],p[2],1e-6*p[1]*p[2]))
    for elem in nearly_planar.elems, ip in elem.ips
        ip.state.σ = Serendip.Vec3(2+0.2ip.coord[1]+0.3ip.coord[2],0.,0.)
    end
    Serendip.update_output_data!(nearly_planar)
    @test nearly_planar.node_fields["σ1′"] ≈
        [2+0.2n.coord[1]+0.3n.coord[2] for n in nearly_planar.nodes] atol=1e-8
    # Three planes share the corner node; a single scalar estimate is retained.
    points = [0. 0. 0.; 1. 0. 0.; 1. 1. 0.; 0. 1. 0.; 1. 0. 1.; 0. 0. 1.; 0. 1. 1.]
    mesh = Mesh(points,[[1,2,3,4],[1,6,5,2],[1,4,7,6]],fill(:quad4,3),quiet=true)
    for cell in mesh.elems
        cell.role = :surface
    end
    Serendip.synchronize(mesh)
    model = FEModel(mesh,RegionModel(MechMembrane,LinearElastic,E=200.,nu=0.3,thickness=0.1),quiet=true)
    for elem in model.elems, ip in elem.ips
        ip.state.σ = Serendip.Vec3(5.,0.,0.)
    end
    Serendip.update_output_data!(model)
    @test model.node_fields["σ1′"] ≈ fill(5.,length(model.nodes)) atol=1e-10
end

@testset "Critical shell layer and midsurface sampling" begin
    for nth in (2,3,4)
        model = surface_recovery_test_model(MechShell,material=VonMises,
            thickness=10.,quadrature=(3,3,nth))
        expected_vm = 0.
        for elem in model.elems
            depths = unique(ip.R[3] for ip in elem.ips)
            selected = depths[2]
            for ip in elem.ips
                if ip.R[3] == selected
                    # Transverse shear makes this the critical layer, even
                    # though another layer has a higher in-plane σ1 and εpa.
                    ip.state.σ = Serendip.Vec6(1.,0.25,0.,10.,0.,0.4)
                    ip.state.εpa = 0.02
                    expected_vm = Serendip.state_values(elem.cmodel,ip.state)[:σvm]
                else
                    ip.state.σ = Serendip.Vec6(5.,0.,0.,0.,0.,0.)
                    ip.state.εpa = 0.1
                end
            end
            samples = Serendip.recovery_samples(elem)
            @test length(samples) == 9
            @test all(sample.coord[3] == 0. for sample in samples)
            @test all(sample.values[:εp′] == 0.02 for sample in samples)
            @test all(Serendip.state_values(elem.cmodel,ip.state)[:εp] == ip.state.εpa for ip in elem.ips)
            @test isempty(Serendip.elem_vals(elem))
        end
        Serendip.update_output_data!(model)
        radius = hypot(0.375,0.4/√2)
        @test model.node_fields["σ1′"] ≈ fill(0.625+radius,length(model.nodes)) atol=1e-10
        @test model.node_fields["σ2′"] ≈ fill(0.625-radius,length(model.nodes)) atol=1e-10
        @test model.node_fields["σvm′"] ≈ fill(expected_vm,length(model.nodes)) atol=1e-10
        @test model.node_fields["εp′"] ≈ fill(0.02,length(model.nodes)) atol=1e-10
        @test isempty(model.elem_fields)
        for elem in model.elems
            first_depth = elem.ips[1].R[3]
            for ip in elem.ips
                ip.state.σ = Serendip.Vec6(ip.R[3]==first_depth ? -4. : 4.,0.,0.,0.,0.,0.)
                ip.state.εpa = ip.R[3]==first_depth ? 0.03 : 0.07
            end
            # Equal σvm retains the first quadrature layer for all fields.
            @test all(s.values[:σ1′] == 0. && s.values[:σ2′] == -4. && s.values[:εp′] == 0.03
                for s in Serendip.recovery_samples(elem))
        end
    end
    elastic = surface_recovery_test_model(MechShell)
    Serendip.update_output_data!(elastic)
    @test !haskey(elastic.node_fields,"εp′")
end

@testset "Surface fields do not merge with solid fields" begin
    geo=GeoModel()
    add_block(geo,[0,0,0],1,1,1,nx=1,ny=1,nz=1,shape=:hex8,tag="solid")
    mesh=Mesh(geo,quiet=true)
    add_boundary_shell_elements(mesh,z==1,tag="membrane",quiet=true)
    mapper=RegionMapper()
    add_mapping(mapper,"solid",MechSolid,VonMises,E=200.,nu=0.3,fy=2.)
    add_mapping(mapper,"membrane",MechMembrane,TensionFieldSoftening,E=200.,nu=0.3,fy=2.,thickness=0.1)
    model=FEModel(mesh,mapper,quiet=true)
    for elem in model.elems, (i,ip) in enumerate(elem.ips)
        if elem.role == :solid
            ip.state.σ=Serendip.Vec6(80.,0.,0.,0.,0.,0.)
            ip.state.εpa=0.1
        else
            ip.state.σ=Serendip.Vec3(3.,1.,0.)
            ip.state.εpa=0.02
            ip.state.ruptured=iseven(i)
        end
    end
    Serendip.update_output_data!(model)
    @test model.node_fields["σvm"] ≈ fill(80.,length(model.nodes)) atol=1e-10
    for n in model.nodes
        if n.coord[3] == 1.
            @test model.node_fields["σ1′"][n.id] ≈ 3.
            @test model.node_fields["σ2′"][n.id] ≈ 1.
            @test model.node_fields["σvm′"][n.id] ≈ √7
            @test model.node_fields["εp′"][n.id] ≈ 0.02
        else
            @test model.node_fields["σvm′"][n.id] == 0.
        end
    end
    membrane=only(select(model,:element,"membrane"))
    @test Serendip.elem_vals(membrane)[:ruptured] ≈ 0.5
    @test collect(keys(Serendip.elem_vals(membrane))) == [:ruptured]
    @test Serendip.state_values(membrane.cmodel,membrane.ips[1].state)[:σ1] == 3.
    NV, node_fields, EV, elem_fields = Serendip.get_node_and_elem_vals(model)
    @test NV[:,findfirst(==(:σvm′),node_fields)] ≈ model.node_fields["σvm′"]
    @test !(:σvm in elem_fields)
    mktempdir() do dir
        filename = joinpath(dir,"surface.vtu")
        save(model,filename,quiet=true)
        exported = Mesh(filename,quiet=true)
        @test exported.node_fields["σ1′"] ≈ model.node_fields["σ1′"]
        @test exported.node_fields["εp′"] ≈ model.node_fields["εp′"]
        @test !haskey(exported.elem_fields,"σvm")
    end
    plot = DomainPlot()
    add_plot(plot,model; field="σ1′")
    Serendip.configure!(plot)
    @test plot.layers[1].field_kind == :node
end
