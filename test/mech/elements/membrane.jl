using Serendip
using Test
using LinearAlgebra

function membrane_test_model(shape; rotation=Matrix{Float64}(LinearAlgebra.I,3,3), thickness=0.2)
    geo=GeoModel()
    add_block(geo,[0,0,0],2,1,0,nx=1,ny=1,shape=shape)
    mesh=Mesh(geo,ndim=3,quiet=true,sort=false)
    model = FEModel(mesh,RegionModel(MechMembrane,LinearElastic,E=200.,nu=0.3,
        thickness=thickness,rho=3.),quiet=true)
    # FEModel copies nodes through the rounded mesh coordinates. Rotate after
    # that copy to keep these planar patch tests planar to machine precision.
    for node in model.nodes
        node.coord=rotation*node.coord
    end
    foreach(Serendip.elem_init,model.elems)
    for elem in model.elems
        C = Serendip.get_coords(elem)
        for ip in elem.ips
            ip.coord = Serendip.Vec3(C'*elem.shape.func(ip.R))
        end
    end
    return model
end

@testset "MechMembrane kinematics, integration, loads and mass" begin
    @test_throws ArgumentError MechMembrane(thickness=0)
    @test_throws ArgumentError MechMembrane(thickness=1,rho=-1)
    Q = Matrix(qr([1. 2. 1.; -1. 1. 2.; 2. -1. 1.]).Q)
    H = [0.001 0.0002 0.; 0.0002 -0.0001 0.; 0. 0. 0.]
    ee = [H[1,1],H[2,2],√2*H[1,2]]
    D = Serendip.membrane_elastic_matrix(200.,0.3)
    σvm = Serendip.membrane_state_values(D*ee,Serendip.Vec3(ee))[:σvm]
    for shape in (:tri3,:tri6,:quad4,:quad8,:quad9), R in (Matrix{Float64}(LinearAlgebra.I,3,3),Q)
        model=membrane_test_model(shape,rotation=R)
        totalenergy=totalarea=0.
        for elem in model.elems
            @test all(length(n.dofs)==3 for n in elem.nodes)
            @test length(elem.ips) == length(Serendip.get_ip_coords(elem.shape))
            K,_,_=Serendip.elem_stiffness(elem)
            @test K ≈ K' atol=1e-12
            localcoords=R'*Serendip.get_coords(elem)'
            U=vec(R*H*localcoords)
            F,_,status=Serendip.elem_internal_forces(elem,U)
            @test status.successful
            @test F ≈ K*U atol=1e-12
            for ip in elem.ips
                @test Serendip.state_values(elem.cmodel,ip.state)[:σvm] ≈ σvm rtol=1e-12
            end
            totalenergy+=dot(U,K*U)
            area=sum(elem.cache.weights)
            totalarea+=area
            # Translation, rigid in-plane rotation, and arbitrary transverse
            # nodal displacements produce no linearized membrane strain.
            for shift in ([1.,0.,0.],[0.,1.,0.],[0.,0.,1.])
                @test norm(K*repeat(shift,length(elem.nodes))) < 1e-10
            end
            normal=R[:,3]
            transverse=vcat([i*normal for i in eachindex(elem.nodes)]...)
            @test norm(K*transverse) < 1e-10
            spin=vec(R*[0. -1. 0.; 1. 0. 0.; 0. 0. 0.]*localcoords)
            @test norm(K*spin) < 1e-10
            M,_,_=Serendip.elem_mass(elem)
            @test M ≈ M' atol=1e-12
            for d in 1:3
                @test sum(M[d:3:end,d:3:end]) ≈ 3*0.2*area
            end
            face=Serendip.Cell(elem.shape,:surface,elem.nodes,owner=elem)
            pressure,_=Serendip.distributed_bc(elem,face,0.,:tn,2.)
            resultant=vec(sum(reshape(pressure,3,:),dims=2))
            @test norm(resultant) ≈ 2area
            body,_=Serendip.body_c(elem,:wz,4.)
            @test sum(body[3:3:end]) ≈ 4*0.2*area
            @test isempty(Serendip.elem_vals(elem))
        end
        Serendip.update_output_data!(model)
        @test all(isapprox(v,σvm; rtol=1e-10) for v in model.node_fields["σvm′"])
        @test totalarea ≈ 2.
        @test totalenergy ≈ 0.2*2*dot(ee,D*ee) rtol=1e-12
    end
    m1=membrane_test_model(:quad4,thickness=0.1)
    m2=membrane_test_model(:quad4,thickness=0.2)
    @test Serendip.elem_stiffness(m2.elems[1])[1] ≈ 2Serendip.elem_stiffness(m1.elems[1])[1]
end

@testset "MechMembrane nodal output names" begin
    model = membrane_test_model(:quad9)
    elem = only(model.elems)
    for ip in elem.ips
        x, y = ip.coord[1:2]
        ip.state.σ = Serendip.Vec3(x, y, √2*x^2)
        ip.state.ε = Serendip.Vec3(2x, 3y, √2*y^2)
    end
    Serendip.update_output_data!(model)
    @test isempty(Serendip.elem_vals(elem))
    @test all(haskey(model.node_fields,k) for k in ("σ1′","σ2′","σvm′"))
    @test !haskey(model.node_fields,"εp′")
    @test all(!haskey(model.elem_fields,k) for k in ("σvm","σ1","σx´x´","εx´x´"))
    @test haskey(Serendip.state_values(elem.cmodel,elem.ips[1].state),:σxx)
end

@testset "MechMembrane static extension" begin
    for material in (LinearElastic,TensionFieldSoftening)
        geo=GeoModel()
        add_block(geo,[0,0,0],2,1,0,nx=2,ny=1,shape=:quad4)
        mesh=Mesh(geo,ndim=3,quiet=true)
        props=material==LinearElastic ? (;E=200.,nu=0.3) : (;E=200.,nu=0.3,fy=2.)
        model=FEModel(mesh,RegionModel(MechMembrane,material;props...,thickness=0.2),quiet=true)
        for elem in model.elems, ip in elem.ips
            material==TensionFieldSoftening && @test ip.state.h ≈ 1.
        end
        mktempdir() do dir
            ana=MechAnalysis(model,outdir=dir)
            stage=add_stage(ana)
            add_bc(stage,:node,z==0,uz=0)
            add_bc(stage,:node,x==0,ux=0)
            add_bc(stage,:node,y==0,uy=0)
            add_bc(stage,:edge,x==2,qx=0.2)
            @test run(ana,quiet=true).successful
            for node in select(model,:node,x==2)
                @test Serendip.get_dof(node,:ux).vals[:ux] ≈ 0.01 rtol=1e-8
            end
        end
    end
    # Post-yield response with free Poisson contraction and global iteration.
    geo=GeoModel()
    add_block(geo,[0,0,0],2,1,0,nx=2,ny=1,shape=:quad4)
    mesh=Mesh(geo,ndim=3,quiet=true)
    model=FEModel(mesh,RegionModel(MechMembrane,TensionFieldSoftening,E=200.,nu=0.3,
        fy=2.,GF=0.2,beta=2.,thickness=0.2),quiet=true)
    mktempdir() do dir
        ana=MechAnalysis(model,outdir=dir)
        stage=add_stage(ana,nincs=20)
        add_bc(stage,:node,z==0,uz=0)
        add_bc(stage,:node,x==0,ux=0)
        add_bc(stage,:node,y==0,uy=0)
        add_bc(stage,:node,x==2,ux=0.04)
        @test run(ana,quiet=true,autoinc=true,tol=1e-10,rtol=0.0).successful
        for elem in model.elems,ip in elem.ips
            @test ip.state.εpa > 0.
            @test !ip.state.ruptured
            @test abs(ip.state.σ[2]) < 1e-7
        end
        @test maximum(model.node_fields["σ1′"]) > 0.0
        @test maximum(model.node_fields["εp′"]) > 0.0
        vals = Serendip.elem_vals(model.elems[1])
        @test collect(keys(vals)) == [:ruptured]
        @test vals[:ruptured] == 0.0
    end
end

@testset "MechMembrane on a solid boundary" begin
    geo=GeoModel()
    add_block(geo,[0,0,0],1,1,1,nx=1,ny=1,nz=1,shape=:hex8,tag="solid")
    mesh=Mesh(geo,quiet=true)
    add_boundary_shell_elements(mesh,z==1,tag="membrane",quiet=true)
    mapper=RegionMapper()
    add_mapping(mapper,"solid",MechSolid,LinearElastic,E=1000.,nu=0.2)
    add_mapping(mapper,"membrane",MechMembrane,LinearElastic,E=200.,nu=0.2,thickness=0.1)
    model=FEModel(mesh,mapper,quiet=true)
    membrane=only(select(model,:element,"membrane"))
    @test all(length(n.dofs)==3 for n in membrane.nodes)
    mktempdir() do dir
        ana=MechAnalysis(model,outdir=dir)
        stage=add_stage(ana)
        add_bc(stage,:node,x==0,ux=0)
        add_bc(stage,:node,y==0,uy=0)
        add_bc(stage,:node,z==0,uz=0)
        add_bc(stage,:node,x==1,ux=0.001)
        @test run(ana,quiet=true).successful
        reaction=sum(Serendip.get_dof(n,:ux).vals[:fx] for n in select(model,:node,x==1))
        @test reaction ≈ 1.02 rtol=1e-8
    end
end

@testset "MechMembrane coupled through a contact interface" begin
    geo=GeoModel()
    add_block(geo,[0,0,0],1,1,1,nx=1,ny=1,nz=1,shape=:hex8,tag="solid")
    mesh=Mesh(geo,quiet=true)
    add_boundary_shell_elements(mesh,z==1,tag="membrane",contact_tag="interface",quiet=true)
    mapper=RegionMapper()
    add_mapping(mapper,"solid",MechSolid,LinearElastic,E=1000.,nu=0.2)
    add_mapping(mapper,"membrane",MechMembrane,TensionFieldSoftening,E=200.,nu=0.2,fy=2.,thickness=0.1)
    add_mapping(mapper,"interface",MechContact,LinearInterface,kn=1e5,ks=1e5)
    model=FEModel(mesh,mapper,quiet=true)
    mktempdir() do dir
        ana=MechAnalysis(model,outdir=dir)
        stage=add_stage(ana)
        add_bc(stage,:node,x==0,ux=0)
        add_bc(stage,:node,y==0,uy=0)
        add_bc(stage,:node,z==0,uz=0)
        add_bc(stage,:face,"membrane",tz=-1.)
        @test run(ana,quiet=true,tol=1e-10,rtol=0.).successful
        interface=only(select(model,:element,"interface"))
        @test all(ip.state.σ[1] ≈ -1.0 for ip in interface.ips)
        membrane=only(select(model,:element,"membrane"))
        @test all(Serendip.get_dof(n,:uz).vals[:uz] < 0.0 for n in membrane.nodes)
    end
end
