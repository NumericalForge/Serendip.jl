using Serendip, Test

# ❱❱❱ Simple cantilever beam under tip load

# ❱❱❱ Parameters
ℓ = 1.0
n = 4
P = -10
E = 1e4
nu = 0.0
G = E/(2*(1+nu))

# ❱❱❱ Rectangular section
b, h = 0.1, 0.2
I = b*h^3/12
αs = 5/6
u_rect = P*ℓ^3/(3*E*I) + P*ℓ/(αs*G*b*h)

# ❱❱❱ Circular section
r = 0.1
I = π*r^4/4
A = π*r^2
αs = 9/10
u_circ = P*ℓ^3/(3*E*I) + P*ℓ/(αs*G*A)


for section in (:rectangular,:circular)
    for ndim in (2,3)

        for shape in (LIN2, LIN3)
            printstyled("section: ", section, "  ndim: ", ndim, "  shape: ", shape.name, "\n", color= :yellow)

            geo = GeoModel()
            add_block(geo, [0,0,0], ℓ, 0, 0, n=n, shape=shape, tag="beam")
            mesh = Mesh(geo, ndim=ndim)

            mapper = RegionMapper()
            if section == :rectangular
                add_map(mapper, "beam", MechBeam, LinearElastic, E=E, nu=0.0, b=b, h=h)
            else
                add_map(mapper, "beam", MechBeam, LinearElastic, E=E, nu=0.0, A=A)
            end

            model = FEModel(mesh, mapper)
            ana   = MechAnalysis(model)

            if ndim==2
                m = add_monitor(ana, :node, x==1, :uy)
                stage = add_stage(ana)
                add_bc(stage, :node, x==0, rz=0, ux=0, uy=0)
                add_bc(stage, :node, x==ℓ, fy=P)
                run(ana, quiet=false)
                u = m.table[:uy][end]
            else
                m = add_monitor(ana, :node, x==1, :uz)
                stage = add_stage(ana)
                add_bc(stage, :node, x==0, rx=0, ry=0, rz=0, ux=0, uy=0, uz=0)
                add_bc(stage, :node, x==ℓ, fz=P)
                run(ana, quiet=false)
                u = m.table[:uz][end]
            end
            
            if section == :rectangular
                T = @test u ≈ u_rect atol=0.1
            else
                T = @test u ≈ u_circ atol=0.1
            end
            println(T)
            
            println()
        end
    end
end
