using Serendip, LinearAlgebra
using Test

@announced_testset "Shape functions" begin
    for shape in Serendip.ALL_ISO_SHAPES
        @announced_testset "$(shape.kind)" begin
            n = shape.npoints
            In = Matrix{Float64}(I, n, n)

            RR = [shape.nat_coords[i, :] for i in 1:n]
            NN = shape.func.(RR)
            II = hcat(NN...)
            @test II ≈ In atol=1e-10

            Q = shape.quadrature[0]
            nip = length(Q)
            RR = [Q[i].coord for i in 1:nip]
            NN = shape.func.(RR)
            @test sum(sum(NN)) ≈ nip atol=1e-10
        end
    end
end

@announced_testset "Shape function derivatives" begin
    for shape in Serendip.ALL_ISO_SHAPES
        @announced_testset "$(shape.kind)" begin
            n = shape.npoints
            ndim = shape.ndim
            RR = [shape.nat_coords[i, :] for i in 1:n]
            f = shape.func
            Id = Matrix{Float64}(I, ndim, ndim)
            δ = 1e-8

            for R in RR
                RI = R .+ Id * δ
                fR = f(R)
                D = zeros(n, ndim)
                for i in 1:ndim
                    Di = 1 / δ * (f(RI[:, i]) - fR)
                    D[:, i] = Di
                end

                @test D ≈ shape.deriv(R) atol=1e-6
            end
        end
    end
end
