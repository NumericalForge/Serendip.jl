using Serendip
using Test

@announced_testset "Shape extrapolation" begin
    for shape in Serendip.ALL_ISO_SHAPES
        @announced_testset "$(shape.kind)" begin
            C = shape.nat_coords
            n = shape.npoints
            V = zeros(n)
            for i in 1:n
                x, y, z = [C[i, :]; 0.0; 0.0]
                V[i] = x + y + z + 1.0
            end

            for (nip, Cip) in shape.quadrature
                nip <= 1 && continue
                (nip == 2 && shape.base_shape == get_shape(:wed6)) && continue

                @announced_testset "nip=$nip" begin
                    P = zeros(nip)
                    for i in 1:nip
                        x, y, z = Cip[i].coord
                        P[i] = x + y + z + 1.0
                    end

                    E = extrapolator(shape, nip)
                    VV = E * P
                    @test V ≈ VV atol=1e-10
                end
            end
        end
    end
end
