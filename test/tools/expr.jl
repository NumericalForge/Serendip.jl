using Serendip
using Test


@testset "Symbolic arithmetic" begin
    @test evaluate(x; x=3.0) == 3.0

    @test evaluate(x + 2; x=3.0) == 5.0
    @test evaluate(2 + x; x=3.0) == 5.0
    @test evaluate(x - 2; x=3.0) == 1.0
    @test evaluate(5 - x; x=3.0) == 2.0
    @test evaluate(x * y; x=3.0, y=4.0) == 12.0
    @test evaluate(2x; x=3.0) == 6.0
    @test evaluate(x / 2; x=3.0) == 1.5
    @test evaluate(6 / x; x=3.0) == 2.0
    @test evaluate(x^2; x=3.0) == 9.0
    @test evaluate(2^x; x=3.0) == 8.0
    @test evaluate(-x; x=3.0) == -3.0

    @test evaluate(abs(-x); x=2.0) == 2.0
    @test evaluate(sin(x); x=pi/2) ≈ 1.0
    @test evaluate(cos(x); x=0.0) ≈ 1.0
    @test evaluate(tan(x); x=0.0) ≈ 0.0
    @test evaluate(log(exp(x)); x=2.0) ≈ 2.0
    @test evaluate(sqrt(x^2); x=3.0) ≈ 3.0

    @test_throws MethodError x + "unsupported"
end


@testset "Symbolic comparisons and logic" begin
    @test evaluate(x < 2; x=1.0)
    @test evaluate(x <= 1; x=1.0)
    @test evaluate(x == 1; x=1.0 + 0.5e-6)
    @test !evaluate(x == 1; x=1.0 + 2e-6)

    @test evaluate(x > 1; x=1.0 + 2e-6)
    @test evaluate(x >= 1; x=1.0)
    @test evaluate(x != 1; x=1.0 + 2e-6)
    @test !evaluate(x != 1; x=1.0 + 0.5e-6)

    @test evaluate(x == y; x=2.0, y=2.0)
    @test evaluate(and(x > 0, y < 2, true); x=1.0, y=1.0)
    @test evaluate(or(x == 0, y == 1); x=2.0, y=1.0)
    @test evaluate(!(x == 1); x=2.0)
end


@testset "Symbolic affine expressions" begin
    constraint = Symbolic(:(2ux - 3uy = 4))
    terms, rhs = get_affine_terms(constraint)

    @test Dict(terms) == Dict(:ux => 2.0, :uy => -3.0)
    @test rhs == 4.0
    @test check_affine(Serendip.getexpr(constraint), terms, rhs)
    @test_throws ErrorException get_affine_terms(:(ux * uy = 0))
    @test_throws ErrorException get_affine_terms(Expr(:(=), 1, 1))

    @test Serendip.getvars(:(0 <= x < 1)) == [:x]
end


@testset "Symbol replacement" begin
    original = :(a + f(a, b))
    replaced = Serendip._replace_symbol(original, :a => :c)

    @test original == :(a + f(a, b))
    @test replaced == :(c + f(c, b))
end


@testset "Symbolic method ambiguities" begin
    ambiguities = Test.detect_ambiguities(Serendip; recursive=true)
    symbolic_ambiguities = filter(ambiguities) do ambiguity
        occursin("Symbolic", sprint(show, ambiguity))
    end
    @test isempty(symbolic_ambiguities)
end
