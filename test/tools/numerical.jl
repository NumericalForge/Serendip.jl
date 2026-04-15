using Serendip, Test

let
    f(x) = exp(x) - 2.0

    x, status = Serendip.findroot_bracket_newton(f, 2.0, tol=1e-10, ftol=1e-12, maxits=80)

    @test Serendip.succeeded(status)
    @test x ≈ log(2.0) atol=1e-8
end


let
    f(x)  = atan(x)
    df(x) = 1.0/(1.0 + x^2)
    
    x, status = Serendip.findroot_bracket_newton(f, df, 10.0, tol=1e-10, ftol=1e-12, maxits=100)
    
    @test Serendip.succeeded(status)
    @test x ≈ 0.0 atol=1e-8
end


let
    f(x)  = x^3 - 2.0*x + 2.0
    df(x) = 3.0*x^2 - 2.0

    _, status = Serendip.findroot_bracket_newton(f, df, 1.0, tol=1e-10, ftol=1e-12, maxits=20)

    @test Serendip.failed(status)
end
