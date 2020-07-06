using SimpleRoots, Test

using SimpleRoots: QuadraticError

@testset "test bracketing" begin
    @testset "sin" begin
        @test findzero(sin, (0.0, pi); atol=1e-30)[1] == 0.0
        @test findzero(sin, Secant(0.0, pi); atol=1e-30)[1] == 0.0
        @test findzero(sin, Brent(0.0, pi); atol=1e-30)[1] == 0.0
        @test findzero(sin, Bisection(0.0, pi); atol=1e-100, max_iter=500)[1] ≈ 0.0 atol = 100
    end

    @testset "sin(x) - x / 2" begin
        func1(x) = sin(x) - x / 2
        @test findzero(func1, [0.5pi, pi]; atol=1e-30)[1] == 1.895494267033981
        @test findzero(func1, Bisection([0.5pi, pi]); atol=1e-30)[1] == 1.895494267033981
        @test findzero(func1, Secant([0.5pi, pi]); atol=1e-30)[1] == 1.895494267033981
        @test findzero(func1, Brent([0.5pi, pi]); atol=1e-30)[1] == 1.895494267033981
    end
end

@testset "test quadratic polynomials" begin
    @test quad(1.0, 3.0, -4.0) == (-4.0, 1.0)
    @test quad(Both(), 1.0, 3.0, -4.0) == (-4.0, 1.0)
    @test quad(Lower(), 1.0, 3.0, -4.0) == -4.0
    @test quad(Upper(), 1.0, 3.0, -4.0) == 1.0
    @test_throws QuadraticError quad(1.0, 1.0, 1.0)
    @test_throws QuadraticError quad(0.0, 0.0, 1.0)
end
