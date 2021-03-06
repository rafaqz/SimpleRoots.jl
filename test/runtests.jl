using SimpleRoots, Test

using SimpleRoots: QuadraticError

@testset "test bracketing" begin
    @testset "sin 0.0 is 0.0" begin
        @test findzero(sin, (-1.0, 1.3); atol=1e-30)[1] ≈ 0.0 atol = 1e-30
        @test findzero(sin, (-1.0, 1.3); atol=1e-30)[2] == true
        @test findzero(sin, Secant(-1.0, 1.3); atol=1e-30)[1] ≈ 0.0 atol = 1e-30
        @test findzero(sin, Secant(-1.0, 1.3); atol=1e-30)[2] == true
        @test findzero(sin, Brent(-1.0, 1.3); atol=1e-30)[1] ≈ 0.0 atol = 1e-30
        @test findzero(sin, Brent(-1.0, 1.3); atol=1e-30)[2] == true
        @test findzero(sin, Bisection(-1.0, 1.3); atol=1e-100, max_iter=500)[1] ≈ 0.0 atol = 1e-30
        @test findzero(sin, Bisection(-1.0, 1.3); atol=1e-100, max_iter=500)[2] == true
    end

    @testset "sin with 1 iteration isn't found" begin
        @test findzero(sin, (-1.0, 1.3); atol=1e-30, max_iter=1)[2] == false
        @test findzero(sin, Secant(-1.0, 1.3); atol=1e-30, max_iter=1)[2] == false
        @test findzero(sin, Brent(-1.0, 1.3); atol=1e-30, max_iter=1)[2] == false
        @test findzero(sin, Bisection(-1.0, 1.3); atol=1e-100, max_iter=1)[2] == false
        # Except bisection when the bracket is symmetrical
        @test findzero(sin, Bisection(-1.0, -1.0); atol=1e-100, max_iter=1)[2] == true
    end

    @testset "sin(x) - x / 2" begin
        func1(x) = sin(x) - x / 2
        @test findzero(func1, [0.5pi, pi]; atol=1e-30)[1] == 1.895494267033981
        @test findzero(func1, Bisection([0.5pi, pi]); atol=1e-30)[1] == 1.895494267033981
        @test findzero(func1, Secant([0.5pi, pi]); atol=1e-30)[1] == 1.895494267033981
        @test findzero(func1, Brent([0.5pi, pi]); atol=1e-30)[1] == 1.895494267033981
    end

    @testset "inference" begin
        @inferred findzero(sin, (-1.0, 1.0); atol=1e-100, max_iter=100)
        @inferred findzero(sin, Secant(-1.0, 1.0); atol=1e-100, max_iter=100)
        @inferred findzero(sin, Brent(-1.0, 1.0); atol=1e-100, max_iter=100)
        @inferred findzero(sin, Bisection(-1.0, 1.0); atol=1e-100, max_iter=100)
    end
end

@testset "test quadratic polynomials" begin
    @test quad(1.0, 3.0, -4.0) == (-4.0, 1.0)
    @test quad(Both(), 1.0, 3.0, -4.0) == (-4.0, 1.0)
    @test quad(Lower(), 1.0, 3.0, -4.0) == -4.0
    @test quad(Upper(), 1.0, 3.0, -4.0) == 1.0

    @test quad(Both(), 0.0, 1.0, 1.0) == (-1.0, -1.0) 
    @test quad(Lower(), 0.0, 1.0, 1.0) == -1.0 
    @test quad(Upper(), 0.0, 1.0, 1.0) == -1.0

    @test quad(Both(), 1.0, 0.0, 0.0) == (0.0, 0.0) 
    @test quad(Lower(), 1.0, 0.0, 0.0) == 0.0 
    @test quad(Upper(), 1.0, 0.0, 0.0) == 0.0

    @test quad(Both(), 0.0, 0.0, 0.0) == (0.0, 0.0) 
    @test quad(Lower(), 0.0, 0.0, 0.0) == 0.0 
    @test quad(Upper(), 0.0, 0.0, 0.0) == 0.0

    @test_throws QuadraticError quad(1.0, 1.0, 1.0)
    @test_throws QuadraticError quad(0.0, 0.0, 1.0)

    @inferred quad(1.0, 3.0, -4.0)
end
