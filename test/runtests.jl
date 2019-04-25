using AffineArithmetic
using IntervalArithmetic
using Polynomials

using Test

@testset "Constructor" begin
    C = Affine(0.0)
    @test C == Affine(0.0, Float64[])

    C = Affine(1.0)
    @test C == Affine(1.0, Float64[])

    C = Affine(1.0, [3.0, 4.0])
    @test C.c == 1.0
    @test C.γ == [3.0, 4.0]

    reset_affine_index()
end

reset_affine_index()

@testset "Construction from intervals" begin

    X = 1..3
    X_a = Affine(X)
    @test X_a == Affine(2.0, [1.0])

    Y = 2..4
    Y_a = Affine(Y)
    @test Y_a == Affine(3.0, [0.0, 1.0])

    @test X_a + Y_a == Affine(5.0, [1.0, 1.0])
    @test X_a * Y_a == Affine(6.0, [3.0, 2.0, 1.0])
end

@testset "Small powers and range" begin

    X = 1..3
    Y = Affine(X)

    f(x) = x^2 - x + 1

    @test range(f(X)) == -1..9
    @test range(f(Y)) == -1..7   # affine is a little better


    g(x) = (x - 1)^3

    @test range(g(X)) == 0..8
    @test range(g(Y)) == -6..8   # affine is significantly worse
end

reset_affine_index()

@testset "Matrix multiplication" begin
    A = 0.5 * [1 2; -1 1]

    X = Affine(-1..1)
    Y = Affine(-1..1)

    XX = [X, Y]
    @test XX == [Affine(0.0, [1.0]), Affine(0.0, [0.0, 1.0])]

    @test A * XX == [Affine(0.0, [0.5, 1.0]), Affine(0.0, [-0.5, 0.5])]
end


reset_affine_index()

@testset "Range of polynomial" begin
    # example from Rump 2015
    p = Poly([-3, 1])   # x - 3
    p2 = p^8

    x = 4 ± 1e-4
    y = Affine(x)

    @test (-70..70) ⊆ range(p2(x))
    @test range(p2(y)) ⊆ 0.998..1.002   # affine is extremely much better!
end
