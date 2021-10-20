using AffineArithmetic, IntervalArithmetic, Polynomials, StaticArrays
using Test

using AffineArithmetic: Aff

@testset "Construction from intervals" begin
    X = 1..3
    X_a = Affine(X)
    X_a isa Affine

    Y = 2..4
    Y_a = Affine(Y)
    @test Y_a.range == 2..4
    @test Y_a.affine == Aff(3.0, SVector{1, Float64}(1.0), 0..0)

    sum_a = X_a + Y_a
    @test sum_a.range == 3..7
    @test sum_a.affine == Aff(5.0, SVector{1, Float64}(2.0), 0..0)

    prod_a = X_a * Y_a
    @test prod_a.range == 2..12
    @test prod_a.affine == Aff(6.0, SVector{1, Float64}(5.0), 0..1)
end

@testset "Small powers and range" begin
    X = 1..3
    Y = Affine(X)

    f(x) = x^2 - x + 1

    @test range(f(X)) == -1..9
    @test range(f(Y)) == 0..7   # affine is a little better


    g(x) = (x - 1)^3

    @test range(g(X)) == 0..8
    @test range(g(Y)) == 0..8   # affine gives the same
end

reset_affine_index()

@testset "Matrix multiplication" begin
    A = 0.5 * [1 2; -1 1]

    X = Affine(-1..1)
    Y = Affine(-1..1)

    XX = [X, Y]
    @test XX == [Affine(Aff(0.0, SVector{1, Float64}(1.0), 0..0), -1..1),
                 Affine(Aff(0.0, SVector{1, Float64}(1.0), 0..0), -1..1)]

    @test A * XX == [Affine(Aff(0.0, SVector{1, Float64}(1.5), 0..0), -1.5 .. 1.5),
                     Affine(Aff(0.0, SVector{1, Float64}(0.0), 0..0), 0..0)]
end

reset_affine_index()

@testset "Range of polynomial" begin
    # example from Rump 2015
    p = Polynomial([-3, 1])   # x - 3
    p2 = p^8

    x = 4 ± 1e-4
    y = Affine(x)

    @test (-70..70) ⊆ range(p2(x))
    @test range(p2(y)) ⊆ 0.998..1.002   # affine is extremely much better!
end

@testset "Constructors with intervals" begin
    X = Affine(1..2)
    @test X isa Affine
    @test X - X == Affine(0..0)

    x, y = affine(1..2, 3..4)
    @test interval(x + y - x - y) == 0..0
end
