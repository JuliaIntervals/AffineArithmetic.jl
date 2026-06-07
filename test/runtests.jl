using AffineArithmetic, IntervalArithmetic, Polynomials, StaticArrays
using Test

using AffineArithmetic: Aff

@testset "Construction from intervals" begin
    X = interval(1, 3)
    X_a = Affine(X)
    X_a isa Affine

    Y = interval(2, 4)
    Y_a = Affine(Y)
    @test isequal_interval(Y_a.range, interval(2, 4))
    @test Y_a.affine == Aff(3.0, SVector{1, Float64}(1.0), interval(0))

    sum_a = X_a + Y_a
    @test isequal_interval(sum_a.range, interval(3, 7))
    @test sum_a.affine == Aff(5.0, SVector{1, Float64}(2.0), interval(0))

    prod_a = X_a * Y_a
    @test isequal_interval(prod_a.range, interval(2, 12))
    @test prod_a.affine == Aff(6.0, SVector{1, Float64}(5.0), interval(0, 1))
end

@testset "Construction from numbers" begin
    @test Affine(1.5) == Affine(interval(1.5))
end

@testset "Basic functionality" begin
    X = Affine(interval(1, 3))
    @test eltype(X) == Float64
    T = typeof(X)
    @test zero(T) == Affine(interval(0))
    @test one(T) == Affine(interval(1))
end

@testset "Small powers and range" begin
    X = interval(1, 3)
    Y = Affine(X)

    f(x) = x^2 - x + 1

    @test isequal_interval(range(f(X)), interval(-1, 9))
    @test isequal_interval(range(f(Y)), interval(0, 7))   # affine is a little better

    @test Y^2 == Affine(Aff(4.0, SVector{1}([4.0]), interval(0, 1)), interval(1, 9))
    @test (Y.affine)^2 == Aff(4.0, SVector{1}([4.0]), interval(0, 1))
    @test_broken !isnai(((Y.affine)^(-2)).Δ)
    @test_broken (Y.affine)^(-1)


    g(x) = (x - 1)^3

    @test isequal_interval(range(g(X)), interval(0, 8))
    @test isequal_interval(range(g(Y)), interval(0, 8))   # affine gives the same
end

reset_affine_index()

@testset "Matrix multiplication" begin
    A = 0.5 * [1 2; -1 1]

    X = Affine(interval(-1, 1))
    Y = Affine(interval(-1, 1))

    XX = [X, Y]
    @test XX == [Affine(Aff(0.0, SVector{1, Float64}(1.0), interval(0)), interval(-1, 1)),
                 Affine(Aff(0.0, SVector{1, Float64}(1.0), interval(0)), interval(-1, 1))]

    @test A * XX == [Affine(Aff(0.0, SVector{1, Float64}(1.5), interval(0)), interval(-1.5, 1.5)),
                     Affine(Aff(0.0, SVector{1, Float64}(0.0), interval(0)), interval(0))]
end

reset_affine_index()

@testset "Range of polynomial" begin
    # example from Rump 2015
    p = Polynomial([-3, 1])   # x - 3
    p2 = p^8

    x = interval(4, 1e-4; format = :midpoint)
    y = Affine(x)

    @test issubset_interval(interval(-70, 70), range(p2(x)))
    @test issubset_interval(range(p2(y)), interval(0.998, 1.002))   # affine is extremely much better!
end

@testset "Other operations" begin
    X = Affine(interval(1, 2))
    @test X isa Affine
    @test X - X == Affine(interval(0))

    x, y = affine(interval(1, 2), interval(3, 4))
    @test isequal_interval(interval(x + y - x - y), interval(0))

    @test -X == Affine(Aff(-1.5, SVector{1}([-0.5]), interval(0)), -interval(1, 2))

    A = X.affine
    B = sqrt(A)
    # Note: numbers may be wrong (just copied the result)
    @test B.c ≈ 1.23743686 && B.γ ≈ SVector{1}([0.176776695]) && -inf(B.Δ) ≈ sup(B.Δ) ≈ 0.060660171
    Y = sqrt(X)
    @test Y.affine == B && inf(Y.range) == 1 && sup(Y.range) == sqrt(2)
end
