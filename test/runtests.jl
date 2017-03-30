using AffineArithmetic
using ValidatedNumerics

using Base.Test

@testset "Constructor" begin
    C = Affine(0.0)
    @test C == Affine(0.0, Float64[])

    C = Affine(1.0)
    @test C == Affine(1.0, Float64[])

    C = Affine(1.0, [3.0, 4.0])
    @test C.c == 1.0
    @test C.γ == [3.0, 4.0]
end

reset_aine_index()

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


X = 1..3
    X_a = Affine(X)

    Y = X_a^2 - 2X_a + one(X_a)
    r1 = range(Y)

    Z = (X_a - one(X_a))^2
    range(Z)

    Y3 = (X - 1)^2

@testset "range" begin
    X = 1..3
    X_a = Affine(X)

    Y1 = X_a^2 - 2X_a + one(X_a)
    r1 = range(Y)

    Y2 = (X_a - one(X_a))^2
    r2 = range(Z)

    @test r1 == (-2..4)
    @test r2 == r1

    Y3 = (X - 1)^2
    @test Y3 == 0..4

    @test Y3 ⊆ range(Y1)

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
