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

reset_affine_index()

@testset "Construction from intervals" begin

    X = 1..3
    X_aff = Affine(X)
    @test X_aff == Affine(2.0, [1.0])

    Y = 2..4
    Y_aff = Affine(Y)
    @test Y_aff == Affine(3.0, [0.0, 1.0])

    @test X_aff + Y_aff == Affine(5.0, [1.0, 1.0])
    @test X_aff * Y_aff == Affine(6.0, [3.0, 2.0, 1.0])
end
