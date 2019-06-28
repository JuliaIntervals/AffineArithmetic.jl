using StaticArrays
using IntervalArithmetic
import IntervalArithmetic: interval

using LinearAlgebra

import Base: +, *, ^, -

"""
Affine form with center `c`, affine components `γ` and error `Δ`.

Variant where Δ is an interval
"""
struct AF{N,T<:AbstractFloat}
    c::T   # mid-point
    γ::SVector{N,T}  # affine terms
    Δ::Interval{T}   # error term
end


function Base.show(io::IO, C::AF{N,T}) where {N,T}
    print(io, "⟨", C.c, "; ", C.γ, "; ", C.Δ, "⟩")
end

# ==(C::Affine, D::Affine) = C.c == D.c && C.γ == D.γ

"""
Make an `AF` based on an interval, which is number `i` of `n` total variables.
"""
function AF(X::Interval{T}, n, i) where {T}
    c = mid(X)
    r = radius(X)

    γ = SVector(ntuple(j->i==j ? r : zero(r), n))

    return AF(c, γ, Interval{T}(0))
end

+(x::AF{N,T}, y::AF{N,T}) where {N,T} = AF(x.c + y.c, x.γ .+ y.γ, x.Δ + y.Δ)

-(x::AF{N,T}, y::AF{N,T}) where {N,T} = AF(x.c - y.c, x.γ .- y.γ, x.Δ - y.Δ)


interval(C::AF) = C.c + sum(abs.(C.γ))*(-1..1) + C.Δ


function *(x::AF{N,T}, y::AF{N,T}) where {N,T}
    c = x.c * y.c

    γ = x.c .* y.γ + y.c .* x.γ

    Δ = (x.γ ⋅ y.γ) * (0..1)  # ϵ_i^2

    if N > 1
        Δ += sum(x.γ[i] * y.γ[j] for i in 1:N, j in 1:N if i ≠ j) * (-1..1)  # ϵ_i * ϵ_j
    end

    Δ += (x.c + sum(abs.(x.γ))*(-1..1)) * y.Δ
    Δ += (y.c + sum(abs.(y.γ))*(-1..1)) * x.Δ

    Δ += x.Δ * y.Δ

    return AF(c, γ, Δ)

end

*(x::AF, α::Real) = AF(α*x.c, α.*x.γ, α*x.Δ)
*(α::Real, x::AF) = x * α

+(x::AF, α::Real) = AF(α+x.c, x.γ, x.Δ)
+(α::Real, x::AF) = x + α

-(x::AF) = AF(-x.c, .-(x.γ), -x.Δ)
-(x::AF, α::Real) = AF(x.c - α, x.γ, x.Δ)
-(α::Real, x::AF) = α + (-x)

/(x::AF, α::Real) = AF(x.c/α, x.γ/α, x.Δ/α)

function ^(x::AF, n::Integer)

    invert = false

    if n < 0
        invert = true
        n = -n
        @show n
    end

    result = Base.power_by_squaring(x, n)

    if invert
        result = inv(result)
    end

    return result
end

Base.literal_pow(::typeof(^), x::AF, ::Val{p}) where {T,p} = x^p

x = AF{2,Float64}(0.0, SVector(1.0, 0.0), 0..0)
y = AF{2,Float64}(0.0, SVector(0.0, 1.0), 0..0)


x = AF(3..5, 2, 1)
y = AF(2..4, 2, 2)
#
# 3-x
# interval(3-x)
#
# x + y
#
#
# interval(x+y)
#
# x * y
# interval(x * y)
#
# interval(x * y)
# interval(x) * interval(y)
#
# z = AF(-1..1, 1, 1)
# z^2
# interval(z^2)
#
# using Polynomials
#
# p = Poly([-3, 1])
# p2 = p^8
#
# x = 4 ± 1e-4
# y = AF(x, 1, 1)
#
# interval(y)
# interval(p2(x))
# interval(p2(y))
#
# @time interval(p2(y))
#
#
# f( (x, y) ) = [x^3 + y, (x - y)^2]
#
# X = IntervalBox(-1..1, -1..1)
#
# f(X)
#
# xx = AF(X[1], 2, 1)
# yy = AF(X[2], 2, 2)
#
# interval.(f((xx, yy)))
#
# f(X)
#
#
#
#
# x = AF(4..6, 1, 1)    # example from Messine
# f(x) = x * (10 - x)
#
# f(x)
# interval(f(x))
#
# interval(10*x - x^2)

"General formula for affine approximation of nonlinear functions"
function affine_approx(x::AF, α, ζ, δ)

    c = α * x.c + ζ
    γ = α .* x.γ
    δ += α * x.Δ  # interval

    return AF(c, γ, δ)
end

function Base.sqrt(x::AF, X=interval(x))

    a, b = X.lo, X.hi

    # @show a, b

    # min-range:  de Figuereido book, pg. 64
    α = 1 / (2*√b)
    ζ = (√b) / 2
    δ = ( (√(b) - √(a))^2 / (2*√b) ) * (-1..1)

    return affine_approx(x, α, ζ, δ)
end

function Base.inv(x::AF, X=interval(x))

    a, b = X.lo, X.hi

    # @show a, b

    # min-range:  de Figuereido book, pg. 70
    α = -1 / (b^2)
    d = interval(1/b - α*b, 1/a - α*a)
    ζ = mid(d)
    δ = radius(d) * (-1..1)

    if a < 0
        ζ = -ζ
    end

    return affine_approx(x, α, ζ, δ)
end
