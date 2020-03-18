using StaticArrays
using IntervalArithmetic
import IntervalArithmetic: interval

using LinearAlgebra

import Base: +, *, ^, -

"""
Affine form with center `c`, affine components `γ` and error `Δ`.

Variant where Δ is an interval
"""
struct Aff{N,T<:AbstractFloat}
    c::T   # mid-point
    γ::SVector{N,T}  # affine terms
    Δ::Interval{T}   # error term
end


function Base.show(io::IO, C::Aff{N,T}) where {N,T}
    print(io, "⟨", C.c, "; ", C.γ, "; ", C.Δ, "⟩")
end

# ==(C::Affine, D::Affine) = C.c == D.c && C.γ == D.γ

"""
Make an `Aff` based on an interval, which is number `i` of `n` total variables.
"""
function Aff(X::Interval{T}, n, i) where {T}
    c = mid(X)
    r = radius(X)

    γ = SVector(ntuple(j->i==j ? r : zero(r), n))

    return Aff(c, γ, Interval{T}(0))
end

+(x::Aff{N,T}, y::Aff{N,T}) where {N,T} = Aff(x.c + y.c, x.γ .+ y.γ, x.Δ + y.Δ)

-(x::Aff{N,T}, y::Aff{N,T}) where {N,T} = Aff(x.c - y.c, x.γ .- y.γ, x.Δ - y.Δ)


interval(C::Aff) = C.c + sum(abs.(C.γ))*(-1..1) + C.Δ


function *(x::Aff{N,T}, y::Aff{N,T}) where {N,T}
    c = x.c * y.c

    γ = x.c .* y.γ + y.c .* x.γ

    Δ = (x.γ ⋅ y.γ) * (0..1)  # ϵ_i^2

    if N > 1
        Δ += sum(x.γ[i] * y.γ[j] for i in 1:N, j in 1:N if i ≠ j) * (-1..1)  # ϵ_i * ϵ_j
    end

    Δ += (x.c + sum(abs.(x.γ))*(-1..1)) * y.Δ
    Δ += (y.c + sum(abs.(y.γ))*(-1..1)) * x.Δ

    Δ += x.Δ * y.Δ

    return Aff(c, γ, Δ)

end

*(x::Aff, α::Real) = Aff(α*x.c, α.*x.γ, α*x.Δ)
*(α::Real, x::Aff) = x * α

+(x::Aff, α::Real) = Aff(α+x.c, x.γ, x.Δ)
+(α::Real, x::Aff) = x + α

-(x::Aff) = Aff(-x.c, .-(x.γ), -x.Δ)
-(x::Aff, α::Real) = Aff(x.c - α, x.γ, x.Δ)
-(α::Real, x::Aff) = α + (-x)

/(x::Aff, α::Real) = Aff(x.c/α, x.γ/α, x.Δ/α)

function ^(x::Aff, n::Integer)

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

Base.literal_pow(::typeof(^), x::Aff, ::Val{p}) where {T,p} = x^p

x = Aff{2,Float64}(0.0, SVector(1.0, 0.0), 0..0)
y = Aff{2,Float64}(0.0, SVector(0.0, 1.0), 0..0)


x = Aff(3..5, 2, 1)
y = Aff(2..4, 2, 2)
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
# z = Aff(-1..1, 1, 1)
# z^2
# interval(z^2)
#
# using Polynomials
#
# p = Poly([-3, 1])
# p2 = p^8
#
# x = 4 ± 1e-4
# y = Aff(x, 1, 1)
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
# xx = Aff(X[1], 2, 1)
# yy = Aff(X[2], 2, 2)
#
# interval.(f((xx, yy)))
#
# f(X)
#
#
#
#
# x = Aff(4..6, 1, 1)    # example from Messine
# f(x) = x * (10 - x)
#
# f(x)
# interval(f(x))
#
# interval(10*x - x^2)

"General formula for affine approximation of nonlinear functions"
function affine_approx(x::Aff, α, ζ, δ)

    c = α * x.c + ζ
    γ = α .* x.γ
    δ += α * x.Δ  # interval

    return Aff(c, γ, δ)
end

function Base.sqrt(x::Aff, X=interval(x))

    a, b = X.lo, X.hi

    # @show a, b

    # min-range:  de Figuereido book, pg. 64
    α = 1 / (2*√b)
    ζ = (√b) / 2
    δ = ( (√(b) - √(a))^2 / (2*√b) ) * (-1..1)

    return affine_approx(x, α, ζ, δ)
end

function Base.inv(x::Aff, X=interval(x))

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


function Base.exp(x::Aff, X=interval(x))

    # min-range approximation:  de Figuereido book, pg. 69

    a, b = X.lo, X.hi

    # @show a, b

    exp_X = exp(X)
    exp_a, exp_b = exp_X.lo, exp_X.hi

    α = exp_a

    if α == 0
        d_min = zero(x.lo)
        d_max = exp_b

    else
        d_max = IntervalArithmetic.@round_up(exp_b - α*b)
        d_min = IntervalArithmetic.@round_down(exp_a - α*a)
    end

    d = interval(d_min, d_max)

    ζ = mid(d)
    δ = radius(d) * (-1..1)

    return affine_approx(x, α, ζ, δ)
end
