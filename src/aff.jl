using StaticArrays
using IntervalArithmetic
import IntervalArithmetic: interval

using LinearAlgebra

import Base: +, *, ^, -

"""
Affine form with center `c`, affine components `γ` and error `Δ`.

Variant where Δ is an interval
"""
struct Aff{N,T<:Interval}
    c::T   # mid-point
    γ::SVector{N,T}  # affine terms
    Δ::T   # error term
end


function Base.show(io::IO, C::Aff{N,T}) where {N,T}
    print(io, "⟨", C.c, "; ", C.γ, "; ", C.Δ, "⟩")
end

# ==(C::Affine, D::Affine) = C.c == D.c && C.γ == D.γ

"""
Make an `Aff` based on an interval, which is number `i` of `n` total variables.
"""
function Aff(X::Interval{T}, n::Int, i) where {T}
    c = interval(mid(X))
    r = interval(radius(X))

    γ = SVector(ntuple(j -> i==j ? r : zero(r), n))

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
