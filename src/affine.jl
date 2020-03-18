using StaticArrays
using IntervalArithmetic
import IntervalArithmetic: interval

using LinearAlgebra

import Base: +, *, ^, -, sqrt, inv, exp

"""
Affine form with center `c`, affine components `γ` and error `Δ`.

Variant where Δ is an interval
"""
struct Affine{N,T<:AbstractFloat}
    affine::Aff{N,T}
    range::Interval{T}
end

interval(x::Affine) = x.range
range(x::Affine) = x.range

eltype(::Affine{N, T}) where {N, T} = T
zero(::Affine{N, T}) where {N, T} = Affine(Interval(zero(T)))
zero(::Type{Affine{T}}) where {N, T} = Affine(Interval(zero(T)))

one(::Affine{T}) where T = Affine(Interval(one(T)))
one(::Type{Affine{T}}) where T = Affine(Interval(one(T)))

function Base.show(io::IO, C::Affine{N,T}) where {N,T}
    print(io, "affine=", C.affine, "; range=", C.range)
end

# ==(C::Affine, D::Affine) = C.c == D.c && C.γ == D.γ

"""
Make an `Affine` based on an interval, which is number `i` of `n` total variables.
"""
function Affine(X::Interval{T}, n, i) where {T}
    return Affine(Aff(X, n, i), X)
end

Affine(X::Interval) = Affine(X, 1, 1)
Affine(X::Number) = Affine(Interval(X), 1, 1)

affine(Xs::Interval...) = Affine.(Xs, length(Xs), 1:length(Xs))

for op in (:+, :*, :-)
    @eval function $op(x::Affine, y::Affine)
        affine = $op(x.affine, y.affine)

        range = $op(x.range, y.range)

        range = range ∩ interval(affine)

        return Affine(affine, range)
    end
end

for op in (:sqrt, :inv, :exp)
    @eval function $op(x::Affine)
        affine = $op(x.affine, x.range)

        range = $op(x.range)

        range = range ∩ interval(affine)

        return Affine(affine, range)
    end
end

*(x::Affine, α::Real) = Affine(α*x.affine, α*x.range)
*(α::Real, x::Affine) = x * α

+(x::Affine, α::Real) = Affine(α+x.affine, α + x.range)
+(α::Real, x::Affine) = x + α

-(x::Affine) = Affine(-x.affine, -x.range)
-(x::Affine, α::Real) = Affine(x.affine - α, x.range - α)
-(α::Real, x::Affine) = α + (-x)

/(x::Affine, α::Real) = Affine(x.affine/α, x.range/α)

function ^(x::Affine, n::Integer)

    invert = false

    if n < 0
        invert = true
        n = -n
        # @show n
    end

    result = Base.power_by_squaring(x, n)

    if invert
        result = inv(result)
    end

    return result
end
Base.literal_pow(::typeof(^), x::Affine, ::Val{p}) where {T,p} = x^p
