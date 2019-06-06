using StaticArrays
using IntervalArithmetic
import IntervalArithmetic: interval

using LinearAlgebra

import Base: +, *, ^, -, sqrt, inv

"""
Affine form with center `c`, affine components `γ` and error `Δ`.

Variant where Δ is an interval
"""
struct AFF{N,T<:AbstractFloat}
    affine::AF{N,T}
    range::Interval{T}
end


function Base.show(io::IO, C::AFF{N,T}) where {N,T}
    print(io, "affine=", C.affine, "; range=", C.range)
end

# ==(C::Affine, D::Affine) = C.c == D.c && C.γ == D.γ

"""
Make an `AF` based on an interval, which is number `i` of `n` total variables.
"""
function AFF(X::Interval{T}, n, i) where {T}
    return AFF(AF(X, n, i), X)
end

for op in (:+, :*, :-)
    @eval function $op(x::AFF, y::AFF)
        affine = $op(x.affine, y.affine)

        range = $op(x.range, y.range)

        range = range ∩ interval(affine)

        return AFF(affine, range)
    end
end

for op in (:sqrt, :inv)
    @eval function $op(x::AFF)
        affine = $op(x.affine)

        range = $op(x.range)

        range = range ∩ interval(affine)

        return AFF(affine, range)
    end
end



#
#
# *(x::AF, α::Real) = AF(α*x.c, α.*x.γ, α*x.Δ)
# *(α::Real, x::AF) = x * α
#
# +(x::AF, α::Real) = AF(α+x.c, x.γ, x.Δ)
# +(α::Real, x::AF) = x + α
#
# -(x::AF) = AF(-x.c, .-(x.γ), -x.Δ)
# -(x::AF, α::Real) = AF(x.c - α, x.γ, x.Δ)
# -(α::Real, x::AF) = α + (-x)
#
^(x::AFF, n::Integer) = Base.power_by_squaring(x, n)

Base.literal_pow(::typeof(^), x::AFF, ::Val{p}) where {T,p} = x^p

x = AF{2,Float64}(0.0, SVector(1.0, 0.0), 0..0)
y = AF{2,Float64}(0.0, SVector(0.0, 1.0), 0..0)

#
# x = AF(3..5, 2, 1)
# y = AF(2..4, 2, 2)
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
