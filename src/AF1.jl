using StaticArrays
using IntervalArithmetic
import IntervalArithmetic: interval

import Base: +, *, ^, -

"""
Affine form with center `c`, affine components `γ` and error `Δ`.
"""
struct AF1{N,T<:AbstractFloat}
    c::T   # mid-point
    γ::SVector{N,T}  # affine terms
    Δ::T   # error term.  Error is Δ.(-1..1)
end


function Base.show(io::IO, C::AF1{N,T}) where {N,T}
    print(io, "⟨", C.c, "; ", C.γ, "; ", C.Δ, "⟩")
end

# ==(C::Affine, D::Affine) = C.c == D.c && C.γ == D.γ

"""
Make an `AF1` based on an interval, which is number `i` of `n` total variables.
"""
function AF1(X::Interval{T}, n, i) where {T}
    c = mid(X)
    r = radius(X)

    γ = SVector(ntuple(j->i==j ? r : zero(r), n))

    return AF1(c, γ, T(0))
end

+(x::AF1{N,T}, y::AF1{N,T}) where {N,T} = AF1(x.c + y.c, x.γ .+ y.γ, x.Δ + y.Δ)

-(x::AF1{N,T}, y::AF1{N,T}) where {N,T} = AF1(x.c - y.c, x.γ .- y.γ, x.Δ - y.Δ)


interval(C::AF1) = C.c + sum(abs.(C.γ))*(-1..1) + (C.Δ) * (-1..1)


function *(x::AF1{N,T}, y::AF1{N,T}) where {N,T}
    c = x.c * y.c

    γ = x.c .* y.γ + y.c .* x.γ

    Δ = abs(x.c) * y.Δ + abs(y.c) * x.Δ
    Δ += ( sum(abs.(x.γ)) + abs(x.Δ) ) * ( sum(abs.(y.γ)) + abs(y.Δ) )

    return AF1(c, γ, Δ)

end

*(x::AF1, α::Real) = AF1(α*x.c, α.*x.γ, α*x.Δ)
*(α::Real, x::AF1) = x * α

+(x::AF1, α::Real) = AF1(α+x.c, x.γ, x.Δ)
+(α::Real, x::AF1) = x + α

-(x::AF1) = AF1(-x.c, .-(x.γ), -x.Δ)
-(x::AF1, α::Real) = AF1(x.c - α, x.γ, x.Δ)
-(α::Real, x::AF1) = α + (-x)

^(x::AF1, n::Integer) = Base.power_by_squaring(x, n)

Base.literal_pow(::typeof(^), x::AF1, ::Val{p}) where {T,p} = x^p




x = AF1(3..5, 2, 1)
y = AF1(2..4, 2, 2)

3-x
interval(3-x)

x + y


interval(x+y)

x * y
interval(x * y)

interval(x * y)
interval(x) * interval(y)

z = AF1(-1..1, 1, 1)
z^2
interval(z^2)

using Polynomials

p = Poly([-3, 1])
p2 = p^8

x = 4 ± 1e-4
y = AF1(x, 1, 1)

interval(y)
interval(p2(x))
interval(p2(y))

@time interval(p2(y))


f( (x, y) ) = [x^3 + y, (x - y)^2]

X = IntervalBox(-1..1, -1..1)

f(X)

xx = AF1(X[1], 2, 1)
yy = AF1(X[2], 2, 2)

interval.(f((xx, yy)))

f(X)




x = AF1(4..6, 1, 1)    # example from Messine
f(x) = x * (10 - x)

f(x)
interval(f(x))

interval(10*x - x^2)


function Base.inv(x::AF1)

    range = interval(x)
    a, b = range.lo, range.hi

    # assumes b > a > 0:
    # min-range approx:

    p = -1 / (b^2)
    q = -p * (a + b)^2 / (2a)
    δ = - p * (a- b)^2 / (2a)

    return AF1(p * x.c + q, p .* x.γ, x.Δ)


end
