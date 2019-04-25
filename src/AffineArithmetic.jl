module AffineArithmetic

using IntervalArithmetic

import Base: +, -, *, /, ^, ==,
            zero, one, range,
            show


export Affine, reset_affine_index


const affine_index = [1]  # which affine vector index to use

reset_affine_index() = affine_index[1] = 1

@doc """
An affine quantity for affine arithmetic.
The usual way to create an affine quantity is from an interval `X`, via `Affine(X)`.
"""
struct Affine{T<:AbstractFloat}
    c::T   # mid-point
    γ::Vector{T}  # error terms
end


function show(io::IO, C::Affine)
    print(io, "⟨", C.c, "; ", C.γ, "⟩")
end

==(C::Affine, D::Affine) = C.c == D.c && C.γ == D.γ

@doc """
    Affine(X::Interval)

Construct a new `Affine` quantity from an `Interval`.
"""
function Affine(X::Interval)
    c = mid(X)
    r = radius(X)

    index = affine_index[1]
    affine_index[1] += 1

    γ = zeros(index)
    γ[end] = r

    return Affine(c, γ)
end

# conversion of numerical constant to affine:
Affine(c::Real) = Affine(c, Float64[])


range(C::Affine) = C.c + sum(abs.(C.γ))*(-1..1)

range(X::Interval) = X

# morally:
# +(C::Affine, D::Afine) = Affine(C.c + D.c, C.γ + D.γ)

for op in (:+, :-)
    @eval function $op(C::Affine, D::Affine)
        k = length(C.γ)
        l = length(D.γ)

        # account for unequal lengths:
        m = min(k, l)
        common_γ = $op(C.γ[1:m], D.γ[1:m])

        if l >= k
            γ = [ common_γ; D.γ[m+1:l] ]
        else
            γ = [ common_γ; C.γ[m+1:k] ]
        end

        Affine($op(C.c, D.c), γ)
    end
end

+(C::Affine, α::Real) = Affine(C.c + α, C.γ)
+(α::Real, C::Affine) = C + α

-(C::Affine, α::Real) = Affine(C.c - α, C.γ)
-(α::Real, C::Affine) = Affine(α - C.c, [-x for x in C.γ])

function *(C::Affine, D::Affine)

    c = C.c
    d = D.c

    k = length(C.γ)
    l = length(D.γ)

    # account for unequal lengths:
    m = min(k, l)
    common_γ = d * C.γ[1:m] + c * D.γ[1:m]

    error_bound = sum(abs.(C.γ)) * sum(abs.(D.γ))

    if l >= k
        γ = [ common_γ; c*D.γ[m+1:l]; error_bound ]
    else
        γ = [ common_γ; d*C.γ[m+1:k]; error_bound ]
    end

    # |γ|₁ = sum(abs.(γ))

    Affine(C.c * D.c, γ)
end

*(α::Real, C::Affine) = Affine(α*C.c, α*C.γ)
*(C::Affine, α::Real) = α * C

^(x::Affine, n::Integer) = Base.power_by_squaring(x, n)

Base.literal_pow(::typeof(^), x::Affine{T}, ::Val{p}) where {T,p} = x^p

eltype(C::Affine{T}) where T = T
zero(C::Affine{T}) where T = Affine(zero(T))
zero(::Type{Affine{T}}) where T = Affine(zero(T))

one(C::Affine{T}) where T = Affine(one(T))
one(::Type{Affine{T}}) where T = Affine(one(T))

#one(C::Affine)

end
