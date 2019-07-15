
const affine_index = [1]  # which affine vector index to use

reset_affine_index() = affine_index[1] = 1

@doc """
An affine quantity for affine arithmetic.
The usual way to create an affine quantity is from an interval `X`, via `FullAffine(X)`.
"""
struct FullAffine{T<:AbstractFloat}
    c::T   # mid-point
    γ::Vector{T}  # error terms
end


function show(io::IO, C::FullAffine)
    print(io, "⟨", C.c, "; ", C.γ, "⟩")
end

==(C::FullAffine, D::FullAffine) = C.c == D.c && C.γ == D.γ

@doc """
    FullAffine(X::Interval)

Construct a new `FullAffine` quantity from an `Interval`.
"""
function FullAffine(X::Interval)
    c = mid(X)
    r = radius(X)

    index = affine_index[1]
    affine_index[1] += 1

    γ = zeros(index)
    γ[end] = r

    return FullAffine(c, γ)
end

# conversion of numerical constant to affine:
FullAffine(c::Real) = FullAffine(c, Float64[])


range(C::FullAffine) = C.c + sum(abs.(C.γ))*(-1..1)

range(X::Interval) = X

# morally:
# +(C::FullAffine, D::Afine) = FullAffine(C.c + D.c, C.γ + D.γ)

for op in (:+, :-)
    @eval function $op(C::FullAffine, D::FullAffine)
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

        FullAffine($op(C.c, D.c), γ)
    end
end

+(C::FullAffine, α::Real) = FullAffine(C.c + α, C.γ)
+(α::Real, C::FullAffine) = C + α

-(C::FullAffine, α::Real) = FullAffine(C.c - α, C.γ)
-(α::Real, C::FullAffine) = FullAffine(α - C.c, [-x for x in C.γ])

function *(C::FullAffine, D::FullAffine)

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

    FullAffine(C.c * D.c, γ)
end

*(α::Real, C::FullAffine) = FullAffine(α*C.c, α*C.γ)
*(C::FullAffine, α::Real) = α * C

^(x::FullAffine, n::Integer) = Base.power_by_squaring(x, n)

Base.literal_pow(::typeof(^), x::FullAffine{T}, ::Val{p}) where {T,p} = x^p

eltype(C::FullAffine{T}) where T = T
zero(C::FullAffine{T}) where T = FullAffine(zero(T))
zero(::Type{FullAffine{T}}) where T = FullAffine(zero(T))

one(C::FullAffine{T}) where T = FullAffine(one(T))
one(::Type{FullAffine{T}}) where T = FullAffine(one(T))

#one(C::FullAffine)
