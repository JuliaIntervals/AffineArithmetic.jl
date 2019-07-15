module AffineArithmetic

using IntervalArithmetic

import Base: +, -, *, /, ^, ==,
            zero, one, range,
            show


include("aff.jl")
include("affine.jl")
include("full_affine.jl")


export Affine
export FullAffine, reset_affine_index

end
