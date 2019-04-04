# Wrapping effect example from Rump paper (2015)

using AffineArithmetic
using ValidatedNumerics
using LinearAlgebra

A = 0.5 * [1 2; -1 1]
eigvals(A)
maximum(abs.(eigvals(A)))  # spectral radius

reset_affine_index()
X = Affine(-1..1)
Y = Affine(-1..1)

range(X)

XX = [X, Y]  # vector

for i in 1:10
    XX = A * XX
    println(i)
    println(XX)
    println(range.(XX))
    println()
end

vertices = []
using Iterators
for (ε1, ε2, ε3, ε4) in product( ([-1, 1] for i in 1:4)... )
    @show (ε1, ε2, ε3, ε4)
    x = 2 + ε1 - 2ε2 + 3ε3 - ε4
    y = 1 + 3ε1 - ε3 + 2ε4

    @show x, y

    push!(vertices, [x, y])
end

vertices
xs = [p[1] for p in vertices]
ys = [p[2] for p in vertices]

using Plots
gr()

scatter(xs, ys, aspect_ratio=:equal)
# Need the convex hull of these points for the zonotope
