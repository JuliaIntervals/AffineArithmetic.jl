using AffineArithmetic
using ValidatedNumerics

X = 1..3
X_aff = Affine(X)

Affine(0.0)

range(X_aff)

X_aff - X_aff

X_aff^2

range(X_aff^2)

X^2 - (X+X)
Y = X_aff^2 - (X_aff + X_aff)
range(Y)

o = Affine(1.0) # constant 1
X_aff - o
(X_aff - o)^2
(X_aff - o)^2 - o

X_aff^2 - (X_aff + X_aff)
(X_aff - o)^2 - o

X_aff^2 - (X_aff + X_aff) == (X_aff - o)^2 - o

range(X_aff^2 - (X_aff + X_aff))

X^2 - 2X
(X-1)^2 - 1   # true image, since only single use of interval

C = Affine(3.0, [4.0, 5.0])

range(C)

C = Affine(2., [1., -2, 3, 2])
D = Affine(1., [3., 0., -1.])

C + D

C - D

C = Affine(2., [2.])

C - C

D = Affine(1., [2.])

C = Affine(2., [1.])
C * C
range(C * C)
range(C)
range(C)^2

C * D
range(C+D)
range(C) + range(D)
