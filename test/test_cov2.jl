using ADerrors, LinearAlgebra # hide
a = uwreal([1.3, 0.01], "Var 1") # 1.3 +/- 0.01
b = uwreal([5.3, 0.23], "Var 2") # 5.3 +/- 0.23
uwerr(a)
uwerr(b)

x = [a+b, a-b]
mat = ADerrors.cov(x)
println("Covariance: ", mat[1,1], " ", mat[1,2])
println("            ", mat[2,1], " ", mat[2,2])
println("Check (should be zero): ",  mat[1,1] - mat[2,2])
println("Check (should be zero): ",  mat[1,2] - (err(a)^2-err(b)^2))
( abs(mat[1,1] - mat[2,2]) < 1.0E-10) && ( abs(mat[1,2] - (err(a)^2-err(b)^2)) < 1.0E-10 )

