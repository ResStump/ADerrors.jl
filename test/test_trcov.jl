using ADerrors, LinearAlgebra # hide
a = uwreal([1.3, 0.01], 1) # 1.3 +/- 0.01
b = uwreal([5.3, 0.23], 2) # 5.3 +/- 0.23
c = uwreal(rand(2000), 3)

x = [a+b+sin(c), a-b+cos(c), c-b/a]
M = [1.0 0.2 0.1
     0.2 2.0 0.3
     0.1 0.3 1.0]

mcov = cov(x)
d = tr(mcov * M)
(abs(d -trcov(M, x)) < 1.0E-10)
