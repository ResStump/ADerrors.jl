using ADerrors, Random # hide
rng = Random.MersenneTwister(1234)
a = uwreal(rand(rng, 2000),   "Ensemble A12")
b = uwreal([1.2, 0.023], "Ensemble XYZ")
c = uwreal([5.2, 0.03],  "Ensemble RRR")
d = a + b - c
uwerr(d)

details(d)

(0 == 0)
