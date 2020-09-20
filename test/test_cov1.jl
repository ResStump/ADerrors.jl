using ADerrors

# Average values and covariance from
# https://inspirehep.net/literature/1477411
# here we try to reproduce the result
# Scale ratio = 21.86(42) Eq.5.6
avg = [16.26, 0.12, -0.0038]
Mcov = [0.478071 -0.176116 0.0135305
        -0.176116 0.0696489 -0.00554431
        0.0135305 -0.00554431 0.000454180]

p = cobs(avg, Mcov, "Beta function fit parameters")
g1s = uwreal([2.6723, 0.0064], 4)
g2s = 11.31

fint(x, p) = - (p[1] + p[2]*x^2 + p[3]*x^4)/x^3
g1 = sqrt(g1s)
g2 = sqrt(g2s)
sint = 2.0*exp(-int_error(fint, g1, g2, p))
uwerr(sint)
print("  From integral evaluation: ")
details(sint)
