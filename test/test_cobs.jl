using ADerrors # hide
# Put some average values and covariance 
avg = [16.26, 0.12, -0.0038]
Mcov = [0.478071 -0.176116 0.0135305
        -0.176116 0.0696489 -0.00554431
        0.0135305 -0.00554431 0.000454180]

# Produce observables with ensemble ID 
# [1, 2001, 32]. Do error analysis
p = cobs(avg, Mcov, "GF beta function parameters")
uwerr.(p)

# Check central values are ok
avg2 = value.(p)
println("Better be zero: ", sum((avg.-avg2).^2))

# Check that the covariance is ok
Mcov2 = cov(p)
println("Better be zero: ", sum((Mcov.-Mcov2).^2))

(sum((Mcov.-Mcov2).^2) < 1.0E-10)
