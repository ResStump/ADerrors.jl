using ADerrors, Distributions # hide

# Generate correlated samples with average 0.1
npt = 12
sig = zeros(npt, npt)
dx  = zeros(npt)
for i in 1:npt
    dx[i]    = 0.01*i
    sig[i,i] = dx[i]^2
    for j in i+1:npt
        sig[i,j] = 0.0001 - 0.000005*abs(i-j)
        sig[j,i] = 0.0001 - 0.000005*abs(i-j)
    end
end
dmv = MvNormal([0.1 for n in 1:npt], sig)
vs  = rand(dmv, 1)

# Create the uwreal data that we want to 
# fit to a constant
dt = cobs(vs[:,1], sig, "Data points")

# Define the chi^2
chisq(p, d) = sum( (d .- p[1]) .^ 2 ./ dx .^2 )

# The result of an uncorrelated fit to a 
# constant is the weighted average
xp = [sum(value.(dt) ./ dx)/sum(1.0 ./ dx)]


# Propagate errors to the fit parameters and
# determine the expected chi^2
(fitp, csqexp) = fit_error(chisq, xp, dt)

uwerr.(fitp)
println(" *** FIT RESULTS ***")
print("Fit parameter:     ")
details.(fitp)
println("chi^2 / chi_exp^2: ", chisq(xp, value.(dt)), " / ", csqexp, "  (dof: ", npt-1, ")")


(0==0)
