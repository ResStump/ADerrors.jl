using ADerrors, Distributions, LinearAlgebra # hide

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

(fitp, csqexp) = fit_error(chisq, xp, dt, C=sig)
chiexp = ADerrors.chiexp(chisq,xp,dt,C=sig)

uwerr.(fitp)
println(" *** UNCORRELATED FIT RESULTS ***")
print("Fit parameter:     ")
details.(fitp)
println("chi^2 / chi_exp^2: ", chisq(xp, value.(dt)), " / ", csqexp, "  (dof: ", npt-1, ")")

# The result of a correlated fit to a constant
W = LinearAlgebra.pinv(sig) |> x-> 0.5*(x+x')
xpc = [sum(W[i,j]*dt[j].mean for i in 1:npt, j in 1:npt)/sum(W)]

# Propagate errors to the fit parameters and
# determine the expected chi^2

(fitp, csqexp) = fit_error(chisq, xpc, dt,W=W,C=sig)
chiexp = ADerrors.chiexp(chisq,xp,dt,C=sig,W=W)

uwerr.(fitp)
println(" *** CORRELATED FIT RESULTS ***")
print("Fit parameter:     ")
details.(fitp)
println("chi^2 / chi_exp^2: ", chisq(xp, value.(dt)), " / ", csqexp, "  (dof: ", npt-1, ")")


(0==0)
