using ADerrors

# Input of uwreal's
a = uwreal(rand(1000), "White noise")
b = uwreal([1.0, 0.1], "Var with error")
p = cobs([1.0, 2.0], [1.0 0.1;
                      0.1 2.0], "Parameters")

# Most common operations.
# You might add something else if
# you use it frequently
for op in (:-, :sin, :cos, :log, :log10, :log2, :sqrt, :exp, :exp2, :exp10, :sinh, :cosh, :tanh)
    @eval c = $op(a)
end
c = 1.0 + b
c = 1.0 - b
c = 1.0 * b
c = 1.0 / b
c = a + 2.0
c = a - 2.0
c = a * 2.0
c = a / 2.0
c = a ^ 3
c = a + b
c = a - b
c = a * b
c = a / b

# Error analysis
uwerr(c)
cov([c, a, b])
trcov([1.0 2.0 3.0;
       2.0 1.0 0.4;
       3.0 0.4 1.0], [c, a, b])
details(c)

# Error analysis of fit parameters
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
y = [0.0802273592699947
  0.09150837606934502
  0.047923648239388834
  0.024851583326401416
  0.01635482325054799
  0.13115281737744588
  0.21013177679604178
  0.002143355151617357
  0.2292183950698425
 -0.05174734852593241
  0.1384913891139784
 -0.05211234898997283]
dt = cobs(y, sig, "Fit data")
chisq(p, d) = sum( (d .- p[1]) .^ 2 ./ dx .^2 )
xp = [sum(value.(dt) ./ dx)/sum(1.0 ./ dx)]
(fitp, csqexp) = fit_error(chisq, xp, dt)
chiexp(chisq, xp, dt)

