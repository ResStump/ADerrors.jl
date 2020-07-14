# ADerrors.jl

This package implementes error analysis of Monte Carlo data with the (``\Gamma``) method and automatic differentiation for error propagation.
```@contents
Pages = ["index.md"]
Depth = 3
```

## Getting started

`ADerrors.jl` is a package for error propagation and analysis of Monte carlo data. At the core of the package is the `uwreal` data type, that is able to store variables with uncertainties
```@repl gs
using ADerrors
a = uwreal([1.0, 0.1], 1) # 1.0 +/- 0.1
```
It can also store MC data
```@repl gs
# Generate some correlated data
eta  = randn(1000);
x    = Vector{Float64}(undef, 1000);
x[1] = 0.0;
for i in 2:1000
    x[i] = x[i-1] + eta[i]
    if abs(x[i]) > 1.0
        x[i] = x[i-1]
    end
end

b = uwreal(x.^2, 200)
c = uwreal(x.^4, 200)
```
Correlations between variables are treated consistently. This requires that each variable that is defined with `uwreal` contains an ensemble `ID` tag. In the previous examples `a` has been measured on ensemble `ID` 1, while both `b` and `c` have been measured on ensemble `ID` 200. This will treat the measurements in `b` and `c` as statistically correlated, while `a` will be uncorrelated with both `b` and `c`.

One can perform operations with `uwreal` variables as if there were normal floats
```@repl gs
d = 2.0 + sin(b)/a + 2.0*log(c)/(a+b+1.0)
```
Error propagation is done automatically, and correlations are consistently taken into account. Note that once complex derived `uwreal` variables are defined, as for example `d` above, they will in general depend on more than one ensemble `ID`. 

In order to perform the error analysis of a variable, one should use the `uwerr` function
```@repl gs
uwerr(d);
println("d: ", d)
```
One can get detailed information on the error of a variables 
```@repl gs
uwerr(d);
println("Details on variable d: ")
details(d)
```
where we can clearly see that there are two ensembles contributing to the uncertainty in `d`. We recognize this as `d` being a function of both `a` (measured on ensemble 1), and `b` and `c`, measured on ensemble 200.

Note that one does not need to use `uwerr` on a variable unless one is interested in the error on that variable. For example
```@repl gs
global x = 1.0
for i in 1:100
	global x = x - (d*cos(x) - x)/(-d*sin(x) - 1.0)
end
uwerr(x)
print("Root of d*cos(x) - x:" )
details(x)
```
determines the root on ``d\cos(x) - x`` by using Newton's method, and propagates the error on `d` into an error on the root. The error on `x` is only determined once, after the 100 iterations are over. 

### Advanced topics

#### Errors in fit parameters

`ADerrors.jl` does not provide an interface to perform fits, but once the minima of the ``\chi^2`` has been found, it can help propagating the error from the data to the fit parameters. 

Here we are going to [repeat the example](https://github.com/JuliaNLSolvers/LsqFit.jl) of the package `LsqFit.jl` using `ADerrors.jl` for error propagation.
```@repl fits
using LsqFit, ADerrors

# a two-parameter exponential model
# x: array of independent variables
# p: array of model parameters
# model(x, p) will accept the full data set as the first argument `x`.
# This means that we need to write our model function so it applies
# the model to the full dataset. We use `@.` to apply the calculations
# across all rows.
@. model(x, p) = p[1]*exp(-x*p[2])

# some example data
# xdata: independent variables
# ydata: dependent variable as uwreal
xdata = range(0, stop=10, length=20);
ydata = Vector{uwreal}(undef, length(xdata));
for i in eachindex(ydata)
	ydata[i] = uwreal([model(xdata[i], [1.0 2.0]) + 0.01*getindex(randn(1),1), 0.01], i)
	uwerr(ydata[i]) # We will need the errors for the weights of the fit
end
p0 = [0.5, 0.5];
```
We are ready to fit `(xdata, ydata)` to our model using `LsqFit.jl`, but for error propagation using `ADerrors.jl` we need the ``\chi^2`` function.
```@repl fits
# ADerrors will need the chi^2 as a function of the fit parameters and 
# the data. This can be constructed easily from the model above
chisq(p, d) = sum( (d .- model(xdata, p)) .^2 ./ 0.01^2)
```
Now we can fit the data and compute the uncertainties in our fit parameters
```@repl fits
fit = curve_fit(model, xdata, value.(ydata), 1.0 ./ err.(ydata).^2, p0) # This is LsqFit.jl
(fitp, cexp) = fit_error(chisq, coef(fit), ydata); # This is error propagation with ADerrors.jl
uwerr.(fitp);

println("chi^2 / chi^2_exp: ", sum(fit.resid .^2), " / ", cexp, " (dof: ", dof(fit), ")")
for i in 1:2
	println("Fit parameter: ", i, ": ", fitp[i])
end
```

!!! note 
    `ADerrors.jl` is completely agnostic about fitting the data. Error propagation is performed once the minimum of the ``\chi^2`` function is known, and it does not matter how this minima is found. One is not constrained to use `LsqFit.jl`, as there are many alternatives in the `Julia` ecosystem: `LeastSquaresOptim.jl`, `MINPACK.jl`, `Optim.jl`, ...


#### Missing measurements in one ensemble

`ADerrors.jl` can deal with observables that are not measured on every configuration, and still take correctly the correlations/autocorrelations into account. Here we show an extreme case where one observable is measured only in the even configurations, while the oher is measured on the odd coficurations.
```@repl gaps
using ADerrors, Plots
pgfplotsx();

# Generate some correlated data
eta  = randn(10000);
x    = Vector{Float64}(undef, 10000);
x[1] = 0.0;
for i in 2:10000
    x[i] = x[i-1] + 0.2*eta[i]
    if abs(x[i]) > 1.0
        x[i] = x[i-1]
    end
end

# x^2 only measured on odd configurations
x2 = uwreal(x[1:2:9999].^2,  1001, collect(1:2:9999),  10000)
# x^4 only measured on even configurations
x4 = uwreal(x[2:2:10000].^4, 1001, collect(2:2:10000), 10000)

rat = x2/x4
uwerr(rat);
```
In this case `uwerr` complains because with the chosen window the variance for ensemble with ID 1001 is negative. Looking at the normalized autocorrelation function we can easily spot the problem
```@repl gaps
iw = window(rat, 1001)
r  = rho(rat, 1001);
dr = drho(rat, 1001);
plot(collect(1:100), 
	r[1:100], 
	yerr = dr[1:100], 
	seriestype = :scatter, title = "Chosen Window: " * string(iw))
savefig("rat_cf.png") # hide
```
![rat plot](rat_cf.png)
The normalized autocorrelation function is oscillating, and the chosen window is 2! We better fix the window to 50 for this case
```@repl gaps
wpm = Dict{Int64, Vector{Float64}}()
wpm[1001] = [50.0, -1.0, -1.0, -1.0]
uwerr(rat, wpm)
println("Ratio: ", rat)
```

Note however that this is very observable dependent
```@repl gaps
prod = x2*x4
uwerr(prod)
iw = window(prod, 1001)
r  = rho(prod, 1001);
dr = drho(prod, 1001);
plot(collect(1:2*iw), 
	r[1:2*iw], 
	yerr = dr[1:2*iw], 
	seriestype = :scatter, title = "Chosen Window: " * string(iw))
savefig("prod_cf.png") # hide
```
![pod plot](prod_cf.png)
In general error analysis with data with arbitrary gaps is possible, and fully supported in `ADerrors.jl`, but it can be tricky, and certaiinly requires to examine the data carefully.

## Creating `uwerr` data types

```@docs
uwreal
cobs
```

## I/O

```@docs
write_uwreal
read_uwreal
```

## Error analysis 

```@docs
uwerr
cov
trcov
```

## Information on error analysis

There are several methods to get information about the error analysis of a `uwreal` variable. With the exception of `value`, these methods require that a proper error analysis has been performed via a call to the `uwerr` method.

```@docs
value
err
derror
taui
dtaui
rho
drho
window
details
neid
```

## Error propagation in iterative algorithms

### Finding root of functions

```@docs
root_error
```

### Fits

`ADerrors.jl` is agnostic about how you approach a fit (i.e. minimize a ``\chi^2`` function), but once the central values of the fit parameters have been found (i.e. using `LsqFit.jl`, `LeastSquaresOptim.jl`, etc... ) `ADerrors.jl` is able to determine the error of the fit parameters, returning them as `uwreal` type. It can also determine the expected value of your ``\chi^2`` given the correlation in your fit parameters.

```@docs
fit_error
chiexp
```

### Integrals

```@docs
int_error
```

## Index 

```@index
```
