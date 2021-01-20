# Getting started

## Basic tutorial
`ADerrors.jl` is a package for error propagation and analysis of Monte
carlo data. At the core of the package is the `uwreal` data type. It
contains variables with uncertainties. In very simple cases this is
just a central value and an error
```@repl gs
using ADerrors
a = uwreal([1.0, 0.1], "Var with error") # 1.0 +/- 0.1
```
But it also supports the case where the uncertainty comes from Monte
Carlo measurements
```@repl gs
# Generate some MC correlated data
eta  = randn(1000);
x    = Vector{Float64}(undef, 1000);
x[1] = 0.0;
for i in 2:1000
    x[i] = x[i-1] + eta[i]
    if abs(x[i]) > 1.0
        x[i] = x[i-1]
    end
end

b = uwreal(x.^2, "Random walk ensemble in [-1,1]")
c = uwreal(x.^4, "Random walk ensemble in [-1,1]")
```
Correlations between variables are treated consistently. Each call to `uwreal` contains an ensemble `ID` tag. In the previous examples `a` has been measured on ensemble `ID` `Var with error`, while both `b` and `c` have been measured on ensemble `ID` `Random walk ensemble in [-1,1]`. Monte Carlo measurements with the same `ID` tag are assumed to be measured on the same configurations: `b` and `c` above are statistically correlated, while `a` will be uncorrelated with both `b` and `c`. 

One can perform operations with `uwreal` variables as if there were normal floats
```@repl gs
d = 2.0 + sin(b)/a + 2.0*log(c)/(a+b+1.0)
```
Error propagation is done automatically, and correlations are consistently taken into account. Note that once complex derived `uwreal` variables are defined, as for example `d` above, they will in general depend on more than one ensemble `ID`. 

The message "(Error not available... maybe run uwerr)" in all the above statements suggests that newly defined `uwreal` variables do not know their error yet. In order to perform the error analysis of a variable, one should use
the `uwerr` function
```@repl gs
uwerr(d);
println("d: ", d)
```
A call to `uwerr` will apply the ``\Gamma``-method to MC ensembles. `ADerrors.jl` can give detailed information on the error of a variable (this requires that `uwerr` has been called)
```@repl gs
println("Details on variable d: ")
details(d)
```
Here we can see that there are two ensembles contributing to the uncertainty in `d`. We recognize this as `d` being a function of both `a` (measured on ensemble `Var with error`), and `b` and `c`, measured on ensemble `Random walk ensemble in [-1,1]`. Most of the error in `d` comes in fact from ensemble `Random walk ensemble in [-1,1]`.

Note that one does not need to use `uwerr` on a variable unless one is interested in the error on that variable. For example, continuing with the previous example
```@repl gs
global x = 1.0
for i in 1:100
	global x = x - (d*cos(x) - x)/(-d*sin(x) - 1.0)
end
uwerr(x)
print("Root of d*cos(x) - x:" )
details(x)
```
determines the root on ``d\cos(x) - x`` by using Newton's method, and propagates the error on `d` into an error on the root. The error on `x` is only determined once, after the 100 iterations are over. This feature is useful becase calls to `uwerr` are typically numerically expensive in comparison with mathematical operations of `uwreal`'s.

## Advanced topics

### Errors in fit parameters

`ADerrors.jl` does not provide an interface to perform fits, but once the minima of the ``\chi^2`` has been found, it can return the fit parameters as `uwreal` variables.

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
	ydata[i] = uwreal([model(xdata[i], [1.0 2.0]) + 0.01*getindex(randn(1),1), 0.01], "Point "*string(i))
	uwerr(ydata[i]) # We will need the errors for the weights of the fit
end
p0 = [0.5, 0.5];
```
We are ready to fit `(xdata, ydata)` to our model using `LsqFit.jl`
```@repl fits
# Fit the data with LsqFit.jl
fit = curve_fit(model, xdata, value.(ydata), 1.0 ./ err.(ydata).^2, p0) 
```
In order to use `ADerrors.jl` we need the ``\chi^2`` function.
```@repl fits
# ADerrors will need the chi^2 as a function of the fit parameters and 
# the data. This can be constructed easily from the model function above
chisq(p, d) = sum( (d .- model(xdata, p)) .^2 ./ 0.01^2)
```
A call to `fit_error` will return the fit parameters as a `Vector{uwreal}` and also the expected value of the ``\chi^2`` as a `Float64`
```@repl fits
# This is error propagation with ADerrors.jl
(fitp, cexp) = fit_error(chisq, coef(fit), ydata);
println("Fit results: \n
	     chi^2 / chi^2_exp: ", sum(fit.resid .^2), " / ", cexp, " (dof: ", dof(fit), ")")

# fitp are the fit parameters in uwreal
# format. One can use them as any other uwreal
# For example to determine the main contribution
# to its uncertainty
for i in 1:2
	uwerr(fitp[i])
	print("Fit parameter: ", i, ": ")
	details(fitp[i])
end
```

!!! note 
    `ADerrors.jl` is completely agnostic about fitting the data. Error propagation is performed once the minimum of the ``\chi^2`` function is known, and it does not matter how this minima is found. One is not constrained to use `LsqFit.jl`, as there are many alternatives in the `Julia` ecosystem: `LeastSquaresOptim.jl`, `MINPACK.jl`, `Optim.jl`, ...


### Missing measurements 

`ADerrors.jl` can deal with observables that are not measured on every configuration, and still take correctly the correlations/autocorrelations into account. Here we show an extreme case where one observable is measured only in the even configurations, while the other is measured on the odd configurations.
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
x2 = uwreal(x[1:2:9999].^2,  "Random walk in [-1,1]", collect(1:2:9999),  10000)
# x^4 only measured on even configurations
x4 = uwreal(x[2:2:10000].^4, "Random walk in [-1,1]", collect(2:2:10000), 10000)
```
the variables `x2` and `x4` are normal `uwreal` variables, and one can work with them as with any other variable. Note however that the automatic determination of the summation window can fail in these cases, because the autocorrelation function with missing measurements can be very complicated. For example
```@repl gaps
rat = x2/x4 # This is perfectly fine
try # Use automatic window: might fail!
    uwerr(rat)
catch e
    println("Automatic window fails")
end
```
In this case `uwerr` complains because with the automatically chosen window the variance for ensemble with ID `"Random walk in [-1,1]"` is negative. Let's have a look at the normalized autocorrelation function
```@repl gaps
iw = window(rat, "Random walk in [-1,1]")
r  = rho(rat, "Random walk in [-1,1]");
dr = drho(rat, "Random walk in [-1,1]");
plot(collect(1:100), 
	r[1:100], 
	yerr = dr[1:100], 
	seriestype = :scatter, title = "Chosen Window: " * string(iw), label="autoCF")
savefig("rat_cf.png") # hide
```
![rat plot](rat_cf.png)
The normalized autocorrelation function is oscillating, and the automatically chosen window is 2, producing a negative variance! We better fix the window to 50 for this case
```@repl gaps
wpm = Dict{String, Vector{Float64}}()
wpm["Random walk in [-1,1]"] = [50.0, -1.0, -1.0, -1.0]
uwerr(rat, wpm) # repeat error analysis with our choice (window=50)
println("Ratio: ", rat)
```
After a visual examination of the autocorrelation function and a reasonable choice of the summation window, the error in `rat` is determined. 
Note however that it is very difficult to anticipate which observables will have a strange autocorrelation function. For example the product of `x2` and `x4` shows a reasonably well behaved autocorrelation function.
```@repl gaps
prod = x2*x4
uwerr(prod)
println("Product: ", prod)
iw = window(prod, "Random walk in [-1,1]")
r  = rho(prod, "Random walk in [-1,1]");
dr = drho(prod, "Random walk in [-1,1]");
plot(collect(1:2*iw), 
	r[1:2*iw], 
	yerr = dr[1:2*iw], 
	seriestype = :scatter, title = "Chosen Window: " * string(iw), label="autoCF")
savefig("prod_cf.png") # hide
```
![pod plot](prod_cf.png)
In general error analysis with data with arbitrary gaps is possible, and fully supported in `ADerrors.jl`. However in cases where many gaps are present, data analysis can be tricky, and certainly requires some care.

### `BDIO` Native format

`ADerrors.jl` can save/read observables (i.e. `uwreal`'s) in `BDIO` records. Here we describe the format of these records.

All observables are stored in a single `BDIO` record of type `BDIO_BIN_GENERIC`. The record has the following data in binary little endian format
- mean (`Float64`): The mean value of the observable.
- neid (`Int32`): The number of ensembles contributing to the observable.
- ndata (`Vector{Int32}(neid)`): The length of each ensemble.
- nrep (`Vector{Int32}(neid)`): The number of replica of each enemble.
- vrep (`Vector{Int32}`): The replica vector for each ensemble.
- ids (`Vector{Int32}(neid)`): The ensemble numeric `ID`'s.
- nt (`Vector{Int32}(neid)`): Half the largest replica of each ensemble.
- zero (`Vector{Float64}(neid)`): just `neid` zeros.
- four (`Vector{Float64}(neid)`): just `neid` fours.
- delta  (`Vector{Float64}`): The fluctuations for each ensemble.
- name (NULL terminated `String`): A description of the observable.
- `ID` tags: A list of `neid` tuples `(Int32, String)` that maps each numeric `ID` to an ensemble tag. All strings are NULL terminated.
- Replica names: A list of `neid` tuples `(Int32, Vector{String})`
  that maps each numeric `ID` to a vector of replica names. All
  strings are NULL terminated. 

!!! alert
    Obviously this weird format is what it is for some legacy reasons, but it is strongly encouraged that new implementations respect this standard with all its weirdness. 




