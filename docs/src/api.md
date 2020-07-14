# API


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

