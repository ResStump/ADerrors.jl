# ADerrors.jl

Error propagation and analysis of Monte Carlo data with the $`\Gamma`$ method and automatic differentiation in `Julia`

- [Features](#features)
- [Installation](#installation)
- [Tutorial](#tutorial)
- [Full documentation](#full-documentation)
- [Performance tips](#performance-tips)
- [How to cite](#how-to-cite)

## Features

- **Exact** linear error propagation, even in iterative algorithms.
- **Exact** linear error propagation in **fit parameters**,
  **integrals** and **roots** of non linear functions. 
- Handles data from **any number of ensembles** (i.e. simulations with
  different parameters).
- Support for **replicas** (i.e. several runs with the same simulation
  parameters). 
- **Irrelgular MC measurements** are handled transparently.

## Installation

The package in not in the general registry. Still one can use the
package manager. `ADerrors.jl` also depends on `BDIO.jl` that is also
not registered and should be installed beforehand:
```julia
julia> import Pkg
(v1.1) pkg> add https://gitlab.ift.uam-csic.es/alberto/bdio.jl
(v1.1) pkg> add https://gitlab.ift.uam-csic.es/alberto/aderrors.jl
```
## Tutorial

Please, have a look at the [Getting started](https://ific.uv.es/~alramos/docs/ADerrors/tutorial/) guide.

## Full documentation

The full documentation of the package is available via the usual
[Julia `REPL` help
mode](https://docs.julialang.org/en/v1/stdlib/REPL/#Help-mode-1) and
online in [HTML format](https://ific.uv.es/~alramos/docs/ADerrors/).

## Performance tips

`Julia` is well known for being slow the first time that you call a
function. This is because `Julia` not only runs the 
code, but also compiles it, making the first call slow.
This problem can be alleviated in general with
[PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl). This
is specially true for the case of `ADerrors` since most functions are
type static. 

As example, the file `extra/typical.jl` contains the most typical
calls to `ADerrors`. One can execute this file telling julia to
annotate the functions that are compiled
```julia
julia --trace-compile=precompile_aderrors.jl typical.jl
```
Now the functions annotated in `precompile_aderrors.jl` can be
compiled and included in a `sysimage` that is automatically loaded
whenever you start Julia
```julia
julia> using PackageCompiler
julia> PackageCompiler.create_sysimage(:ADerrors; precompile_statements_file="precompile_aderrors.jl", replace_default=true)
```

This will make `ADerrors` fast from the first call. Obviously you can tune
the file `typical.jl` to your usage, or add other packages. Please
note that packages included in the sysimage are locked to the versions
of the sysimage. If you update `ADerrors` make sure to re-generate the
sysimage. Probably is better to read [the documentation of
PackageCompiler](https://julialang.github.io/PackageCompiler.jl/dev/sysimages/)
in order to fully understand the drawbacks.

## How to cite

This work is an implementation of several ideas in data analysis. If you use this package for your scientific work, please consider citing:
- U. Wolff, "Monte Carlo errors with less errors".
  Comput.Phys.Commun. 156 (2004) 143-153. 
- F. Virotta, "Critical slowing down and error analysis of lattice QCD simulations." PhD thesis.
- Stefan Schaefer, Rainer Sommer, Francesco Virotta, "Critical slowing
  down and error analysis in lattice QCD simulations". Nucl.Phys.B 845 (2011) 93-119.
- A. Ramos, "Automatic differentiation for error analysis of Monte Carlo data". Comput.Phys.Commun. 238 (2019) 19-35. 
- M. Bruno, R. Sommer, In preparation.



